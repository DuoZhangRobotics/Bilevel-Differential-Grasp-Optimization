#include "PhysicsRegistration.h"
#include "GravityEnergy.h"
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/MultiPrecisionLQP.h>
#include <Articulated/ArticulatedLoader.h>
#include <CommonFile/geom/BVHBuilder.h>
#include <Environment/Environment.h>
#include <Utils/RotationUtil.h>
#include <Utils/Utils.h>
#include <Eigen/Eigen>

USE_PRJ_NAMESPACE

//PhysicsRegistrationParameter
PhysicsRegistrationParameter::PhysicsRegistrationParameter(Options& ops)
{
  REGISTER_FLOAT_TYPE("coefPotential",PhysicsRegistrationParameter,scalarD,t._coefPotential)
  REGISTER_FLOAT_TYPE("coefRegister",PhysicsRegistrationParameter,scalarD,t._coefRegister)
  REGISTER_FLOAT_TYPE("tau",PhysicsRegistrationParameter,scalarD,t._tau)
  REGISTER_FLOAT_TYPE("g0",PhysicsRegistrationParameter,scalarD,t._g[0])
  REGISTER_FLOAT_TYPE("g1",PhysicsRegistrationParameter,scalarD,t._g[1])
  REGISTER_FLOAT_TYPE("g2",PhysicsRegistrationParameter,scalarD,t._g[2])
  //solver
  REGISTER_FLOAT_TYPE("sigma0",PhysicsRegistrationParameter,scalarD,t._sigma0)
  REGISTER_FLOAT_TYPE("cThres",PhysicsRegistrationParameter,scalarD,t._cThres)
  REGISTER_FLOAT_TYPE("alphaThres",PhysicsRegistrationParameter,scalarD,t._alphaThres)
  REGISTER_FLOAT_TYPE("newtonThres",PhysicsRegistrationParameter,scalarD,t._newtonThres)
  REGISTER_FLOAT_TYPE("initPenalty",PhysicsRegistrationParameter,scalarD,t._initPenalty)
  REGISTER_FLOAT_TYPE("maxPenalty",PhysicsRegistrationParameter,scalarD,t._maxPenalty)
  REGISTER_FLOAT_TYPE("lineSearchThres",PhysicsRegistrationParameter,scalarD,t._lineSearchThres)
  REGISTER_BOOL_TYPE("useAugLag",PhysicsRegistrationParameter,bool,t._useAugLag)
  REGISTER_BOOL_TYPE("callback",PhysicsRegistrationParameter,bool,t._callback)
  REGISTER_BOOL_TYPE("sparse",PhysicsRegistrationParameter,bool,t._sparse)
  REGISTER_INT_TYPE("maxIterNewton",PhysicsRegistrationParameter,sizeType,t._maxIterNewton)
  REGISTER_INT_TYPE("maxIterAugLag",PhysicsRegistrationParameter,sizeType,t._maxIterAugLag)
  reset(ops);
}
void PhysicsRegistrationParameter::reset(Options& ops)
{
  PhysicsRegistrationParameter::initOptions(*this);
  ops.setOptions(this);
}
void PhysicsRegistrationParameter::initOptions(PhysicsRegistrationParameter& sol)
{
  sol._coefPotential=0;
  sol._coefRegister=1;
  sol._tau=1;
  sol._g=Vec3d(0,0,-9.81f);
  //solver
  sol._sigma0=1;
  sol._cThres=1e-10f;
  sol._alphaThres=1e-20f;
  sol._newtonThres=1e-10f;
  sol._initPenalty=1e3f;
  sol._maxPenalty=1e5f;
  sol._lineSearchThres=2.0f;
  sol._useAugLag=true;
  sol._callback=true;
  sol._sparse=true;
  sol._maxIterNewton=2000;
  sol._maxIterAugLag=1000;
}
//Penetration
template <typename T>
bool PhysicsRegistration<T>::Penetration::operator<(const Penetration& other) const
{
  if(_oid<other._oid)
    return true;
  else if(_oid>other._oid)
    return false;

  if(_oidOther<other._oidOther)
    return true;
  else if(_oidOther>other._oidOther)
    return false;

  return false;
}
template <typename T>
bool PhysicsRegistration<T>::Penetration::operator!=(const Penetration& other) const {
  return *this<other || other<*this;
}
template <typename T>
bool PhysicsRegistration<T>::Penetration::operator==(const Penetration& other) const {
  return !(*this<other) && !(other<*this);
}
template <typename T>
void PhysicsRegistration<T>::Penetration::writeVTK(const std::string& path) const {
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  for(const std::tuple<sizeType,Vec3T,T>& t:_penetratedPoints)
    vss.push_back(std::get<1>(t).unaryExpr([&](const T& in) {
    return (scalar)std::to_double(in);
  }));
  VTKWriter<scalar> os("particles",path,true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size(),0,0),
                 VTKWriter<scalar>::POINT);
}
//PhysicsRegistration
template <typename T>
PhysicsRegistration<T>::PhysicsRegistration()
{
  _indexModifier=[&](sizeType,const Vec&) {};
}
template <typename T>
void PhysicsRegistration<T>::reset(const std::vector<ObjMesh>& objs,T rad,bool convex,T SDFRes,T SDFExtension,bool SDFRational)
{
  //build objects
  _body=ArticulatedBody();
  std::vector<ArticulatedBody> bodies;
  for(const ObjMesh& obj:objs) {
    bodies.push_back(ArticulatedLoader::createMesh(obj));
    ArticulatedUtils(bodies.back()).addBase(3,Vec3d::Zero());
  }
  ArticulatedUtils(_body).combine(bodies);
  ArticulatedUtils(_body).simplify(10);
  GraspPlanner<T>::reset(rad,convex,SDFRes,SDFExtension,SDFRational,false);

  _bvhss.resize(_pnss.size());
  _pLMax.resize(_pnss.size());
  for(sizeType i=0; i<(sizeType)_pnss.size(); i++) {
    //construct BVH for each object
    _bvhss[i].assign(_pnss[i].first.cols(),Node<sizeType,BBox<scalarD>>());
    if(_bvhss[i].empty())
      continue;
    for(sizeType j=0; j<_pnss[i].first.cols(); j++) {
      _bvhss[i][j]=Node<sizeType,BBox<scalarD>>();
      _bvhss[i][j]._cell=j;
      _bvhss[i][j]._nrCell=1;
      _bvhss[i][j]._bb.setUnion(_pnss[i].first.col(j).unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }));
    }
    buildBVH<sizeType>(_bvhss[i],3,-1);

    //compute _pLMax
    _pLMax[i]=0;
    for(sizeType j=0; j<_pnss[i].first.cols(); j++)
      _pLMax[i]=std::max(_pLMax[i],std::to_double(std::sqrt(_pnss[i].first.col(j).squaredNorm())));
  }
}
template <typename T>
bool PhysicsRegistration<T>::read(std::istream& is,IOData* dat)
{
  _indexModifier=[&](sizeType,const Vec&) {};
  GraspPlanner<T>::read(is,dat);
  readBinaryData(_bvhss,is,dat);
  readBinaryData(_pLMax,is,dat);
  return is.good();
}
template <typename T>
bool PhysicsRegistration<T>::write(std::ostream& os,IOData* dat) const
{
  GraspPlanner<T>::write(os,dat);
  writeBinaryData(_bvhss,os,dat);
  writeBinaryData(_pLMax,os,dat);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> PhysicsRegistration<T>::copy() const
{
  return std::shared_ptr<SerializableBase>(new PhysicsRegistration<T>);
}
template <typename T>
std::string PhysicsRegistration<T>::type() const
{
  return typeid(PhysicsRegistration<T>).name();
}
//index set helper
template <typename T>
void PhysicsRegistration<T>::setIndexModifier(std::function<void(sizeType,const Vec&)> func)
{
  _indexModifier=func;
}
template <typename T>
bool PhysicsRegistration<T>::existContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid) const
{
  if(!_objs.template getComponent<ContactIndexConstraint<T>>())
    return false;
  return _objs.template getComponent<ContactIndexConstraint<T>>()->existContactIndexMap(oid,oidOther,pid);
}
template <typename T>
bool PhysicsRegistration<T>::existPointCloudMap(sizeType pid,sizeType oid) const
{
  if(!_objs.template getComponent<PointCloudRegistrationEnergy<T>>())
    return false;
  return _objs.template getComponent<PointCloudRegistrationEnergy<T>>()->existPointCloudMap(pid,oid);
}
template <typename T>
void PhysicsRegistration<T>::addContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid,T mu,sizeType nDir)
{
  if(!_objs.template getComponent<ContactIndexConstraint<T>>())
    return;
  _objs.template getComponent<ContactIndexConstraint<T>>()->addContactIndexMap(oid,oidOther,pid,mu,nDir);
}
template <typename T>
void PhysicsRegistration<T>::addPointCloudMap(sizeType pid,sizeType oid,const Vec3T& objPos)
{
  if(!_objs.template getComponent<PointCloudRegistrationEnergy<T>>())
    return;
  _objs.template getComponent<PointCloudRegistrationEnergy<T>>()->addPointCloudMap(pid,oid,objPos);
}
template <typename T>
void PhysicsRegistration<T>::clearContactIndexMap()
{
  if(!_objs.template getComponent<ContactIndexConstraint<T>>())
    return;
  _objs.template getComponent<ContactIndexConstraint<T>>()->clearContactIndexMap();
}
template <typename T>
void PhysicsRegistration<T>::clearPointCloudMap()
{
  if(!_objs.template getComponent<PointCloudRegistrationEnergy<T>>())
    return;
  _objs.template getComponent<PointCloudRegistrationEnergy<T>>()->clearPointCloudMap();
}
template <typename T>
void PhysicsRegistration<T>::writeContactVTK(const Vec& x,const std::string& path) const
{
  const PBDArticulatedGradientInfo<T> info(_body,x.segment(0,_body.nrDOF()));
  _objs.template getComponent<ContactIndexConstraint<T>>()->writeContactVTK(x,info,path);
}
template <typename T>
void PhysicsRegistration<T>::getDeepestPenetration(std::set<Penetration>& pss) const
{
  pss.clear();
  for(sizeType i=-1; i<nrObject(); i++)
    for(sizeType j=0; j<nrObject(); j++) {
      Penetration p;
      p._oid=i;
      p._oidOther=j;
      if(i==j)
        continue;
      getDeepestPenetration(p);
      if(!p._penetratedPoints.empty())
        pss.insert(p);
    }
}
template <typename T>
void PhysicsRegistration<T>::getDeepestPenetration(Penetration& p) const
{
  const PointCloudObject<T>& obj=std::dynamic_pointer_cast<ArticulatedObjective<T>>(_objs.components().begin()->second)->object();
  const PBDArticulatedGradientInfo<T>& info=std::dynamic_pointer_cast<ArticulatedObjective<T>>(_objs.components().begin()->second)->info();
  const std::vector<Node<sizeType,BBox<scalarD>>>& bvh=p._oid<0?obj.getBVH():_bvhss[p._oid*2+2];
  const Environment<T>& e=env(p._oidOther*2+2);
  OBBTpl<scalarD,3> bbOther(ROTI(info._TM,p._oidOther*2+2).unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  }),CTRI(info._TM,p._oidOther*2+2).unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  }),e.getBB());

  p._penetratedPoints.clear();
  std::get<2>(p._deepestPenetration)=ScalarUtil<T>::scalar_max();
  std::stack<sizeType> ss;
  ss.push(bvh.size()-1);
  while(!ss.empty()) {
    sizeType pid=ss.top();
    ss.pop();
    OBBTpl<scalarD,3> bb;
    if(p._oid>=0) {
      bb=OBBTpl<scalarD,3>(ROTI(info._TM,p._oid*2+2).unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),CTRI(info._TM,p._oid*2+2).unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),bvh[pid]._bb);
    } else bb=bvh[pid]._bb;
    if(!bb.intersect(bbOther))
      continue;
    else if(bvh[pid]._cell>=0) {
      Vec3T pt,ptG,ptL;

      if(p._oid>=0) {
        pt=_pnss[p._oid*2+2].first.col(bvh[pid]._cell);
        ptG=ROTI(info._TM,p._oid*2+2)*pt+CTRI(info._TM,p._oid*2+2);
      } else pt=ptG=obj.pss().col(bvh[pid]._cell);

      T phi=e.phi(ROTI(info._TM,p._oidOther*2+2).transpose()*(ptG-CTRI(info._TM,p._oidOther*2+2)));
      if(phi<0) {
        p._penetratedPoints.push_back(std::make_tuple(bvh[pid]._cell,ptG,phi));
        if(std::get<2>(p._deepestPenetration)>phi)
          p._deepestPenetration=p._penetratedPoints.back();
      }
    } else {
      ss.push(bvh[pid]._l);
      ss.push(bvh[pid]._r);
    }
  }
}
template <typename T>
void PhysicsRegistration<T>::debugPenetration(sizeType iter,T scale)
{
  std::set<Penetration> pss;
  for(sizeType it=0; it<iter;) {
    const PointCloudObject<T>& obj=std::dynamic_pointer_cast<ArticulatedObjective<T>>(_objs.components().begin()->second)->object();
    PBDArticulatedGradientInfo<T>& info=std::dynamic_pointer_cast<ArticulatedObjective<T>>(_objs.components().begin()->second)->info();
    info.reset(_body,Vec::Random(_body.nrDOF())*scale);
    getDeepestPenetration(pss);
    if(!pss.empty()) {
      recreate("penetration"+std::to_string(it));
      _body.writeVTK(info._TM.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),"penetration"+std::to_string(it)+"/body.vtk",Joint::MESH);
      obj.writeVTK("penetration"+std::to_string(it)+"/object",0);
      for(const Penetration& p:pss)
        p.writeVTK("penetration"+std::to_string(it)+"/O"+std::to_string(p._oid)+"OO"+std::to_string(p._oidOther)+".vtk");
      it++;
    }
  }
}
template <typename T>
sizeType PhysicsRegistration<T>::nrObject() const
{
  return _body.nrJ()/2;
}
//optimize
template <typename T>
typename PhysicsRegistration<T>::Vec PhysicsRegistration<T>::optimize(bool debug,const Vec& init,const PointCloudObject<T>& env,const PointCloudObject<T>* obj,PhysicsRegistrationParameter& ops)
{
  //build objective function
  _objs=DSSQPObjectiveCompound<T>();
  _info=PBDArticulatedGradientInfo<T>();
  if(ops._tau>0)
    _objs.addComponent(std::shared_ptr<ContactIndexConstraint<T>>(new ContactIndexConstraint<T>(_objs,_info,*this,env,ops._g.template cast<T>(),ops._tau)));
  if(ops._coefRegister>0 && obj)
    _objs.addComponent(std::shared_ptr<PointCloudRegistrationEnergy<T>>(new PointCloudRegistrationEnergy<T>(_objs,_info,*this,*obj,ops._coefRegister)));
  if(ops._coefPotential>0)
    _objs.addComponent(std::shared_ptr<GravityEnergy<T>>(new GravityEnergy<T>(_objs,_info,*this,env,ops._g.template cast<T>(),ops._coefPotential)));

  Vec x;
  SolveNewton<T>::template solveNewton<Vec>(_A.transpose()*_A,_A.transpose()*(_b-init),x,true);
  clearContactIndexMap();
  //clearPointCloudMap();

  _b.setZero(x.size());
  _A=MatT::Identity(x.size(),x.size()).sparseView();
  _l=_objs.lb();
  _u=_objs.ub();
  _gl=_objs.gl();
  _gu=_objs.gu();

  if(debug)
    debugSystemAugLag(x);
  else if(ops._useAugLag)
    x=optimizeAugLag(x,ops);
  else x=optimizeSQP(x,ops);
  return x;
}
template <typename T>
bool PhysicsRegistration<T>::assembleAugLag(Vec x,const std::vector<T>& lambda,T penalty,bool update,T& e,Vec* g,MatT* h)
{
  Vec c;
  MatT cjac;
  bool succ=GraspPlanner<T>::assemble(x,update,e,g,h,&c,g?&cjac:NULL);
  if(!succ)
    return false;
  for(sizeType i=0; i<c.size(); i++) {
    T CL=c[i]-_gl[i];
    T CU=_gu[i]-c[i];
    if(CL<lambda[i]/penalty) {
      e+=-lambda[i]*CL+0.5f*penalty*CL*CL;
      if(g)
        *g+=cjac.row(i).transpose()*(penalty*CL-lambda[i]);
      if(h)
        *h+=cjac.row(i).transpose()*penalty*cjac.row(i);
    } else if(CU<lambda[i]/penalty) {
      e+=-lambda[i]*CU+0.5f*penalty*CU*CU;
      if(g)
        *g-=cjac.row(i).transpose()*(penalty*CU-lambda[i]);
      if(h)
        *h+=cjac.row(i).transpose()*penalty*cjac.row(i);
    } else {
      e+=-0.5f*lambda[i]*lambda[i]/penalty;
    }
  }
  return true;
}
template <typename T>
bool PhysicsRegistration<T>::assembleAugLag(Vec x,const std::vector<T>& lambda,T penalty,bool update,T& e,Vec* g,SMat* h)
{
  Vec c;
  SMat cjac;
  bool succ=GraspPlanner<T>::assemble(x,update,e,g,h,&c,g?&cjac:NULL);
  if(!succ)
    return false;
  for(sizeType i=0; i<c.size(); i++) {
    T CL=c[i]-_gl[i];
    T CU=_gu[i]-c[i];
    if(CL<lambda[i]/penalty) {
      e+=-lambda[i]*CL+0.5f*penalty*CL*CL;
      if(g)
        *g+=cjac.row(i).transpose()*(penalty*CL-lambda[i]);
      if(h)
        *h+=cjac.row(i).transpose()*penalty*cjac.row(i);
    } else if(CU<lambda[i]/penalty) {
      e+=-lambda[i]*CU+0.5f*penalty*CU*CU;
      if(g)
        *g-=cjac.row(i).transpose()*(penalty*CU-lambda[i]);
      if(h)
        *h+=cjac.row(i).transpose()*penalty*cjac.row(i);
    } else {
      e+=-0.5f*lambda[i]*lambda[i]/penalty;
    }
  }
  return true;
}
template <typename T>
bool PhysicsRegistration<T>::lineSearchThresViolated(const Vec& x,const Vec& xNew,const PhysicsRegistrationParameter& ops) const
{
  PBDArticulatedGradientInfo<T> info(_body,x);
  PBDArticulatedGradientInfo<T> infoNew(_body,xNew);
  for(sizeType i=2; i<_body.nrJ(); i+=2) {
    Vec3T dCtr=CTRI(infoNew._TM,i)-CTRI(info._TM,i);
    Vec3T dRot=invExpW<T>(ROTI(infoNew._TM,i).transpose()*ROTI(info._TM,i));
    T delta=std::sqrt(dCtr.squaredNorm())+std::sqrt(dRot.squaredNorm())*_pLMax[i];
    if(delta>ops._lineSearchThres*rad())
      return true;
  }
  return false;
}
template <typename T>
typename PhysicsRegistration<T>::Vec PhysicsRegistration<T>::optimizeNewton(Vec x,std::vector<T>& lambda,T penalty,PhysicsRegistrationParameter& ops)
{
  Vec d;
  T e,e2;
  MatT hD;
  SMat hS;
  Vec g,xTmp;
  T dNorm,alphaDec=0.5f,alphaInc=1.5f,coefWolfe=0.1f,alpha=1,reg=0;
  _gl=_objs.gl(),_gu=_objs.gu();

  for(sizeType it=0; it<ops._maxIterNewton; it++) {
    //update index
    for(typename std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>::const_iterator beg=_objs.components().begin(),end=_objs.components().end(); beg!=end; beg++)
      beg->second->setUpdateCache(x,true);
    _indexModifier(it,x);
    bool updated=false;
    while((sizeType)lambda.size()<_objs.values()) {
      updated=true;
      lambda.push_back(0);
    }
    if(updated) {
      if(_objs.inputs()>x.size())
        x=concat(x,Vec::Zero(_objs.inputs()-x.size()));
      _b.setZero(x.size());
      _A=MatT::Identity(x.size(),x.size()).sparseView();
      _l=_objs.lb();
      _u=_objs.ub();
      _gl=_objs.gl();
      _gu=_objs.gu();
    }
    //solve system
    //debugSystemAugLag(x);
    if(ops._sparse) {
      if(!assembleAugLag(x,lambda,penalty,true,e,&g,&hS)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(invalid configuration)",it)
        }
        return Vec::Zero(0);
      }
      if(reg==0)
        reg=std::max<T>(1e-3f,hS.toDense().diagonal().unaryExpr([&](const T& in) {
        return (scalarD)std::abs(in);
      }).maxCoeff());
      if(!solveSparseQP(d,x,g,hS,NULL,NULL,0,0,reg)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(qp failed)",it)
        }
        break;
      }
    } else {
      if(!assembleAugLag(x,lambda,penalty,true,e,&g,&hD)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(invalid configuration)",it)
        }
        return Vec::Zero(0);
      }
      if(!solveDenseQP(d,x,g,hD,NULL,NULL,0,0)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(qp failed)",it)
        }
        break;
      }
    }
    //termination & callback
    dNorm=std::sqrt(d.squaredNorm());
    if(dNorm<ops._newtonThres) {
      if(ops._callback) {
        INFOV("Iter=%d succeed(dNorm=%f<thres=%f)",it,std::to_double(dNorm),std::to_double(ops._newtonThres))
      }
      break;
    } else if(ops._callback) {
      INFOV("Iter=%d E=%f dNorm=%f alpha=%f (inputs=%d,values=%d)",it,std::to_double(e),std::to_double(dNorm),std::to_double(alpha),_objs.inputs(),_objs.values())
    }
    //line search
    while(alpha>ops._alphaThres) {
      xTmp=x+d*alpha;
      if(lineSearchThresViolated(x,xTmp,ops)) {
        alpha*=alphaDec;
        continue;
      }
      if(!assembleAugLag(xTmp,lambda,penalty,false,e2,NULL,(DMat*)NULL)) {
        alpha*=alphaDec;
        continue;
      }
      if(e2>=e+g.dot(d.template cast<T>())*alpha*coefWolfe) {
        alpha*=alphaDec;
        continue;
      } else {
        alpha=std::min<T>(alpha*alphaInc,1);
        x=xTmp;
        break;
      }
    }
    if(alpha<ops._alphaThres) {
      if(ops._callback) {
        INFOV("Iter=%d failed(alpha=%f<alphaThres=%f)",it,std::to_double(alpha),std::to_double(ops._alphaThres))
      }
      break;
    }
  }
  return x;
}
template <typename T>
typename PhysicsRegistration<T>::Vec PhysicsRegistration<T>::optimizeAugLag(Vec x,PhysicsRegistrationParameter& ops)
{
  Vec c;
  T penalty=ops._initPenalty,e;
  std::vector<T> lambda(_objs.values(),0);
  T lastCNorm=ScalarUtil<T>::scalar_max(),cNorm;
  for(sizeType it=0; it<ops._maxIterAugLag; it++) {
    x=optimizeNewton(x,lambda,penalty,ops);
    if(x.size()==0) {
      if(ops._callback) {
        INFOV("AugLagIter=%d failed(optimizeNewton failed)",it)
      }
      break;
    }
    //acquire constraint
    if(!GraspPlanner<T>::assemble(x,true,e,NULL,(SMat*)NULL,&c,NULL)) {
      if(ops._callback) {
        INFOV("Iter=%d failed(invalid configuration)",it)
      }
      return Vec::Zero(0);
    }
    cNorm=0;
    for(sizeType i=0; i<c.size(); i++)
      if(c[i]<_gl[i])
        cNorm=std::min<T>(cNorm,_gl[i]-c[i]);
      else if(c[i]>_gu[i])
        cNorm=std::max<T>(cNorm,c[i]-_gu[i]);
    if(cNorm<ops._cThres) {
      if(ops._callback) {
        INFOV("AugLagIter=%d succeed(cNorm=%f<thres=%f)",it,std::to_double(cNorm),std::to_double(ops._cThres))
      }
      break;
    }
    //update penalty
    if(cNorm>lastCNorm*0.25f)
      penalty=std::max<T>(std::pow(it+1,2),penalty*10);
    lastCNorm=cNorm;
    if(ops._callback) {
      INFOV("AugLagIter=%d cNorm=%f penalty=%f",it,std::to_double(cNorm),std::to_double(penalty))
    }
    if(penalty>ops._maxPenalty) {
      if(ops._callback) {
        INFOV("Iter=%d failed(max penalty=%f)",it,ops._maxPenalty)
      }
      break;
    }
    //update lambda
    for(sizeType i=0; i<c.size(); i++) {
      T CL=c[i]-_gl[i];
      T CU=_gu[i]-c[i];
      if(CL<0)
        lambda[i]-=CL*penalty;
      else if(CU<0)
        lambda[i]-=CU*penalty;
      else lambda[i]=0;
    }
  }
  return x;
}
template <typename T>
typename PhysicsRegistration<T>::Vec PhysicsRegistration<T>::optimizeSQP(Vec x,PhysicsRegistrationParameter& ops)
{
  Vec d;
  T e,e2;
  MatT hD,cjacD;
  SMat hS,cjacS;
  Vec g,c,c2,xTmp;
  T realReduction,predReduction;
  T cNorm,dNorm,sigma=ops._sigma0,reg=0,TR=rad()*x.size(),rho=0,TRDec=0.5f;
  _gl=_objs.gl(),_gu=_objs.gu();

  for(sizeType it=0; it<ops._maxIterNewton; it++) {
    //solve system
    if(ops._sparse) {
      if(!GraspPlanner<T>::assemble(x,true,e,&g,&hS,&c,&cjacS)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(invalid configuration)",it)
        }
        return Vec::Zero(0);
      }
      if(reg==0)
        reg=std::max<T>(1e-3f,hS.toDense().diagonal().unaryExpr([&](const T& in) {
        return (scalarD)std::abs(in);
      }).maxCoeff());
      if(!solveSparseQP(d,x,g,hS,&c,&cjacS,TR,sigma,reg)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(qp failed)",it)
        }
        break;
      }
    } else {
      if(!GraspPlanner<T>::assemble(x,true,e,&g,&hD,&c,&cjacD)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(invalid configuration)",it)
        }
        return Vec::Zero(0);
      }
      if(!solveDenseQP(d,x,g,hD,&c,&cjacD,TR,sigma)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(qp failed)",it)
        }
        break;
      }
    }
    //termination & callback
    if(ops._sparse)
      cNorm=computeLinearizedCInf(d,c,cjacS);
    else cNorm=computeLinearizedCInf(d,c,cjacD);
    dNorm=std::sqrt(d.squaredNorm());
    if(cNorm<ops._cThres) {
      //QP is consistent, keep sigma unchanged
      if(dNorm<ops._newtonThres) {
        if(ops._callback) {
          INFOV("Iter=%d succeed(dNorm=%f<thres=%f)",it,std::to_double(dNorm),std::to_double(ops._newtonThres))
        }
        break;
      }
    } else {
      //QP is inconsistent, we increase sigma
      sigma=std::min<T>(sigma*2,ops._maxPenalty);
    }
    //compute rho
    {
      GraspPlanner<T>::assemble(x+d,false,e2,(Vec*)NULL,(DMat*)NULL,&c2);
      realReduction=computePhi(e,c,sigma)-computePhi(e2,c2,sigma),predReduction=0;
      if(ops._sparse)
        predReduction=computePhiPred(Vec::Zero(g.size()),e,hS,g,c,cjacS,sigma)-computePhiPred(d,e,hS,g,c,cjacS,sigma);
      else
        predReduction=computePhiPred(Vec::Zero(g.size()),e,hD,g,c,cjacD,sigma)-computePhiPred(d,e,hD,g,c,cjacD,sigma);
      rho=realReduction/predReduction;
    }
    //update trust region
    if(lineSearchThresViolated(x,x+d,ops)) {
      TR*=TRDec;    //avoid tunneling
      if(ops._callback) {
        INFOV("Iter=%d E=%f cNorm=%f dNorm=%f sigma=%f predReduction=%f realReduction=%f TR=%f (inputs=%d,values=%d,lineSearchThresViolated)",
              it,std::to_double(e2),std::to_double(cNorm),std::to_double(dNorm),std::to_double(sigma),std::to_double(predReduction),std::to_double(realReduction),std::to_double(TR),_objs.inputs(),_objs.values())
      }
    } else if(realReduction>0 && predReduction>0) {
      x+=d; //accept, adjust trust region
      TR*=std::min<T>(std::max<T>(std::pow(2*rho-1,3)+1,0.25f),2);
      if(ops._callback) {
        INFOV("Iter=%d E=%f cNorm=%f dNorm=%f sigma=%f predReduction=%f realReduction=%f TR=%f (inputs=%d,values=%d,accepted)",
              it,std::to_double(e2),std::to_double(cNorm),std::to_double(dNorm),std::to_double(sigma),std::to_double(predReduction),std::to_double(realReduction),std::to_double(TR),_objs.inputs(),_objs.values())
      }
      //update index is only allowed after a successful trust region step
      for(typename std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>::const_iterator beg=_objs.components().begin(),end=_objs.components().end(); beg!=end; beg++)
        beg->second->setUpdateCache(x,true);
      _indexModifier(it,x);
      bool updated=false;
      if(_gl.size()<_objs.values())
        updated=true;
      if(updated) {
        if(_objs.inputs()>x.size())
          x=concat(x,Vec::Zero(_objs.inputs()-x.size()));
        _b.setZero(x.size());
        _A=MatT::Identity(x.size(),x.size()).sparseView();
        _l=_objs.lb();
        _u=_objs.ub();
        _gl=_objs.gl();
        _gu=_objs.gu();
      }
    } else {
      //reject, shrink trust region
      TR*=std::min<T>(std::max<T>(std::pow(2*rho-1,3)+1,0.25f),2);
      if(ops._callback) {
        INFOV("Iter=%d E=%f cNorm=%f dNorm=%f sigma=%f predReduction=%f realReduction=%f TR=%f (inputs=%d,values=%d,rejected)",
              it,std::to_double(e2),std::to_double(cNorm),std::to_double(dNorm),std::to_double(sigma),std::to_double(predReduction),std::to_double(realReduction),std::to_double(TR),_objs.inputs(),_objs.values())
      }
    }
    if(TR<ops._alphaThres) {
      if(ops._callback) {
        INFOV("Iter=%d failed(TR=%f<alphaThres=%f)",it,std::to_double(TR),std::to_double(ops._alphaThres))
      }
      break;
    }
  }
  return x;
}
template <typename T>
T PhysicsRegistration<T>::computePhi(T e,const Vec& c,T sigma) const
{
  for(sizeType i=0; i<c.size(); i++)
    if(c[i]<_gl[i])
      e+=sigma*(_gl[i]-c[i]);
    else if(c[i]>_gu[i])
      e+=sigma*(c[i]-_gu[i]);
  return e;
}
template <typename T>
template <typename MAT>
T PhysicsRegistration<T>::computePhiPred(const Vec& d,T e,const MAT& h,const Vec& g,Vec c,const MAT& cjac,T sigma) const
{
  e+=d.dot(h*d/2+g);
  c+=cjac*d;
  return computePhi(e,c,sigma);
}
template <typename T>
template <typename MAT>
T PhysicsRegistration<T>::computeLinearizedCInf(const Vec& d,const Vec& c,const MAT& cjac) const
{
  T val=0;
  for(sizeType i=0; i<c.size(); i++) {
    T ci=(c+cjac*d)[i];
    if(ci<_gl[i])
      val=std::max<T>(val,_gl[i]-ci);
    else if(ci>_gu[i])
      val=std::max<T>(val,ci-_gu[i]);
  }
  return val;
}
template <typename T>
void PhysicsRegistration<T>::debugSystemAugLag(const Vec& x)
{
  DEFINE_NUMERIC_DELTA_T(T)
  Vec dx=Vec::Random(x.size());
  _gl=_objs.gl(),_gu=_objs.gu();
  //first evaluate
  T e;
  Vec g,c;
  MatT hD;
  SMat hS;
  T penalty=RandEngine::randR(0,1);
  std::vector<T> lambda(_objs.values());
  for(sizeType i=0; i<_objs.values(); i++)
    lambda[i]=RandEngine::randR(0,1);
  assembleAugLag(x,lambda,penalty,true,e,&g,&hD);
  assembleAugLag(x,lambda,penalty,true,e,&g,&hS);
  //second evaluate
  T e2;
  Vec g2,c2;
  assembleAugLag(x+dx*DELTA,lambda,penalty,false,e2,&g2,(DMat*)NULL);
  //compare
  DEBUG_GRADIENT("PhysicsRegistration-G",g.dot(dx),g.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT("PhysicsRegistration-HDense",std::sqrt((hD*dx).squaredNorm()),std::sqrt((hD*dx-(g2-g)/DELTA).squaredNorm()))
  DEBUG_GRADIENT("PhysicsRegistration-HSparse",std::sqrt((hS*dx).squaredNorm()),std::sqrt((hS*dx-(g2-g)/DELTA).squaredNorm()))
}
//instance
PRJ_BEGIN
template class PhysicsRegistration<double>;
#ifdef ALL_TYPES
template class PhysicsRegistration<__float128>;
template class PhysicsRegistration<mpfr::mpreal>;
#endif
PRJ_END
