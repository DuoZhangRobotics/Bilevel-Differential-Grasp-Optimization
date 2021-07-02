#include <Utils/Scalar.h>
#include "GraspPlanner.h"
#include <Utils/Utils.h>
#include <Utils/SparseUtils.h>
#include <Utils/DebugGradient.h>
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/MultiPrecisionLQP.h>
#include <CommonFile/ParallelPoissonDiskSampling.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/Environment.h>
#include <CommonFile/Timing.h>
#include <Eigen/Sparse>
#include <Eigen/Eigen>
//energy
#include <Optimizer/DSSQPObjective.h>
#include "ConvexLogBarrierSelfEnergy.h"
#include "PrimalDualQInfMetricEnergyFGT.h"
#include "PrimalDualQInfMetricEnergy.h"
#include "CentroidClosednessEnergy.h"
#include "ObjectClosednessEnergy.h"
#include "LogBarrierObjEnergy.h"
#include "MetricEnergy.h"

USE_PRJ_NAMESPACE

//GraspPlannerParameter
GraspPlannerParameter::GraspPlannerParameter(Options& ops)
{
  REGISTER_FLOAT_TYPE("d0",GraspPlannerParameter,scalarD,t._d0)
  REGISTER_FLOAT_TYPE("alpha",GraspPlannerParameter,scalarD,t._alpha)
  REGISTER_INT_TYPE("metric",GraspPlannerParameter,sizeType,t._metric)
  REGISTER_INT_TYPE("activation",GraspPlannerParameter,sizeType,t._activation)
  REGISTER_FLOAT_TYPE("normalExtrude",GraspPlannerParameter,scalarD,t._normalExtrude)
  REGISTER_FLOAT_TYPE("FGTThres",GraspPlannerParameter,scalarD,t._FGTThres)
  REGISTER_FLOAT_TYPE("coefM",GraspPlannerParameter,scalarD,t._coefM)
  REGISTER_FLOAT_TYPE("coefOC",GraspPlannerParameter,scalarD,t._coefOC)
  REGISTER_FLOAT_TYPE("coefCC",GraspPlannerParameter,scalarD,t._coefCC)
  REGISTER_FLOAT_TYPE("coefO",GraspPlannerParameter,scalarD,t._coefO)
  REGISTER_FLOAT_TYPE("coefS",GraspPlannerParameter,scalarD,t._coefS)
  REGISTER_FLOAT_TYPE("useGJK",GraspPlannerParameter,bool,t._useGJK)
  //solver
  REGISTER_FLOAT_TYPE("rho0",GraspPlannerParameter,scalarD,t._rho0)
  REGISTER_FLOAT_TYPE("thres",GraspPlannerParameter,scalarD,t._thres)
  REGISTER_FLOAT_TYPE("alphaThres",GraspPlannerParameter,scalarD,t._alphaThres)
  REGISTER_BOOL_TYPE("callback",GraspPlannerParameter,bool,t._callback)
  REGISTER_BOOL_TYPE("sparse",GraspPlannerParameter,bool,t._sparse)
  REGISTER_INT_TYPE("maxIter",GraspPlannerParameter,sizeType,t._maxIter)
  reset(ops);
}
void GraspPlannerParameter::reset(Options& ops)
{
  GraspPlannerParameter::initOptions(*this);
  ops.setOptions(this);
}
void GraspPlannerParameter::initOptions(GraspPlannerParameter& sol)
{
  sol._d0=1;
  sol._alpha=1e-3f;
  sol._metric=Q_INF_CONSTRAINT;
  sol._activation=SQR_EXP_ACTIVATION;
  sol._normalExtrude=1;
  sol._FGTThres=1e-6f;
  sol._coefM=-1;
  sol._coefOC=0;
  sol._coefCC=0;
  sol._coefO=100;
  sol._coefS=1;
  sol._useGJK=false;
  //solver
  sol._rho0=1;
  sol._thres=1e-10f;
  sol._alphaThres=1e-20f;
  sol._callback=true;
  sol._sparse=false;
  sol._maxIter=2000;
}
//GraspPlanner
template <typename T>
GraspPlanner<T>::GraspPlanner() {}
template <typename T>
void GraspPlanner<T>::reset(T rad,bool convex,T SDFRes,T SDFExtension,bool SDFRational,bool checkValid)
{
  _pnss.resize(_body.nrJ());
  _rad=rad;
  //reset geom
  ObjMesh m;
  _body.getGeom().clear();
  for(sizeType i=0; i<_body.nrJ(); i++) {
    std::shared_ptr<StaticGeomCell> cell=_body.joint(i).getGeomPtr();
    if(!cell)
      cell.reset(new ObjMeshGeomCell(Mat4::Identity(),m,0));
    _body.getGeom().addGeomCell(cell);
  }
  _body.getGeom().assemble();
  //mimic
  typename ArticulatedBody::MatT A;
  typename ArticulatedBody::Vec b,l,u;
  _body.mimic(A,b,l,u);
  for(sizeType i=0; i<l.size(); i++) {
    if(!std::isfinite(l[i]))
      l[i]=-std::to_double(DSSQPObjective<T>::infty());
    if(!std::isfinite(u[i]))
      u[i]= std::to_double(DSSQPObjective<T>::infty());
  }
  _A=A.template cast<T>().sparseView();
  _b=b.template cast<T>();
  _l=l.template cast<T>();
  _u=u.template cast<T>();
  //distExact
  _env.resize(_body.nrJ());
  for(sizeType i=0; i<_body.nrJ(); i++)
    if(SDFRes>0) {
      if(SDFRational)
        _env[i].reset(new EnvironmentCubic<T>(ObjMeshGeomCellExact(dynamic_cast<const ObjMeshGeomCell&>(_body.getGeom().getG(i))),std::to_double(SDFRes),std::to_double(SDFExtension)));
      else _env[i].reset(new EnvironmentCubic<T>(dynamic_cast<const ObjMeshGeomCell&>(_body.getGeom().getG(i)),std::to_double(SDFRes),std::to_double(SDFExtension)));
    } else if(convex)
      _env[i].reset(new EnvironmentExact<T>(ConvexHullExact(dynamic_cast<const ObjMeshGeomCell&>(_body.getGeom().getG(i)))));
    else _env[i].reset(new EnvironmentExact<T>(ObjMeshGeomCellExact(dynamic_cast<const ObjMeshGeomCell&>(_body.getGeom().getG(i)))));
  //sample
  PBDArticulatedGradientInfo<T> info(_body,Vec::Zero(_body.nrDOF()));
  for(sizeType i=0; i<_body.nrJ(); i++) {
    _body.getGeom().getG(i).getMesh(m);
    if(m.getV().empty()) {
      continue;
    }

    //rescale mesh to gain accuracy
    ParallelPoissonDiskSampling sampler(3);
    sampler.setRadius(std::to_double(_rad));
    sampler.sample(m);
    //pss
    sizeType k=0;
    _pnss[i].first.resize(3,sampler.getPSet().size());
    _pnss[i].second.resize(3,sampler.getPSet().size());
    for(sizeType j=0; j<sampler.getPSet().size(); j++) {
      Vec3T v=sampler.getPSet()[j]._pos.template cast<T>();
      Vec3T n=sampler.getPSet()[j]._normal.template cast<T>();
      if(checkValid && !validSample(i,info,ROTI(info._TM,i)*v+CTRI(info._TM,i)))
        continue;
      _pnss[i].first.col(k)=v;
      _pnss[i].second.col(k)=n;
      k++;
    }
    _pnss[i].first=_pnss[i].first.block(0,0,3,k).eval();
    _pnss[i].second=_pnss[i].second.block(0,0,3,k).eval();
  }
}
template <typename T>
void GraspPlanner<T>::reset(const std::string& path,T rad,bool convex,T SDFRes,T SDFExtension,bool SDFRational)
{
  _body=ArticulatedLoader().readURDF(path,convex,true);
  ArticulatedUtils(_body).addBase(3,Vec3d::Zero());
  ArticulatedUtils(_body).simplify(10);
  reset(rad,convex,SDFRes,SDFExtension,SDFRational);
}
template <typename T>
void GraspPlanner<T>::fliterSample(std::function<bool(sizeType lid,const Vec3T& p,const Vec3T& n)> f)
{
  for(sizeType i=0; i<(sizeType)_pnss.size(); i++) {
    sizeType k=0;
    std::pair<Mat3XT,Mat3XT> pn=_pnss[i];
    for(sizeType j=0; j<pn.first.cols(); j++)
      if(f(i,pn.first.col(j),pn.second.col(j))) {
        pn.first.col(k)=pn.first.col(j);
        pn.second.col(k)=pn.second.col(j);
        k++;
      }
    _pnss[i]=std::make_pair<Mat3XT,Mat3XT>(pn.first.block(0,0,3,k),pn.second.block(0,0,3,k));
  }
}
template <typename T>
bool GraspPlanner<T>::read(std::istream& is,IOData* dat)
{
  registerType<EnvironmentCubic<T>>(dat);
  registerType<EnvironmentExact<T>>(dat);
  registerType<GraspPlanner<T>>(dat);
  readBinaryData(_env,is,dat);
  _body.read(is,dat);
  //mimic
  readBinaryData(_A,is);
  readBinaryData(_b,is);
  readBinaryData(_l,is);
  readBinaryData(_u,is);
  //do not serialize _gl,_gu, this can change due to number of constraints user added
  //readBinaryData(_gl,is);
  //readBinaryData(_gu,is);
  //sample
  readBinaryData(_pnss,is);
  readBinaryData(_rad,is);
  return is.good();
}
template <typename T>
bool GraspPlanner<T>::write(std::ostream& os,IOData* dat) const
{
  registerType<EnvironmentCubic<T>>(dat);
  registerType<EnvironmentExact<T>>(dat);
  registerType<GraspPlanner<T>>(dat);
  writeBinaryData(_env,os,dat);
  _body.write(os,dat);
  //mimic
  writeBinaryData(_A,os);
  writeBinaryData(_b,os);
  writeBinaryData(_l,os);
  writeBinaryData(_u,os);
  //do not serialize _gl,_gu, this can change due to number of constraints user added
  //writeBinaryData(_gl,os);
  //writeBinaryData(_gu,os);
  //sample
  writeBinaryData(_pnss,os);
  writeBinaryData(_rad,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> GraspPlanner<T>::copy() const
{
  return std::shared_ptr<SerializableBase>(new GraspPlanner<T>);
}
template <typename T>
std::string GraspPlanner<T>::type() const
{
  return typeid(GraspPlanner<T>).name();
}
template <typename T>
ArticulatedBody& GraspPlanner<T>::body()
{
  return _body;
}
template <typename T>
const ArticulatedBody& GraspPlanner<T>::body() const
{
  return _body;
}
template <typename T>
const Environment<T>& GraspPlanner<T>::env(sizeType jid) const
{
  return *(_env.at(jid));
}
template <typename T>
const ObjMeshGeomCellExact& GraspPlanner<T>::dist(sizeType jid) const
{
  return std::dynamic_pointer_cast<EnvironmentExact<T>>(_env.at(jid))->getObj();
}
template <typename T>
const typename GraspPlanner<T>::PNSS& GraspPlanner<T>::pnss() const
{
  return _pnss;
}
template <typename T>
void GraspPlanner<T>::writeVTK(const Vec& x,const std::string& path,T len) const
{
  sizeType lid=0;
  PBDArticulatedGradientInfo<T> info(_body,x);
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  for(const std::pair<Mat3XT,Mat3XT>& pn:_pnss) {
    for(sizeType i=0; i<pn.first.cols(); i++) {
      Vec3T p=ROTI(info._TM,lid)*pn.first.col(i)+CTRI(info._TM,lid);
      Vec3T n=ROTI(info._TM,lid)*pn.second.col(i);
      vss.push_back(p.unaryExpr([&](const T& in) {
        return (scalar)std::to_double(in);
      }));
      vss.push_back((p+n*len*_rad).unaryExpr([&](const T& in) {
        return (scalar)std::to_double(in);
      }));
    }
    lid++;
  }
  create(path);
  VTKWriter<scalar> os("particles",path+"/sample.vtk",true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                 VTKWriter<scalar>::POINT);
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                 VTKWriter<scalar>::LINE);
  _body.writeVTK(info._TM.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  }),path+"/body.vtk",Joint::MESH);
}
template <typename T>
void GraspPlanner<T>::writeLocalVTK(const std::string& path,T len) const
{
  create(path);
  for(sizeType j=0; j<_body.nrJ(); j++)
    if(_pnss[j].first.cols()>0) {
      std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
      for(sizeType i=0; i<_pnss[j].first.cols(); i++) {
        vss.push_back(_pnss[j].first.col(i).unaryExpr([&](const T& in) {
          return (scalar)std::to_double(in);
        }));
        vss.push_back((_pnss[j].first.col(i)+_pnss[j].second.col(i)*len*_rad).unaryExpr([&](const T& in) {
          return (scalar)std::to_double(in);
        }));
      }
      VTKWriter<scalar> os("particles",path+"/joint"+std::to_string(j)+".vtk",true);
      os.appendPoints(vss.begin(),vss.end());
      os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                     VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                     VTKWriter<scalar>::POINT);
      os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                     VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                     VTKWriter<scalar>::LINE);
      _env[j]->getMesh().writeVTK(path+"/jointMesh"+std::to_string(j)+".vtk",true);
    }
}
template <typename T>
void GraspPlanner<T>::writeLimitsVTK(const std::string& path) const
{
  create(path);
  for(sizeType i=0; i<_l.size(); i++) {
    if(std::isfinite(_l[i])) {
      Vec x=Vec::Unit(_l.size(),i)*_l[i];
      PBDArticulatedGradientInfo<T> info(_body,_A*x+_b);
      _body.writeVTK(info._TM.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),path+"/lower"+std::to_string(i)+".vtk",Joint::MESH);
    }
    if(std::isfinite(_u[i])) {
      Vec x=Vec::Unit(_u.size(),i)*_u[i];
      PBDArticulatedGradientInfo<T> info(_body,_A*x+_b);
      _body.writeVTK(info._TM.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),path+"/upper"+std::to_string(i)+".vtk",Joint::MESH);
    }
  }
}
template <typename T>
typename GraspPlanner<T>::Vec GraspPlanner<T>::optimize(bool debug,const Vec& init,PointCloudObject<T>& object,GraspPlannerParameter& ops)
{
  _objs=DSSQPObjectiveCompound<T>();
  _info=PBDArticulatedGradientInfo<T>();
  if(ops._metric==Q_1 || ops._metric==Q_INF || ops._metric==Q_INF_BARRIER)
    _objs.addComponent(std::shared_ptr<ArticulatedObjective<T>>(new MetricEnergy<T>(_objs,_info,*this,object,ops._d0,ops._alpha,ops._coefM,(METRIC_TYPE)ops._metric,(METRIC_ACTIVATION)ops._activation,_rad*ops._normalExtrude)));
  if(ops._metric==Q_INF_CONSTRAINT)
    _objs.addComponent(std::shared_ptr<PrimalDualQInfMetricEnergy<T>>(new PrimalDualQInfMetricEnergy<T>(_objs,_info,*this,object,ops._alpha,ops._coefM,(METRIC_ACTIVATION)ops._activation,_rad*ops._normalExtrude)));
  if(ops._metric==Q_INF_CONSTRAINT_FGT)
    _objs.addComponent(std::shared_ptr<PrimalDualQInfMetricEnergyFGT<T>>(new PrimalDualQInfMetricEnergyFGT<T>(_objs,_info,*this,object,ops._alpha,ops._coefM,_rad*ops._normalExtrude,ops._FGTThres)));
  if(ops._coefOC>0)
    _objs.addComponent(std::shared_ptr<ArticulatedObjective<T>>(new ObjectClosednessEnergy<T>(_objs,_info,*this,object,ops._coefOC)));
  if(ops._coefCC>0)
    _objs.addComponent(std::shared_ptr<ArticulatedObjective<T>>(new CentroidClosednessEnergy<T>(_objs,_info,*this,object,ops._coefCC)));
  if(ops._coefO>0)
    _objs.addComponent(std::shared_ptr<ArticulatedObjective<T>>(new LogBarrierObjEnergy<T>(_objs,_info,*this,object,_rad*ops._d0,ops._coefO,ops._useGJK)));
  if(ops._coefS>0)
    _objs.addComponent(std::shared_ptr<ArticulatedObjective<T>>(new ConvexLogBarrierSelfEnergy<T>(_objs,_info,*this,object,_rad*ops._d0,ops._coefS)));

  Vec x;

  SolveNewton<T>::template solveNewton<Vec>(_A.transpose()*_A,_A.transpose()*(_b-init),x,true);
  sizeType nAdd=_objs.inputs()-init.size();
  std::cout << "nAdd = " << nAdd << std::endl;
  if(nAdd>0) {
    x=concat<Vec,Vec>(x,Vec::Zero(nAdd));
    _b=concat<Vec,Vec>(_b,Vec::Zero(nAdd));
    _A=concatDiag<T,0,sizeType>(_A,MatT::Identity(nAdd,nAdd).eval().sparseView());
    _l=concat<Vec,Vec>(_l,Vec::Constant(nAdd,-DSSQPObjective<T>::infty()));
    _u=concat<Vec,Vec>(_u,Vec::Constant(nAdd, DSSQPObjective<T>::infty()));
  }
  sizeType it;
  TBEG();
  if(debug)
    debugSystem(x);
  else x=optimizeSQP(x,ops,it);
  scalarD time=TENDV();
  INFOV("OptimizeSQP %d iterations, average time=%f",it,time/it)
  if(nAdd>0) {
    _b=_b.segment(0,_b.size()-nAdd).eval();
    _A=_A.block(0,0,_A.rows()-nAdd,_A.cols()-nAdd).eval();
    _l=_l.segment(0,_l.size()-nAdd).eval();
    _u=_u.segment(0,_u.size()-nAdd).eval();
  }
  return _A*x.segment(0,_A.cols())+_b;
}
template <typename T>
bool GraspPlanner<T>::solveDenseQP(Vec& d, const Vec& x,const Vec& g,MatT& h,const Vec* c,const MatT* cjac,T TR,T rho)
{
  scalarD maxConditionNumber=1e5f,minDiagonalValue=1e-5f;
  Eigen::SelfAdjointEigenSolver<Matd> eig(h.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  }),Eigen::ComputeEigenvectors);
  scalarD minEv=std::max<scalarD>(eig.eigenvalues().cwiseAbs().maxCoeff()/maxConditionNumber,minDiagonalValue);
  Cold ev=eig.eigenvalues().array().max(minEv).matrix();
  Matd hAdjusted=eig.eigenvectors()*ev.asDiagonal()*eig.eigenvectors().transpose();
  h=hAdjusted.template cast<T>();

  Vec lb=_l-x;
  Vec ub=_u-x;
  bool succ=false;
  if(c && c->size()>0) {
    //0.5*(x-x0)^T*H*(x-x0)+g^T*(x-x0)=
    //0.5*x^T*H*x-x0^T*H*x+0.5f*x0^T*H*x0+g^T*x-g^T*x0=
    //C+0.5*x^T*H*x-x0^T*H*x+g^T*x
    //
    //gl<=c+cjac*(x-x0)<=gu
    //gl-c+cjac*x0<=cjac*x<=gu-c+cjac*x0
    Vec gl=_gl-*c;
    Vec gu=_gu-*c;
    if(TR<=0)
      succ=_sol.solveQP(d=x,h,g,cjac,&lb,&ub,&gl,&gu,_objs.getQCones())==QCQPSolver<T>::SOLVED;
    else succ=_sol.solveL1QP(d=x,h,g,cjac,&lb,&ub,&gl,&gu,TR,rho,_objs.getQCones())==QCQPSolver<T>::SOLVED;
  } else {
    if(TR<=0)
      succ=_sol.solveQP(d=x,h,g,NULL,&lb,&ub,NULL,NULL,_objs.getQCones())==QCQPSolver<T>::SOLVED;
    else succ=_sol.solveL1QP(d=x,h,g,NULL,&lb,&ub,NULL,NULL,TR,rho,_objs.getQCones())==QCQPSolver<T>::SOLVED;
  }
  return succ;
}
template <typename T>
bool GraspPlanner<T>::solveSparseQP(Vec& d,const Vec& x,const Vec& g,SMat& h,const Vec* c,const SMat* cjac,T TR,T rho,T& reg)
{
  scalarD maxRegularization=1e5f,regInc=10.0f,regDec=0.9f;
  SMat Id=MatT::Identity(h.rows(),h.cols()).sparseView();

  Vec lb=_l-x;
  Vec ub=_u-x;
  bool succ=false;
  while(true) {
    SMat hReg=h+Id*reg;
    if(c && c->size()>0) {
      //0.5*(x-x0)^T*H*(x-x0)+g^T*(x-x0)=
      //0.5*x^T*H*x-x0^T*H*x+0.5f*x0^T*H*x0+g^T*x-g^T*x0=
      //C+0.5*x^T*H*x-x0^T*H*x+g^T*x
      //
      //gl<=c+cjac*(x-x0)<=gu
      //gl-c+cjac*x0<=cjac*x<=gu-c+cjac*x0
      Vec gl=_gl-*c;
      Vec gu=_gu-*c;
      if(TR<=0)
        succ=_sol.solveQP(d=x,hReg,g,cjac,&lb,&ub,&gl,&gu,_objs.getQCones())==QCQPSolver<T>::SOLVED;
      else succ=_sol.solveL1QP(d=x,hReg,g,cjac,&lb,&ub,&gl,&gu,TR,rho,_objs.getQCones())==QCQPSolver<T>::SOLVED;
    } else {
      if(TR<=0)
        succ=_sol.solveQP(d=x,hReg,g,NULL,&lb,&ub,NULL,NULL,_objs.getQCones())==QCQPSolver<T>::SOLVED;
      else succ=_sol.solveL1QP(d=x,hReg,g,NULL,&lb,&ub,NULL,NULL,TR,rho,_objs.getQCones())==QCQPSolver<T>::SOLVED;
    }
    //update regularization
    if(!succ) {
      reg=reg*regInc;
      if(reg>maxRegularization)
        return false;
    } else {
      reg=std::max<T>(1e-3f,reg*regDec);
      h=hReg;
      return true;
    }
  }
  return false;
}
template <typename T>
bool GraspPlanner<T>::assemble(Vec x,bool update,T& e,Vec* g,MatT* h,Vec* c,MatT* cjac)
{
  x=_A*x+_b;
  sizeType nCons=_objs.values();
  ParallelMatrix<Mat3XT> G;
  ParallelMatrix<Mat12XT> H;
  ParallelMatrix<T> E(0);
  if(g) {
    G.assign(Mat3XT::Zero(3,_body.nrJ()*4));
    g->setZero(x.size());
  }
  if(h) {
    H.assign(Mat12XT::Zero(12,_body.nrJ()*12));
    h->setZero(x.size(),x.size());
  }
  bool valid=true;
  for(typename std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>::const_iterator beg=_objs.components().begin(),end=_objs.components().end(); beg!=end; beg++) {
    beg->second->setUpdateCache(x,update);
    // std::cout << beg->second->_name << " " << std::dynamic_pointer_cast<ArticulatedObjective<T>>(beg->second)->operator()(x,E,g?&G:NULL,h?&H:NULL,g,h) << std::endl;
    // std::cout << "E value" << E.getValue() << std::endl;
    if(std::dynamic_pointer_cast<ArticulatedObjective<T>>(beg->second)->operator()(x,E,g?&G:NULL,h?&H:NULL,g,h)<0) {
      valid=false;
      // std::cout << beg->second->_name << " " << std::dynamic_pointer_cast<ArticulatedObjective<T>>(beg->second)->operator()(x,E,g?&G:NULL,h?&H:NULL,g,h) << std::endl;
      break;
    }
  }
  if(!valid)
    return false;
  //assemble body gradient / hessian

  Mat3XT tmpG;
  Mat12XT tmpH;
  e=E.getValue();
  if(g) {
    tmpG=G.getMatrix();
    _info.DTG(_body,mapM(tmpG),mapV(*g));
//   for (int i = 0; i < g->size(); i++)
//   {/* code */
//     std::cout << (*g)[i] << " ";
//   }
//   std::cout <<"In assemble"<< std::endl;
    *g=_A.transpose()**g;
    //  for (int i = 0; i < g->size(); i++)
    //   {/* code */
    //    std::cout << (*g)[i] << " ";
    //  }
    //  std::cout <<"After assemble"<< std::endl;
  }

  if(h) {
    tmpH=H.getMatrix();
    Eigen::Map<const MatT,0,Eigen::OuterStride<>> HMap(tmpH.data(),tmpH.rows(),tmpH.cols(),tmpH.outerStride());
    _info.toolAB(_body,HMap,mapM(tmpG=G.getMatrix()),mapM(*h));
    *h=_A.transpose()*(*h*_A);
  }
  //assemble constraint (jacobian)

  if(c || cjac) {
    if(c)
      c->setZero(nCons);
    if(cjac)
      cjac->setZero(nCons,x.size());
    if(nCons>0)
      for(typename std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>::const_iterator beg=_objs.components().begin(),end=_objs.components().end(); beg!=end; beg++)
        beg->second->setUpdateCache(x,update);
    if(_objs.DSSQPObjective<T>::operator()(x,*c,cjac)<0)
      return false;
    if(cjac)
      *cjac*=_A;
  }
  return true;
}
template <typename T>
bool GraspPlanner<T>::assemble(Vec x,bool update,T& e,Vec* g,SMat* h,Vec* c,SMat* cjac)
{
  x=_A*x+_b;
  sizeType nCons=_objs.values();
  ParallelMatrix<Mat3XT> G;
  ParallelMatrix<Mat12XT> H;
  ParallelMatrix<T> E(0);
  if(g) {
    G.assign(Mat3XT::Zero(3,_body.nrJ()*4));
    g->setZero(x.size());
  }
  if(h) {
    H.assign(Mat12XT::Zero(12,_body.nrJ()*12));
    h->resize(x.size(),x.size());
  }
  bool valid=true;
  for(typename std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>::const_iterator beg=_objs.components().begin(),end=_objs.components().end(); beg!=end; beg++) {
    beg->second->setUpdateCache(x,update);
    if(std::dynamic_pointer_cast<ArticulatedObjective<T>>(beg->second)->operator()(x,E,g?&G:NULL,h?&H:NULL,g,h)<0) {
      valid=false;
      break;
    }
  }
  if(!valid)
    return false;
  //assemble body gradient / hessian
  Mat3XT tmpG;
  Mat12XT tmpH;
  e=E.getValue();
  if(g) {
    tmpG=G.getMatrix();
    _info.DTG(_body,mapM(tmpG),mapV(*g));
    *g=_A.transpose()**g;
  }
  if(h) {
    MatT hDense;
    tmpH=H.getMatrix();
    hDense.setZero(x.size(),x.size());
    Eigen::Map<const MatT,0,Eigen::OuterStride<>> HMap(tmpH.data(),tmpH.rows(),tmpH.cols(),tmpH.outerStride());
    _info.toolAB(_body,HMap,mapM(tmpG=G.getMatrix()),mapM(hDense));
    *h=MatT(_A.transpose()*(hDense*_A)).sparseView();
  }
  //assemble constraint (jacobian)
  if(c || cjac) {
    if(c)
      c->setZero(nCons);
    if(cjac)
      cjac->resize(nCons,x.size());
    if(nCons>0)
      for(typename std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>::const_iterator beg=_objs.components().begin(),end=_objs.components().end(); beg!=end; beg++)
        beg->second->setUpdateCache(x,update);
    if(_objs.DSSQPObjective<T>::operator()(x,*c,cjac)<0)
      return false;
    if(cjac)
      *cjac=*cjac*_A;
  }
  return true;
}
template <typename T>
typename GraspPlanner<T>::Vec GraspPlanner<T>::optimizeSQP(Vec x,GraspPlannerParameter& ops,sizeType& it)
{
  Vec d;
  T e,e2,m,m2;
  MatT hD,cjacD;
  SMat hS,cjacS;
  Vec g,c,c2,xTmp, xTmpTmp;
  T dNorm,cNorm,cNorm2,alphaDec=0.5f,alphaInc=1.5f,coefWolfe=0.1f,alpha=1,rho=ops._rho0,gamma=0.1f,reg=0;
  _gl=_objs.gl(),_gu=_objs.gu();
  bool tmpUseGJK=ops._useGJK;

  for(it=0; it<ops._maxIter; it++) {
    if(ops._sparse) {
      if(!assemble(x,true,e,&g,&hS,&c,&cjacS)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(invalid configuration)",it)
        }
        return Vec::Zero(0);
      }
      if(reg==0)
        reg=std::max<T>(1e-3f,hS.toDense().diagonal().unaryExpr([&](const T& in) {
        return (scalarD)std::abs(in);
      }).maxCoeff());
      if(!solveSparseQP(d,x,g,hS,&c,&cjacS,0,0,reg)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(qp failed)",it)
        }
        break;
      }
    } else {
      if(!assemble(x,true,e,&g,&hD,&c,&cjacD)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(invalid configuration)",it)
        }
        return Vec::Zero(0);
      }
      if(!solveDenseQP(d,x,g,hD,&c,&cjacD,0,0)) {
        if(ops._callback) {
          INFOV("Iter=%d failed(qp failed)",it)
        }
        break;
      }
    }
    //termination & callback
    dNorm=std::sqrt(d.squaredNorm());
    cNorm=-c.cwiseMin(0).sum();
    if(dNorm<ops._thres && cNorm<ops._thres) {
      if(ops._callback) {
        INFOV("Iter=%d succeed(dNorm=%f<thres=%f,cNorm=%f<thres=%f)",it,std::to_double(dNorm),std::to_double(ops._thres),std::to_double(cNorm),std::to_double(ops._thres))
      }
      break;
    } else if(ops._callback) {
      INFOV("Iter=%d E=%f dNorm=%f cNorm=%f alpha=%f rho=%f",it,std::to_double(e),std::to_double(dNorm),std::to_double(cNorm),std::to_double(alpha),std::to_double(rho))
    }
    //merit-function parameter
    if(c.size()>0) {
      T D=d.dot(g);//+0.5f*d.dot(hAdjusted*d);
      if(D>0)
        rho=std::max(rho,D/((1-gamma)*cNorm));
      //replace g with directional derivative
      for(sizeType i=0; i<c.size(); i++)
        if(c[i]<0) {
          if(ops._sparse)
            g-=cjacS.row(i).transpose()*rho;
          else g-=cjacD.row(i).transpose()*rho;
        }
    }
    //line search
    m=e+cNorm*rho;
    ops._useGJK=true;
    while(alpha>ops._alphaThres) {
      xTmp=x+d*alpha;
      if(!assemble(xTmp,false,e2,(Vec*)NULL,(DMat*)NULL,&c2)) {
        alpha*=alphaDec;
        continue;
      }

      cNorm2=-c2.cwiseMin(0).sum();
      m2=e2+cNorm2*rho;
      if(m2>=m+g.dot(d.template cast<T>())*alpha*coefWolfe) {
        alpha*=alphaDec;
        continue;
      } else {
        alpha=std::min<T>(alpha*alphaInc,1);
        // alpha *= alphaInc;
        xTmpTmp=xTmp;
        break;
      }

    }
    ops._useGJK=tmpUseGJK;
    if(alpha<ops._alphaThres) {
      if(ops._callback) {
        INFOV("Iter=%d failed(alpha=%f<alphaThres=%f)",it,std::to_double(alpha),std::to_double(ops._alphaThres))
      }
      break;
    }
    if(!assemble(xTmpTmp,false,e2,(Vec*)NULL,(DMat*)NULL,&c2)) {
      std::cout << "Cannot use GJK.\n Redoing the line search..." << std::endl;
      while(alpha>ops._alphaThres) {
        xTmp=x+d*alpha;
        if(!assemble(xTmp,false,e2,(Vec*)NULL,(DMat*)NULL,&c2)) {
          alpha*=alphaDec;
          continue;
        }

        cNorm2=-c2.cwiseMin(0).sum();
        m2=e2+cNorm2*rho;
        if(m2>=m+g.dot(d.template cast<T>())*alpha*coefWolfe) {
          alpha*=alphaDec;
          continue;
        } else {
          alpha=std::min<T>(alpha*alphaInc,1);
          // alpha *= alphaInc;
          x=xTmp;
          break;
        }

      }
    }
    else x=xTmpTmp;
    //update plane
    for(const std::pair<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>& p:_objs.components()) {
      std::shared_ptr<ConvexLogBarrierSelfEnergy<T>> ESelf=std::dynamic_pointer_cast<ConvexLogBarrierSelfEnergy<T>>(p.second);
      if(ESelf) {
        ESelf->updatePlanes();
      }
    }
    if(!assemble(x,false,e2,(Vec*)NULL,(DMat*)NULL,&c2)) {
      std::cout << "after updating plane goes wrong" << std::endl;
    }
  }
  return x;
}
template <typename T>
void GraspPlanner<T>::evaluateQInf( Vec& x, PointCloudObject<T>& object,GraspPlannerParameter& ops)
{
  _objs=DSSQPObjectiveCompound<T>();
  _info=PBDArticulatedGradientInfo<T>();
  ParallelMatrix<T> E(0);
  _objs.addComponent(std::shared_ptr<ArticulatedObjective<T>>(new MetricEnergy<T>(_objs,_info,*this,object,ops._d0,ops._alpha,ops._coefM,(METRIC_TYPE)ops._metric,(METRIC_ACTIVATION)ops._activation,_rad*ops._normalExtrude)));
  sizeType nAdd=_objs.inputs()-x.size();
  if(nAdd>0) {
    x=concat<Vec,Vec>(x,Vec::Zero(nAdd));
    _b=concat<Vec,Vec>(_b,Vec::Zero(nAdd));
    _A=concatDiag<T,0,sizeType>(_A,MatT::Identity(nAdd,nAdd).eval().sparseView());
    _l=concat<Vec,Vec>(_l,Vec::Constant(nAdd,-DSSQPObjective<T>::infty()));
    _u=concat<Vec,Vec>(_u,Vec::Constant(nAdd, DSSQPObjective<T>::infty()));
  }
  x=_A*x+_b;
  for(typename std::unordered_map<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>::const_iterator beg=_objs.components().begin(),end=_objs.components().end(); beg!=end; beg++) {
    beg->second->setUpdateCache(x,true);
    std::dynamic_pointer_cast<ArticulatedObjective<T>>(beg->second)->operator()(x,E,NULL,NULL,(Vec*)NULL,(DMat*)NULL);

//    std::cout << beg->second->_name << " " << std::dynamic_pointer_cast<ArticulatedObjective<T>>(beg->second)->Quality(x)<< std::endl;
  }


}
template <typename T>
void GraspPlanner<T>::debugSystem(const Vec& x)
{
  DEFINE_NUMERIC_DELTA_T(T)
  Vec dx=Vec::Random(x.size());
  //first evaluate
  T e;
  Vec g,c;
  MatT hD,cjacD;
  SMat hS,cjacS;
  assemble(x,true,e,&g,&hD,&c,&cjacD);
  assemble(x,true,e,&g,&hS,&c,&cjacS);
  //second evaluate
  T e2;
  Vec g2,c2;
  assemble(x+dx*DELTA,false,e2,&g2,(DMat*)NULL,&c2,NULL);
  //compare
  DEBUG_GRADIENT("GraspPlanner-G",g.dot(dx),g.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT("GraspPlanner-HDense",std::sqrt((hD*dx).squaredNorm()),std::sqrt((hD*dx-(g2-g)/DELTA).squaredNorm()))
  DEBUG_GRADIENT("GraspPlanner-HSparse",std::sqrt((hS*dx).squaredNorm()),std::sqrt((hS*dx-(g2-g)/DELTA).squaredNorm()))
  if(_objs.values()>0) {
    DEBUG_GRADIENT("GraspPlanner-CJacDense",std::sqrt((cjacD*dx).squaredNorm()),std::sqrt((cjacD*dx-(c2-c)/DELTA).squaredNorm()))
    DEBUG_GRADIENT("GraspPlanner-CJacSparse",std::sqrt((cjacS*dx).squaredNorm()),std::sqrt((cjacS*dx-(c2-c)/DELTA).squaredNorm()))
  }
}
template <typename T>
const typename GraspPlanner<T>::SMat& GraspPlanner<T>::A() const
{
  return _A;
}
template <typename T>
const typename GraspPlanner<T>::Vec& GraspPlanner<T>::b() const
{
  return _b;
}
template <typename T>
T GraspPlanner<T>::area() const
{
  return _rad*_rad*M_PI;
}
template <typename T>
T GraspPlanner<T>::rad() const
{
  return _rad;
}
template <typename T>
bool GraspPlanner<T>::validSample(sizeType l,const PBDArticulatedGradientInfo<T>& info,const Vec3T& p) const
{
  sizeType j;
  for(sizeType i=0; i<_body.nrJ(); i++) {
    ASSERT_MSG(l>=0 && l<_body.nrJ(),"Invalid joint id")
    if(i==l)
      continue;
    bool isDirectParent=false;
    //find from i to l
    j=_body.joint(i)._parent;
    while(j>=0) {
      if(_body.joint(j)._M>0)
        break;
      j=_body.joint(j)._parent;
    }
    if(j==l)
      isDirectParent=true;
    //find from l to i
    j=_body.joint(l)._parent;
    while(j>=0) {
      if(_body.joint(j)._M>0)
        break;
      j=_body.joint(j)._parent;
    }
    if(j==i)
      isDirectParent=true;
    //check
    if(isDirectParent)
      if(!env(i).empty() && env(i).phi(ROTI(info._TM,i).transpose()*(p-CTRI(info._TM,i)))<=0)
        return false;
  }
  return true;
}
//instance
PRJ_BEGIN
template class GraspPlanner<double>;
#ifdef ALL_TYPES
template class GraspPlanner<__float128>;
template class GraspPlanner<mpfr::mpreal>;
#endif
PRJ_END
