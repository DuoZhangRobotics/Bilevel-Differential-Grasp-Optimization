#include "GraspPlanner.h"
#include <Utils/Utils.h>
#include <Utils/DebugGradient.h>
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/MultiPrecisionLQP.h>
#include <CommonFile/geom/ObjMeshGeomCell.h>
#include <CommonFile/ParallelPoissonDiskSampling.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <Utils/ArticulatedBodyPragma.h>
#include "CentroidClosednessEnergy.h"
#include "ObjectClosednessEnergy.h"
#include "LogBarrierObjEnergy.h"
#include "LogBarrierSelfEnergy.h"
#include "MetricEnergy.h"
#include <Eigen/Eigen>
#include <qpOASES.hpp>

USE_PRJ_NAMESPACE

//GraspPlanner
template <typename T>
GraspPlanner<T>::GraspPlanner() {}
template <typename T>
void GraspPlanner<T>::reset(const std::string& path,T rad)
{
  _body=ArticulatedLoader().readURDF(path,true,true);
  ArticulatedUtils(_body).addBase(3,Vec3d::Zero());
  ArticulatedUtils(_body).simplify(10);
  _distExact.resize(_body.nrJ());
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
      l[i]=-qpOASES::INFTY;
    if(!std::isfinite(u[i]))
      u[i]= qpOASES::INFTY;
  }
  _A=A.template cast<T>();
  _b=b.template cast<T>();
  _l=l.template cast<T>();
  _u=u.template cast<T>();
  //distExact
  for(sizeType i=0; i<(sizeType)_distExact.size(); i++) {
    _distExact[i].reset(new ObjMeshGeomCellExact(dynamic_cast<const ObjMeshGeomCell&>(_body.getGeom().getG(i))));
    //sample
    _body.getGeom().getG(i).getMesh(m);
    if(m.getV().empty())
      continue;
    ParallelPoissonDiskSampling sampler(3);
    sampler.setRadius(std::to_double(_rad));
    sampler.sample(m);
    //pss
    _pnss[i].first.resize(3,sampler.getPSet().size());
    _pnss[i].second.resize(3,sampler.getPSet().size());
    for(sizeType j=0; j<sampler.getPSet().size(); j++) {
      _pnss[i].first.col(j)=sampler.getPSet()[j]._pos.template cast<T>();
      _pnss[i].second.col(j)=sampler.getPSet()[j]._normal.template cast<T>();
    }
  }
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
  registerType<ObjMeshGeomCellExact>(dat);
  registerType<GraspPlanner<T>>(dat);
  readBinaryData(_distExact,is,dat);
  _body.read(is,dat);
  //mimic
  Matd A;
  readBinaryData(A,is);
  _A=A.template cast<T>();
  Cold b;
  readBinaryData(b,is);
  _b=b.template cast<T>();
  Cold l;
  readBinaryData(l,is);
  _l=l.template cast<T>();
  Cold u;
  readBinaryData(u,is);
  _u=u.template cast<T>();
  //sample
  std::vector<std::pair<Mat3Xd,Mat3Xd>> pnss;
  readBinaryData(pnss,is);
  _pnss.resize(pnss.size());
  for(sizeType i=0; i<(sizeType)pnss.size(); i++) {
    _pnss[i].first=pnss[i].first.template cast<T>();
    _pnss[i].second=pnss[i].second.template cast<T>();
  }
  scalarD rad;
  readBinaryData(rad,is);
  _rad=rad;
  return is.good();
}
template <typename T>
bool GraspPlanner<T>::write(std::ostream& os,IOData* dat) const
{
  registerType<ObjMeshGeomCellExact>(dat);
  registerType<GraspPlanner<T>>(dat);
  writeBinaryData(_distExact,os,dat);
  _body.write(os,dat);
  //mimic
  Matd A=_A.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(A,os);
  Cold b=_b.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(b,os);
  Cold l=_l.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(l,os);
  Cold u=_u.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(u,os);
  //sample
  std::vector<std::pair<Mat3Xd,Mat3Xd>> pnss(_pnss.size());
  for(sizeType i=0; i<(sizeType)pnss.size(); i++) {
    pnss[i].first=_pnss[i].first.unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    });
    pnss[i].second=_pnss[i].second.unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    });
  }
  writeBinaryData(pnss,os);
  scalarD rad=std::to_double(_rad);
  writeBinaryData(rad,os);
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
const ObjMeshGeomCellExact& GraspPlanner<T>::dist(sizeType jid) const
{
  return *(_distExact.at(jid));
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
      _distExact[j]->getMesh().writeVTK(path+"/jointMesh"+std::to_string(j)+".vtk",true);
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
typename GraspPlanner<T>::Vec GraspPlanner<T>::optimize(const Vec& init,GraspQualityMetric<T>& obj,T d0,T alpha,METRIC_TYPE m,T coefM,T coefOC,T coefCC,T coefO,T coefS) const
{
  std::vector<std::shared_ptr<ArticulatedObjective<T>>> objs;
  if(m!=NO_METRIC)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new MetricEnergy<T>(*this,obj,alpha,m,coefM)));
  if(coefOC>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new ObjectClosednessEnergy<T>(*this,obj,coefOC)));
  if(coefCC>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new CentroidClosednessEnergy<T>(*this,obj,coefCC)));
  if(coefO>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new LogBarrierObjEnergy<T>(*this,obj,_rad*d0,coefO)));
  if(coefS>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new LogBarrierSelfEnergy<T>(*this,obj,_rad*d0,coefS)));
  Vec x;
  SolveNewton<T>::template solveNewton<Vec>(_A.transpose()*_A,_A.transpose()*(_b-init),x,true);
  x=optimizeNewton(x,objs);
  return _A*x+_b;
}
template <typename T>
typename GraspPlanner<T>::Vec GraspPlanner<T>::optimizeNewton(Vec x,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs,sizeType maxIter,T gThres,T alphaThres,bool callback) const
{
  x=x.cwiseMin(_u).cwiseMax(_l);
  T E,E2,gNorm,alphaDec=0.5f,alphaInc=1.5f,coefWolfe=0.1f,alpha=1;
  scalarD maxConditionNumber=1e5f;
  Mat3XT G,tmpG;
  MatT h,H,tmpH;
  Vec g,gRes,xTmp;
  Cold d;
  for(sizeType it=0; it<maxIter; it++) {
    //assemble
    {
      PBDArticulatedGradientInfo<T> info(_body,_A*x+_b);
      bool valid=true;
      E=0;
      G.setZero(3,_body.nrJ()*4);
      H.setZero(12,_body.nrJ()*12);
      for(sizeType i=0; i<(sizeType)objs.size() && valid; i++)
        if(objs[i]->operator()(info,E,&G,&H)<0) {
          valid=false;
          break;
        }
      if(!valid) {
        if(callback) {
          INFOV("Iter=%d failed(invalid configuration)",it)
        }
        return Vec::Zero(0);
      }
      g.setZero(_body.nrDOF());
      h.setZero(_body.nrDOF(),_body.nrDOF());
      info.DTG(_body,mapM(tmpG=G),mapV(g));
      info.toolAB(_body,mapCM(tmpH=H),mapM(tmpG=G),mapM(h));
      //handle mimic
      g=_A.transpose()*g;
      h=_A.transpose()*(h*_A);
    }
    //termination & callback
    gRes=g; //account for joint limits
    for(sizeType i=0; i<gRes.size(); i++)
      if(x[i]<=_l[i]) {
        if(gRes[i]>0)
          gRes[i]=0;
      } else if(x[i]>=_u[i]) {
        if(gRes[i]<0)
          gRes[i]=0;
      }
    gNorm=std::sqrt(gRes.squaredNorm());
    if(gNorm<gThres) {
      if(callback) {
        INFOV("Iter=%d succeed(gNorm=%f<gThres=%f)",it,std::to_double(gNorm),std::to_double(gThres))
      }
      break;
    } else if(callback) {
      //we print the base for debug
      INFOV("Iter=%d E=%f gNorm=%f base=%f,%f,%f",
            it,std::to_double(E),std::to_double(gNorm),
            std::to_double(x[0]),std::to_double(x[1]),std::to_double(x[2]))
    }
    //solve
    {
      Eigen::SelfAdjointEigenSolver<Matd> eig(h.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),Eigen::ComputeEigenvectors);
      scalarD minEv=eig.eigenvalues().cwiseAbs().maxCoeff()/maxConditionNumber;
      Cold ev=eig.eigenvalues().array().max(minEv).matrix();
      Matd hAdjusted=eig.eigenvectors()*ev.asDiagonal()*eig.eigenvectors().transpose();
      //respect joint limits
      Cold ld=(_l-x).unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),ud=(_u-x).unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),gAdjusted=g.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      });
      qpOASES::int_t nWSR=10000;
      qpOASES::SQProblem prob(g.size(),0);
      prob.init(hAdjusted.data(),gAdjusted.data(),NULL,ld.data(),ud.data(),NULL,NULL,nWSR);
      if(!prob.isSolved()) {
        if(callback) {
          INFOV("Iter=%d failed(qpOASES failed)",it)
        }
      } else {
        d.resize(x.size());
        prob.getPrimalSolution(d.data());
        //T moveBase=std::sqrt(d.segment<3>(0).squaredNorm());
        //if(moveBase>_rad)
        //  d*=std::to_double(_rad/moveBase);
      }
    }
    //line search
    while(alpha>alphaThres) {
      xTmp=x+d.template cast<T>()*alpha;
      PBDArticulatedGradientInfo<T> info(_body,_A*xTmp+_b);
      bool valid=true;
      E2=0;
      G.setZero(3,_body.nrJ()*4);
      for(sizeType i=0; i<(sizeType)objs.size() && valid; i++)
        if(objs[i]->operator()(info,E2,NULL,NULL)<0) {
          valid=false;
          break;
        }
      if(!valid || E2>=E+g.dot(d.template cast<T>())*alpha*coefWolfe) {
        alpha*=alphaDec;
        continue;
      } else {
        alpha=std::min<T>(alpha*alphaInc,1);
        x=xTmp;
        break;
      }
    }
    if(alpha<alphaThres) {
      if(callback) {
        INFOV("Iter=%d failed(alpha=%f<alphaThres=%f)",it,std::to_double(alpha),std::to_double(alphaThres))
      }
      break;
    }
    //update plane
    {
      PBDArticulatedGradientInfo<T> info(_body,_A*x+_b);
      for(sizeType i=0; i<(sizeType)objs.size(); i++) {
        LogBarrierSelfEnergy<T>* ESelf=dynamic_cast<LogBarrierSelfEnergy<T>*>(objs[i].get());
        if(ESelf)
          ESelf->updatePlanes(info);
      }
    }
  }
  return x;
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
//ArticulatedObjective
template <typename T>
ArticulatedObjective<T>::ArticulatedObjective(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj):_planner(planner),_obj(obj) {}
template <typename T>
int ArticulatedObjective<T>::operator()(const Vec& x,T& e,Vec* g,MatT* h)
{
  PBDArticulatedGradientInfo<T> info(_planner.body(),x.segment(0,_planner.body().nrDOF()));
  Mat3XT G,tmpG;
  MatT H,tmpH;
  if(g)
    G.setZero(3,_planner.body().nrJ()*4);
  if(h)
    H.setZero(12,_planner.body().nrJ()*12);
  int ret=operator()(info,e,g?&G:NULL,h?&H:NULL);
  if(ret<0)
    return ret;
  if(g)
    info.DTG(_planner.body(),mapM(tmpG=G),mapV(*g));
  if(h)
    info.toolAB(_planner.body(),mapCM(tmpH=H),mapM(tmpG=G),mapM(h));
  return ret;
}
template <typename T>
void ArticulatedObjective<T>::debug(const Vec& x)
{
  sizeType nDOF=_planner.body().nrDOF();
  DEFINE_NUMERIC_DELTA_T(T)
  T e=0,e2=0;
  MatT h=MatT::Zero(nDOF,nDOF);
  Vec g=Vec::Zero(nDOF),g2=Vec::Zero(nDOF);
  Vec dx=Vec::Random(nDOF);
  if(operator()(x,e,&g,&h)<0) {
    std::ostringstream oss;
    for(sizeType i=0; i<x.size(); i++)
      oss << x[i] << " ";
    WARNINGV("Invalid debug position (x=%s)",oss.str().c_str())
    return;
  }
  operator()(x+dx*DELTA,e2,&g2);
  DEBUG_GRADIENT("ArticulatedObjective-G",g.dot(dx),g.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT("ArticulatedObjective-H",std::sqrt((h*dx).squaredNorm()),std::sqrt((h*dx-(g2-g)/DELTA).squaredNorm()))
}
template <typename T>
std::vector<KDOP18<scalar>> ArticulatedObjective<T>::updateBVH(const PBDArticulatedGradientInfo<T>& info) const
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  std::vector<KDOP18<scalar>> ret(bvhHand.size());
  for(sizeType i=0; i<(sizeType)bvhHand.size(); i++)
    if(bvhHand[i]._cell) {
      Vec3 pos;
      BBox<scalar> bb;
      Mat3 R=ROTI(info._TM,i).unaryExpr([&](const T& in) {
        return (scalar)std::to_double(in);
      });
      Vec3 t=CTRI(info._TM,i).unaryExpr([&](const T& in) {
        return (scalar)std::to_double(in);
      });
      for(sizeType x=0; x<2; x++) {
        pos[0]=x==0?bvhHand[i]._bb._minC[0]:bvhHand[i]._bb._maxC[0];
        for(sizeType y=0; y<2; y++) {
          pos[1]=y==0?bvhHand[i]._bb._minC[1]:bvhHand[i]._bb._maxC[1];
          for(sizeType z=0; z<2; z++) {
            pos[2]=z==0?bvhHand[i]._bb._minC[2]:bvhHand[i]._bb._maxC[2];
            ret[i].setUnion(R*pos+t);
          }
        }
      }
    } else {
      ret[i]=ret[bvhHand[i]._l];
      ret[i].setUnion(ret[bvhHand[i]._r]);
    }
  return ret;
}
//instance
PRJ_BEGIN
template class GraspPlanner<double>;
template class ArticulatedObjective<double>;
#ifdef ALL_TYPES
template class GraspPlanner<__float128>;
template class ArticulatedObjective<__float128>;
template class GraspPlanner<mpfr::mpreal>;
template class ArticulatedObjective<mpfr::mpreal>;
#endif
PRJ_END
