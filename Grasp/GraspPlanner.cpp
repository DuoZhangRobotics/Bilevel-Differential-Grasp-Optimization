#include "GraspPlanner.h"
#include <Utils/Utils.h>
#include <Utils/DebugGradient.h>
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/MultiPrecisionLQP.h>
#include <CommonFile/ParallelPoissonDiskSampling.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <Environment/ConvexHullExact.h>
#include <Eigen/Eigen>
#include <qpOASES.hpp>
//energy
#include "PrimalDualQInfMetricEnergyFGT.h"
#include "PrimalDualQInfMetricEnergy.h"
#include "CentroidClosednessEnergy.h"
#include "ObjectClosednessEnergy.h"
#include "LogBarrierObjEnergy.h"
#include "LogBarrierSelfEnergy.h"
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
  sol._coefM=-100;
  sol._coefOC=0;
  sol._coefCC=0;
  sol._coefO=10;
  sol._coefS=1;
  sol._useGJK=false;
  //solver
  sol._rho0=1;
  sol._thres=1e-10f;
  sol._alphaThres=1e-20f;
  sol._callback=true;
  sol._maxIter=2000;
}
//GraspPlanner
template <typename T>
GraspPlanner<T>::GraspPlanner() {}
template <typename T>
void GraspPlanner<T>::reset(const std::string& path,T rad,bool convex)
{
  _body=ArticulatedLoader().readURDF(path,convex,true);
  std::cout <<"BODY READ: " << _body.nrDOF() <<std::endl;
  ArticulatedUtils(_body).addBase(3,Vec3d::Zero());
  ArticulatedUtils(_body).simplify(10);
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
  _distExact.resize(_body.nrJ());
  for(sizeType i=0; i<_body.nrJ(); i++)
    if(convex)
      _distExact[i].reset(new ConvexHullExact(dynamic_cast<const ObjMeshGeomCell&>(_body.getGeom().getG(i))));
    else _distExact[i].reset(new ObjMeshGeomCellExact(dynamic_cast<const ObjMeshGeomCell&>(_body.getGeom().getG(i))));
  //sample

  PBDArticulatedGradientInfo<T> info(_body,Vec::Zero(_body.nrDOF()));
  for(sizeType i=0; i<_body.nrJ(); i++) {
    _body.getGeom().getG(i).getMesh(m);
    if(m.getV().empty())
      continue;
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
      if(!validSample(i,info,ROTI(info._TM,i)*v+CTRI(info._TM,i)))
        continue;
      _pnss[i].first.col(k)=v;
      _pnss[i].second.col(k)=n;
      k++;
    }
    _pnss[i].first=_pnss[i].first.block(0,0,3,k).eval();
    _pnss[i].second=_pnss[i].second.block(0,0,3,k).eval();
  }
  std::cout << "Done" << std::endl;

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
  registerType<ConvexHullExact>(dat);
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
  registerType<ConvexHullExact>(dat);
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
typename GraspPlanner<T>::Vec GraspPlanner<T>::optimize(bool debug,const Vec& init,GraspQualityMetric<T>& obj,GraspPlannerParameter& ops)
{
  std::vector<std::shared_ptr<ArticulatedObjective<T>>> objs;
  if(ops._metric==Q_1 || ops._metric==Q_INF || ops._metric==Q_INF_BARRIER)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new MetricEnergy<T>(*this,obj,ops._d0,ops._alpha,ops._coefM,(METRIC_TYPE)ops._metric,(METRIC_ACTIVATION)ops._activation,_rad*ops._normalExtrude)));
  if(ops._metric==Q_INF_CONSTRAINT)
    objs.push_back(std::shared_ptr<PrimalDualQInfMetricEnergy<T>>(new PrimalDualQInfMetricEnergy<T>(*this,obj,ops._alpha,ops._coefM,(METRIC_ACTIVATION)ops._activation,_rad*ops._normalExtrude)));
  if(ops._metric==Q_INF_CONSTRAINT_FGT)
    objs.push_back(std::shared_ptr<PrimalDualQInfMetricEnergyFGT<T>>(new PrimalDualQInfMetricEnergyFGT<T>(*this,obj,ops._alpha,ops._coefM,_rad*ops._normalExtrude,ops._FGTThres)));
  if(ops._coefOC>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new ObjectClosednessEnergy<T>(*this,obj,ops._coefOC)));
  if(ops._coefCC>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new CentroidClosednessEnergy<T>(*this,obj,ops._coefCC)));
  if(ops._coefO>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new LogBarrierObjEnergy<T>(*this,obj,_rad*ops._d0,ops._coefO,ops._useGJK)));
  if(ops._coefS>0)
    objs.push_back(std::shared_ptr<ArticulatedObjective<T>>(new LogBarrierSelfEnergy<T>(*this,obj,_rad*ops._d0,ops._coefS)));

  Vec x;
  SolveNewton<T>::template solveNewton<Vec>(_A.transpose()*_A,_A.transpose()*(_b-init),x,true);
  sizeType nAdd=nrAdditionalDOF(objs);
  if(nAdd>0) {
    x=concat<Vec,Vec>(x,Vec::Zero(nAdd));
    _b=concat<Vec,Vec>(_b,Vec::Zero(nAdd));
    _A=concatDiag<T,0,sizeType>(_A.sparseView(),MatT::Identity(nAdd,nAdd).eval().sparseView());
    _l=concat<Vec,Vec>(_l,Vec::Constant(nAdd,-qpOASES::INFTY));
    _u=concat<Vec,Vec>(_u,Vec::Constant(nAdd,qpOASES::INFTY));
  }
  if(debug)
    debugSystem(x,objs);
  else x=optimizeNewton(x,objs,ops);
  if(nAdd>0) {
    _b=_b.segment(0,_b.size()-nAdd).eval();
    _A=_A.block(0,0,_A.rows()-nAdd,_A.cols()-nAdd).eval();
    _l=_l.segment(0,_l.size()-nAdd).eval();
    _u=_u.segment(0,_u.size()-nAdd).eval();
  }
  return _A*x.segment(0,_A.cols())+_b;
}
template <typename T>
typename GraspPlanner<T>::Vec GraspPlanner<T>::optimizeNewton(Vec x,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs,GraspPlannerParameter& ops) const
{
  Cold d;
  T e,e2,m,m2;
  MatT h,cjac;
  Vec g,c,c2,xTmp;
  scalarD maxConditionNumber=1e5f,minDiagonalValue=1e-5f;
  Eigen::Matrix<scalarD,-1,-1,Eigen::RowMajor> hAdjusted;
  T dNorm,cNorm,cNorm2,alphaDec=0.5f,alphaInc=1.5f,coefWolfe=0.1f,alpha=1,rho=ops._rho0,gamma=0.1f;
  PBDArticulatedGradientInfo<T> info;
  bool tmpUseGJK=ops._useGJK;

  for(sizeType it=0; it<ops._maxIter; it++) {
    //assemble
    if(!assemble(x,info,true,objs,e,&g,&h,&c,&cjac)) {
      if(ops._callback) {
        INFOV("Iter=%d failed(invalid configuration)",it)
      }
      return Vec::Zero(0);
    }
    //solve
    {
      Cold cld;
      Eigen::Matrix<scalarD,-1,-1,Eigen::RowMajor> cjacd;
      Eigen::SelfAdjointEigenSolver<Matd> eig(h.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),Eigen::ComputeEigenvectors);
      scalarD minEv=std::max<scalarD>(eig.eigenvalues().cwiseAbs().maxCoeff()/maxConditionNumber,minDiagonalValue);
      Cold ev=eig.eigenvalues().array().max(minEv).matrix();
      hAdjusted=eig.eigenvectors()*ev.asDiagonal()*eig.eigenvectors().transpose();
      if(c.size()>0) {
        cjacd=cjac.unaryExpr([&](const T& in) {
          return (scalarD)std::to_double(in);
        });
        cld=-c.unaryExpr([&](const T& in) {
          return (scalarD)std::to_double(in);
        });
      }
      //respect joint limits
      Cold ld=(_l-x).unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),ud=(_u-x).unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),gAdjusted=g.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      });
      qpOASES::int_t nWSR=10000;
      qpOASES::SQProblem prob(g.size(),c.size());
      if(c.size()>0)
        prob.init(hAdjusted.data(),gAdjusted.data(),cjacd.data(),ld.data(),ud.data(),cld.data(),NULL,nWSR);
      else prob.init(hAdjusted.data(),gAdjusted.data(),NULL,ld.data(),ud.data(),NULL,NULL,nWSR);
      if(!prob.isSolved()) {
        if(ops._callback) {
          INFOV("Iter=%d failed(qpOASES failed)",it)
        }
        break;
      } else {
        d.resize(x.size());
        prob.getPrimalSolution(d.data());
      }
    }
    //termination & callback
    dNorm=d.norm();
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
      T D=d.template cast<T>().dot(g);//+0.5f*d.dot(hAdjusted*d);
      if(D>0)
        rho=std::max(rho,D/((1-gamma)*cNorm));
      //replace g with directional derivative
      for(sizeType i=0; i<c.size(); i++)
        if(c[i]<0)
          g-=cjac.row(i).transpose()*rho;
    }
    //line search
    m=e+cNorm*rho;
    ops._useGJK=true;
    while(alpha>ops._alphaThres) {
      xTmp=x+d.template cast<T>()*alpha;
      if(!assemble(xTmp,info,false,objs,e2,NULL,NULL,&c2)) {
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
        x=xTmp;
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
    //update plane
    for(sizeType i=0; i<(sizeType)objs.size(); i++) {
      LogBarrierSelfEnergy<T>* ESelf=dynamic_cast<LogBarrierSelfEnergy<T>*>(objs[i].get());
      if(ESelf)
        ESelf->updatePlanes(info);
    }
  }
  return x;
}
template <typename T>
bool GraspPlanner<T>::assemble(Vec x,PBDArticulatedGradientInfo<T>& info,bool update,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs,T& e,Vec* g,MatT* h,Vec* c,MatT* cjac) const
{
  x=_A*x+_b;
  info.reset(_body,x);
  sizeType nCons=nrCons(objs),nDOF=_body.nrDOF();
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
  for(sizeType i=0,offC=nDOF; i<(sizeType)objs.size() && valid; i++) {
    objs[i]->setUpdateCache(update);
    if(objs[i]->operator()(x,info,offC,E,g?&G:NULL,h?&H:NULL,g,h)<0) {
      valid=false;
      break;
    }
    offC+=objs[i]->nrAdditionalDOF();
  }
  if(!valid)
    return false;
  //assemble body gradient / hessian
  Mat3XT tmpG;
  Mat12XT tmpH;
  e=E.getValue();
  if(g) {
    tmpG=G.getMatrix();
    info.DTG(_body,mapM(tmpG),mapV(*g));
    *g=_A.transpose()**g;
  }
  if(h) {
    tmpH=H.getMatrix();
    Eigen::Map<const MatT,0,Eigen::OuterStride<>> HMap(tmpH.data(),tmpH.rows(),tmpH.cols(),tmpH.outerStride());
    info.toolAB(_body,HMap,mapM(tmpG=G.getMatrix()),mapM(*h));
    *h=_A.transpose()*(*h*_A);
  }
  //assemble constraint (jacobian)
  if(c || cjac) {
    if(c)
      c->setZero(nCons);
    if(cjac)
      cjac->setZero(nCons,x.size());
    if(nCons>0)
      for(sizeType i=0,offR=0,offC=nDOF; i<(sizeType)objs.size(); i++) {
        objs[i]->setUpdateCache(update);
        objs[i]->cons(x,info,offR,offC,*c,cjac);
        offC+=objs[i]->nrAdditionalDOF();
        offR+=objs[i]->nrCons();
      }
    if(cjac)
      *cjac*=_A;
  }
  return true;
}
template <typename T>
void GraspPlanner<T>::debugSystem(const Vec& x,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs) const
{
  DEFINE_NUMERIC_DELTA_T(T)
  Vec dx=Vec::Random(x.size());
  PBDArticulatedGradientInfo<T> info;
  //first evaluate
  T e;
  Vec g,c;
  MatT h,cjac;
  assemble(x,info,true,objs,e,&g,&h,&c,&cjac);
  //second evaluate
  T e2;
  Vec g2,c2;
  assemble(x+dx*DELTA,info,false,objs,e2,&g2,NULL,&c2,NULL);
  //compare
  DEBUG_GRADIENT("GraspPlanner-G",g.dot(dx),g.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT("GraspPlanner-H",std::sqrt((h*dx).squaredNorm()),std::sqrt((h*dx-(g2-g)/DELTA).squaredNorm()))
  if(nrCons(objs)>0) {
    DEBUG_GRADIENT("GraspPlanner-CJac",std::sqrt((cjac*dx).squaredNorm()),std::sqrt((cjac*dx-(c2-c)/DELTA).squaredNorm()))
  }
}
template <typename T>
sizeType GraspPlanner<T>::nrAdditionalDOF(std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs) const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)objs.size(); i++)
    ret+=objs[i]->nrAdditionalDOF();
  return ret;
}
template <typename T>
sizeType GraspPlanner<T>::nrCons(std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs) const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)objs.size(); i++)
    ret+=objs[i]->nrCons();
  return ret;
}
template <typename T>
const typename GraspPlanner<T>::MatT& GraspPlanner<T>::A() const
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
  Vec3T n,normal;
  Mat3T hessian;
  Vec2i feat;
  for(sizeType i=0; i<_body.nrJ(); i++) {
    ASSERT_MSG(l>=0 && l<_body.nrJ(),"Invalid joint id")
    if(_body.joint(i)._parent!=l && _body.joint(l)._parent!=i)
      continue;
    Vec3T pi=ROTI(info._TM,i).transpose()*(p-CTRI(info._TM,i));
    if(!_distExact[i]->empty() && _distExact[i]->template closest<T>(pi,n,normal,hessian,feat)<=0)
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
