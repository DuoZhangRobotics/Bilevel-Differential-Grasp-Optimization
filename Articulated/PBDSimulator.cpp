#ifdef ENVIRONMENT_SUPPORT
#include "PBDSimulator.h"
#include "MultiPrecisionLQP.h"
#include <Optimizer/KNInterface.h>
#include <Optimizer/IPOPTInterface.h>
#include <CommonFile/Timing.h>
#include <Utils/DebugGradient.h>
#include <Utils/RotationUtil.h>
#include <Utils/Utils.h>
#include "PDTarget.h"

PRJ_BEGIN

//parameters
void PBDSimulatorTraits<double>::registerOptions(Options& ops)
{
  REGISTER_FLOAT_TYPE("alphaPBTO",PBDSimulator<double>,double,t._alphaPBTO)
  REGISTER_FLOAT_TYPE("betaPBTO",PBDSimulator<double>,double,t._betaPBTO)
  REGISTER_INT_TYPE("maxIter",PBDSimulator<double>,sizeType,t._maxIter)
  REGISTER_FLOAT_TYPE("tolG",PBDSimulator<double>,double,t._tolG)
  REGISTER_FLOAT_TYPE("tolLS",PBDSimulator<double>,double,t._tolLS)
  REGISTER_FLOAT_TYPE("tolK",PBDSimulator<double>,double,t._tolK)
  REGISTER_FLOAT_TYPE("minDt",PBDSimulator<double>,double,t._minDt)
  REGISTER_FLOAT_TYPE("betaMin",PBDSimulator<double>,double,t._betaMin)
  REGISTER_FLOAT_TYPE("alphaBetaInc",PBDSimulator<double>,double,t._alphaBetaInc)
  REGISTER_BOOL_TYPE("callback",PBDSimulator<double>,bool,t._callback)
  REGISTER_BOOL_TYPE("useKEE",PBDSimulator<double>,bool,t._useKEE)
}
PBDSimulatorTraits<double>::Vec PBDSimulatorTraits<double>::solveQP(qpOASES::SQProblem& prob,const Eigen::Matrix<scalarD,-1,-1,Eigen::RowMajor>& Ad,MatT HId,const Vec& gd,const Vec& wd,T beta)
{
  qpOASES::int_t nWSR=10000;
  HId.diagonal().array()+=std::to_double(beta);
  Cold lbd=-wd,ubAd=Cold::Ones(prob.getNC())-Ad*wd;
  qpOASES::returnValue ret;
  if(!prob.isInitialised())
    ret=prob.init(HId.data(),gd.data(),Ad.data(),lbd.data(),NULL,NULL,ubAd.data(),nWSR);
  else ret=prob.hotstart(HId.data(),gd.data(),Ad.data(),lbd.data(),NULL,NULL,ubAd.data(),nWSR);
  if(!prob.isSolved()) {
    INFOV("qpOASES failed: %d!",ret)
    return Vec::Zero(0);
  } else {
    Vec dwd;
    dwd.resize(wd.size());
    prob.getPrimalSolution(dwd.data());
    return dwd;
  }
}
void PBDSimulatorTraits<__float128>::registerOptions(Options& ops)
{
  REGISTER_FLOAT128_TYPE("alphaPBTO",PBDSimulator<__float128>,__float128,t._alphaPBTO)
  REGISTER_FLOAT128_TYPE("betaPBTO",PBDSimulator<__float128>,__float128,t._betaPBTO)
  REGISTER_INT_TYPE("maxIter",PBDSimulator<__float128>,sizeType,t._maxIter)
  REGISTER_FLOAT128_TYPE("tolG",PBDSimulator<__float128>,__float128,t._tolG)
  REGISTER_FLOAT128_TYPE("tolLS",PBDSimulator<__float128>,__float128,t._tolLS)
  REGISTER_FLOAT128_TYPE("tolK",PBDSimulator<__float128>,__float128,t._tolK)
  REGISTER_FLOAT128_TYPE("minDt",PBDSimulator<__float128>,__float128,t._minDt)
  REGISTER_FLOAT128_TYPE("betaMin",PBDSimulator<__float128>,__float128,t._betaMin)
  REGISTER_FLOAT128_TYPE("alphaBetaInc",PBDSimulator<__float128>,__float128,t._alphaBetaInc)
  REGISTER_BOOL_TYPE("callback",PBDSimulator<__float128>,bool,t._callback)
  REGISTER_BOOL_TYPE("useKEE",PBDSimulator<__float128>,bool,t._useKEE)
}
typename PBDSimulatorTraits<__float128>::Vec PBDSimulatorTraits<__float128>::solveQP(qpOASES::SQProblem& prob,const Eigen::Matrix<scalarD,-1,-1,Eigen::RowMajor>& Ad,MatT HI,const Vec& g,const Vec& w,T beta)
{
  qpOASES::int_t nWSR=10000;
  HI.diagonal().array()+=beta;
  Matd HId=HI.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Cold gd=g.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Cold wd=w.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Cold lbd=-wd,ubAd=Cold::Ones(prob.getNC())-Ad*wd;
  qpOASES::returnValue ret;
  if(!prob.isInitialised())
    ret=prob.init(HId.data(),gd.data(),Ad.data(),lbd.data(),NULL,NULL,ubAd.data(),nWSR);
  else ret=prob.hotstart(HId.data(),gd.data(),Ad.data(),lbd.data(),NULL,NULL,ubAd.data(),nWSR);
  if(!prob.isSolved()) {
    INFOV("qpOASES failed: %d!",ret)
    return Vec::Zero(0);
  } else {
    Cold dwd;
    dwd.resize(wd.size());
    prob.getPrimalSolution(dwd.data());
    return dwd.template cast<T>();
  }
}
void PBDSimulatorTraits<mpfr::mpreal>::registerOptions(Options& ops)
{
  REGISTER_MPFR_TYPE("alphaPBTO",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._alphaPBTO)
  REGISTER_MPFR_TYPE("betaPBTO",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._betaPBTO)
  REGISTER_INT_TYPE("maxIter",PBDSimulator<mpfr::mpreal>,sizeType,t._maxIter)
  REGISTER_MPFR_TYPE("tolG",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._tolG)
  REGISTER_MPFR_TYPE("tolLS",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._tolLS)
  REGISTER_MPFR_TYPE("tolK",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._tolK)
  REGISTER_MPFR_TYPE("minDt",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._minDt)
  REGISTER_MPFR_TYPE("betaMin",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._betaMin)
  REGISTER_MPFR_TYPE("alphaBetaInc",PBDSimulator<mpfr::mpreal>,mpfr::mpreal,t._alphaBetaInc)
  REGISTER_BOOL_TYPE("callback",PBDSimulator<mpfr::mpreal>,bool,t._callback)
  REGISTER_BOOL_TYPE("useKEE",PBDSimulator<mpfr::mpreal>,bool,t._useKEE)
}
typename PBDSimulatorTraits<mpfr::mpreal>::Vec PBDSimulatorTraits<mpfr::mpreal>::solveQP(qpOASES::SQProblem& prob,const Eigen::Matrix<scalarD,-1,-1,Eigen::RowMajor>& Ad,MatT HI,const Vec& g,const Vec& w,T beta)
{
  qpOASES::int_t nWSR=10000;
  HI.diagonal().array()+=beta;
  Matd HId=HI.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Cold gd=g.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Cold wd=w.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Cold lbd=-wd,ubAd=Cold::Ones(prob.getNC())-Ad*wd;
  qpOASES::returnValue ret;
  if(!prob.isInitialised())
    ret=prob.init(HId.data(),gd.data(),Ad.data(),lbd.data(),NULL,NULL,ubAd.data(),nWSR);
  else ret=prob.hotstart(HId.data(),gd.data(),Ad.data(),lbd.data(),NULL,NULL,ubAd.data(),nWSR);
  if(!prob.isSolved()) {
    INFOV("qpOASES failed: %d!",ret)
    return Vec::Zero(0);
  } else {
    Cold dwd;
    dwd.resize(wd.size());
    prob.getPrimalSolution(dwd.data());
    return dwd.template cast<T>();
  }
}
//PBDSimulator
template <typename T>
PBDSimulator<T>::PBDSimulator(const ArticulatedBody& body,Options& ops,const Vec3T& g,PBD_SIMULATOR_MODE integrator):_body(body),_mode(integrator)
{
  if(!ops.hasType<PBDSimulator<T>>())
    PBDSimulatorTraits<T>::registerOptions(ops);
  reset(ops);
  _JRCF.setZero(3,_body.nrJ()*4);
  for(sizeType i=0; i<(sizeType)_body.nrJ(); i++)
    if(_body.joint(i)._M > 0) {
      ROTI(_JRCF,i)=-g*_body.joint(i)._MC.transpose().template cast<T>();
      CTRI(_JRCF,i)=-_body.joint(i)._M*g;
    }
}
template <typename T>
PBDSimulator<T>::PBDSimulator(const PBDSimulator<T>& other):_body(other._body)
{
  operator=(other);
}
template <typename T>
PBDSimulator<T>& PBDSimulator<T>::operator=(const PBDSimulator<T>& other)
{
  _PDTarget=other._PDTarget;
  _wrench=other._wrench;
  _externalWrench=other._externalWrench;
  _mode=other._mode;
  _JRCF=other._JRCF;

  Options ops;
  ops.setOptions(const_cast<PBDSimulator<T>*>(&other));
  reset(ops);
  return *this;
}
template <typename T>
void PBDSimulator<T>::setWrenchConstructor(std::shared_ptr<C0EnvWrenchConstructor<T>> wrench)
{
  _wrench=wrench;
}
template <typename T>
void PBDSimulator<T>::setPDTarget(std::shared_ptr<PDTarget> PDTarget)
{
  _PDTarget=PDTarget;
  if(_PDTarget)
    _PDTarget->reset();
}
template <typename T>
void PBDSimulator<T>::writeVTKSeq(const std::string& path,T horizon,T interval,T dt,std::shared_ptr<PDTarget> PDTarget,const std::map<std::string,std::set<sizeType>>* jointMask)
{
  recreate(path);
  setPDTarget(PDTarget);
  T betaInit=_betaMin,lastDt=dt;
  sizeType id=0,N=_body.nrDOF();
  sizeType intervalFrm=sizeType(interval/dt);
  std::vector<sizeType> iter,iterG;
  std::vector<scalarD> frametime,kinetic;
  std::vector<Vec,Eigen::aligned_allocator<Vec>> qss;
  ASSERT_MSGV(intervalFrm>0,"IntervalFrm(%ld)<=0 is not allowed!",intervalFrm)
  Vec s=concat<Vec>(_qqMqMM[1]._xM,_qqMqMM[2]._xM);
  for(T t=0; t<horizon; t+=dt) {
    TBEG("");
    bool succ=step(dt,lastDt,NULL,s,&betaInit);
    iter.push_back(_iter);
    iterG.push_back(_iterG);
    frametime.push_back(TENDV());
    _qqMqMM[0].reset(_body,s.segment(0,N));
    _qqMqMM[1].reset(_body,s.segment(0,N));
    kinetic.push_back(std::to_double(K(lastDt,mapV((Vec*)NULL),mapM((MatT*)NULL))));
    qss.push_back(s.segment(0,N));
    if(!succ) {
      INFOV("Frm%d--failed!",id)
      break;
    }
    if((id%intervalFrm)==0) {
      INFOV("Frm%d--lastDt=%f--writeVTK!",id,std::to_double(lastDt))
      writeVTK(s.segment(0,N),path+"/frm"+std::to_string(id/intervalFrm)+".vtk",jointMask);
    } else {
      INFOV("Frm%d--lastDt=%f!",id,std::to_double(lastDt))
    }
    id++;
  }
  {
    std::ofstream os(path+"/data.dat");
    for(sizeType i=0; i<(sizeType)iter.size(); i++)
      os << iter[i] << " " << iterG[i] << " " << frametime[i] << " " << kinetic[i] << std::endl;
  }
  {
    std::ofstream os(path+"/pose.dat");
    for(sizeType i=0; i<(sizeType)qss.size(); i++)
      os << qss[i].transpose().unaryExpr([&](const T& in) {
      return (double)std::to_double(in);
    }) << std::endl;
  }
}
template <typename T>
void PBDSimulator<T>::writeVTK(const Vec& q,const std::string& path,const std::map<std::string,std::set<sizeType>>* jointMask) const
{
  ASSERT_MSGV(endsWith(path,".vtk"),"writeVTK path(%s) must ends with .vtk!",path.c_str())
  PBDArticulatedGradientInfo<T> qInfo(_body,q);
  Mat3Xd trans=qInfo._TM.unaryExpr([&](const T& in) {
    return (double)std::to_double(in);
  });
  if(!jointMask)
    _body.writeVTK(trans,path,Joint::MESH);
  else {
    for(const std::pair<std::string,std::set<sizeType>>& item:*jointMask) {
      std::experimental::filesystem::v1::path p(path);
      p=p.parent_path()/(item.first+p.filename().string());
      _body.writeVTK(trans,p.string(),Joint::MESH,&(item.second));
    }
  }
}
template <typename T>
void PBDSimulator<T>::setState(const Vec& qM,const Vec& qMM)
{
  _qqMqMM[1].reset(_body,qM);
  _qqMqMM[2].reset(_body,qMM);
}
template <typename T>
bool PBDSimulator<T>::step(T dt,T& lastDt,const Vec* tau0,Vec& qqMqMM,T* betaInit)
{
  return step(dt,lastDt,mapCV(tau0),mapV(qqMqMM),betaInit);
}
template <typename T>
bool PBDSimulator<T>::step(T dt,T& lastDt,VecCM tau0,VecM qMqMM,T* betaInit)
{
  bool succ=true;
  if(_mode==PBTO) {
    sizeType N=qMqMM.size()/2;
    Vec ret=stepPBTO(dt,lastDt,tau0,mapV2CV(qMqMM));
    qMqMM=concat<Vec>(ret,qMqMM.segment(0,N));
  } else if(_mode==NPBD_SQP) {
    sizeType N=qMqMM.size()/2;
    Vec ret=stepNMDPSQP(dt,lastDt,tau0,mapV2CV(qMqMM));
    qMqMM=concat<Vec>(ret,qMqMM.segment(0,N));
  } else if(_mode==NPBD_PGM || _mode==NPBD_ZOPGM)
    succ=stepNMDPPGMAdaptive(dt,lastDt,tau0,qMqMM,betaInit);
  else {
    ASSERT(false)
    succ=false;
  }
  if(_PDTarget)
    _PDTarget->advance(std::to_double(dt));
  return succ;
}
template <typename T>
bool PBDSimulator<T>::eval(const PBDArticulatedGradientInfo<T>& x,const PBDArticulatedGradientInfo<T>& xM,const PBDArticulatedGradientInfo<T>& xMM,const Vec& ctrl,T dt,T lastDt,T* e,Vec* g,TRIPS* H,sizeType off) const
{
  sizeType iterG=0;
  Vec gV=Vec::Zero(_body.nrDOF());
  MatT h=MatT::Zero(_body.nrDOF(),_body.nrDOF());
  std::vector<ExternalWrench<T>> externalWrench;
  gV=G(dt,lastDt,mapCV(ctrl),mapM(H?&h:NULL),NULL,mapM((MatT*)NULL),x,xM,xMM,externalWrench,iterG);
  if(e)
    *e+=E(dt,dt,mapCV(ctrl),x,xM,xMM);
  if(g)
    g->segment(off,gV.size())+=gV;
  if(H)
    addBlock(*H,off,off,h);
  return true;
}
template <typename T>
void PBDSimulator<T>::debug(T dt,sizeType nrIter,T scale)
{
  Vec3T pos[2],f,f2;
  Mat3T DFDPos,DFDPosN;
  bool useKEE=_useKEE;
  sizeType N=_body.nrDOF();
  T EVal,EVal2,KVal,KVal2,phi0=0.001f;
  Vec GVal,GVal2,w,w2,DKDTheta=Vec::Zero(N),DKDTheta2=Vec::Zero(N),delta,tau0;
  MatT DGDTheta=MatT::Zero(N,N),DDKDDTheta=MatT::Zero(N,N),DGDw;
  std::vector<EndEffectorBounds> EE=_wrench->_externalForces;
  DEFINE_NUMERIC_DELTA_T(T)
  std::vector<std::string> modes;
  modes.push_back("PBD_NO_FORCE");
  modes.push_back("PBTO");
  modes.push_back("NPBD_SQP");
  modes.push_back("NPBD_PGM");
  modes.push_back("NPBD_ZOPGM");
  INFOV("-------------------------------------------------------------DebugPBDSimulator: scale=%s!",std::to_string(scale).c_str())
  for(sizeType pass=0; pass<(sizeType)modes.size(); pass++) {
    _mode=(PBD_SIMULATOR_MODE)pass;
    for(sizeType i=0; i<nrIter; i++) {
      T lastDt=RandEngine::randR(std::to_double(dt*0.5f),std::to_double(dt*1.5f));
      delta=Vec::Random(N);
      tau0=Vec::Random(N);
      _wrench->_externalForces.clear();
      for(sizeType ei=0; ei<_body.nrJ(); ei++) {
        if(RandEngine::randR01()>0.5f) {
          _wrench->_externalForces.push_back(EndEffectorBounds());
          _wrench->_externalForces.back()._localPos.setRandom();
          _wrench->_externalForces.back()._JID.push_back(ei);
          _wrench->_externalForces.back()._phi0=std::to_double(phi0);
        }
      }
      _PDTarget.reset(new PDTarget(PDTarget::Vec::Random(N),PDTarget::Vec::Random(N),PDTarget::Vec::Random(2*N)));
      _qqMqMM[0].reset(_body,Vec::Random(N)*scale);
      _qqMqMM[1].reset(_body,Vec::Random(N)*scale);
      _qqMqMM[2].reset(_body,Vec::Random(N)*scale);
      (*_wrench)(_externalWrench,_qqMqMM[0]);
      for(sizeType KEEVal=0; KEEVal<2; KEEVal++) {
        _useKEE=KEEVal;
        KVal=K(dt,mapV(DKDTheta),mapM(DDKDDTheta));
        EVal=E(dt,lastDt,mapCV(tau0));
        w.setRandom(nW());
        DGDw.resize(N,nW());
        GVal=G(dt,lastDt,mapCV(tau0),mapM(DGDTheta),&w,mapM(DGDw));
        //DGDw
        if(_mode==NPBD_SQP || _mode==NPBD_PGM || _mode==NPBD_ZOPGM) {
          delta.setRandom(nW());
          w2=w+delta*DELTA;
          GVal2=G(dt,lastDt,mapCV(tau0),mapM((MatT*)NULL),&w2,mapM((MatT*)NULL));
          DEBUG_GRADIENT(modes[pass]+"-DGDw",std::sqrt((DGDw*delta).squaredNorm()),std::sqrt((DGDw*delta-(GVal2-GVal)/DELTA).squaredNorm()))
        }
        //G
        delta.setRandom(N);
        _qqMqMM[0].reset(_body,_qqMqMM[0]._xM+delta*DELTA);
        if(_mode==PBD_NO_FORCE) {
          EVal2=E(dt,lastDt,mapCV(tau0));
          DEBUG_GRADIENT(modes[pass]+"-G",GVal.dot(delta),GVal.dot(delta)-(EVal2-EVal)/DELTA)
        }
        //K
        KVal2=K(dt,mapV(DKDTheta2),mapM((MatT*)NULL));
        DEBUG_GRADIENT(modes[pass]+"-DKDTheta-useKEE="+std::to_string(_useKEE?"1":"0"),DKDTheta.dot(delta),DKDTheta.dot(delta)-(KVal2-KVal)/DELTA)
        DEBUG_GRADIENT(modes[pass]+"-DDKDDTheta-useKEE="+std::to_string(_useKEE?"1":"0"),std::sqrt((DDKDDTheta*delta).squaredNorm()),std::sqrt((DDKDDTheta*delta-(DKDTheta2-DKDTheta)/DELTA).squaredNorm()))
        //DGDTheta
        GVal2=G(dt,lastDt,mapCV(tau0),mapM((MatT*)NULL),&w,mapM((MatT*)NULL));
        DEBUG_GRADIENT(modes[pass]+"-DGDTheta",std::sqrt((DGDTheta*delta).squaredNorm()),std::sqrt((DGDTheta*delta-(GVal2-GVal)/DELTA).squaredNorm()))
      }
      //PBTO
      if(_mode==PBTO) {
        pos[0]=Vec3T::Random();
        pos[1]=Vec3T::Random();
        f=forcePBTO(dt,phi0,pos,&DFDPos,&DFDPosN);
        delta=Vec3T::Random();
        pos[0]+=delta*DELTA;
        f2=forcePBTO(dt,phi0,pos,NULL,NULL);
        DEBUG_GRADIENT(modes[pass]+"-DFDPos",std::sqrt((DFDPos*delta).squaredNorm()),std::sqrt((DFDPos*delta-(f2-f)/DELTA).squaredNorm()))

        f=forcePBTO(dt,phi0,pos,&DFDPos,&DFDPosN);
        delta=Vec3T::Random();
        pos[1]+=delta*DELTA;
        f2=forcePBTO(dt,phi0,pos,NULL,NULL);
        DEBUG_GRADIENT(modes[pass]+"-DFDPosN",std::sqrt((DFDPosN*delta).squaredNorm()),std::sqrt((DFDPosN*delta-(f2-f)/DELTA).squaredNorm()))

        T e=0,e2=0;
        TRIPS h;
        Vec ctrl=Vec::Random(N);
        sizeType off=RandEngine::randI(1,10);
        Vec g=Vec::Zero(N+off),g2=g,dx=Vec::Random(N+off);
        PBDArticulatedGradientInfo<T> q(_body,Vec::Random(N)*scale);
        PBDArticulatedGradientInfo<T> q2(_body,q._xM+dx.segment(off,N)*DELTA);
        PBDArticulatedGradientInfo<T> qM(_body,Vec::Random(N)*scale);
        PBDArticulatedGradientInfo<T> qMM(_body,Vec::Random(N)*scale);
        eval(q,qM,qMM,ctrl,dt,dt/2,&e,&g,&h,off);
        eval(q2,qM,qMM,ctrl,dt,dt/2,&e2,&g2,NULL,off);
        Eigen::SparseMatrix<T,0,sizeType> H;
        H.resize(off+N,off+N);
        H.setFromTriplets(h.begin(),h.end());
        DEBUG_GRADIENT(modes[pass]+"-evalE",g.dot(dx),g.dot(dx)-(e2-e)/DELTA)
        DEBUG_GRADIENT(modes[pass]+"-evalG",std::sqrt((H*dx).squaredNorm()),std::sqrt((H*dx-(g2-g)/DELTA).squaredNorm()))
      }
    }
  }
  _wrench->_externalForces=EE;
  _useKEE=useKEE;
}
template <typename T>
void PBDSimulator<T>::reset(Options& ops)
{
  _alphaPBTO=1000;
  _betaPBTO=1000;
  _tolG=1E-8f;
  _tolLS=1E-6f;
  _tolK=1E-6f;
  _minDt=1E-5f;
  _betaMin=1e-4f;
  _alphaBetaInc=1.5f;
  _maxIter=1e6;
  _callback=true;
  _useKEE=true;
  ops.setOptions(this);
}
template <typename T>
void PBDSimulator<T>::setMode(PBD_SIMULATOR_MODE m)
{
  _mode=m;
}
template <typename T>
PBD_SIMULATOR_MODE PBDSimulator<T>::getMode() const
{
  return _mode;
}
//helper
template <typename T>
std::shared_ptr<qpOASES::DenseMatrix> createAQP(sizeType nW,const std::vector<ExternalWrench<T>>& f)
{
  scalarD* val=new scalarD[f.size()*nW];
  for(sizeType i=0; i<(sizeType)f.size()*nW; i++)
    val[i]=0;
  for(sizeType i=0,offW=0; i<(sizeType)f.size(); i++) {
    for(sizeType j=0; j<f[i]._B.cols(); j++)
      val[i*nW+offW+j]=1;
    offW+=f[i]._B.cols();
  }
  std::shared_ptr<qpOASES::DenseMatrix> Ad(new qpOASES::DenseMatrix(f.size(),nW,nW,val));
  Ad->doFreeMemory();
  return Ad;
}
template <typename T>
typename PBDSimulator<T>::Vec PBDSimulator<T>::stepPBTO(T dt,T lastDt,VecCM tau0,VecCM qMqMM)
{
  _iterG=0;
  sizeType N=_body.nrDOF();
  Vec init=qMqMM.segment(0,N);
  _qqMqMM[1].reset(_body,qMqMM.segment(0,N));
  _qqMqMM[2].reset(_body,qMqMM.segment(N,N));
  return manifoldProjection(dt,lastDt,tau0,mapCV(init),NULL);
}
template <typename T>
typename PBDSimulator<T>::Vec PBDSimulator<T>::stepNMDPSQP(T dt,T lastDt,VecCM tau0,VecCM qMqMM)
{
#ifdef OPTIMIZER_SUPPORT
  _iterG=0;
  sizeType N=_body.nrDOF();
  _qqMqMM[1].reset(_body,qMqMM.segment(0,N));
  _qqMqMM[2].reset(_body,qMqMM.segment(N,N));
  (*_wrench)(_externalWrench,_qqMqMM[1]);
  PBDNMDPObjective<T> obj(dt,lastDt,tau0,*this);

  Vec x=Vec::Zero(obj.inputs());
#if defined(KNITRO_SUPPORT)
  KNInterface<T> KNSol(obj);
  KNSol.solve(false,0,std::to_double(_tolG),std::to_double(_tolG),true,_maxIter,1.0f,x);
  return x.segment(0,_body.nrDOF());
#elif defined(IPOPT_SUPPORT)
  IPOPTInterface<T>::optimize(x,obj,std::to_double(_tolG),_maxIter,0,1.0f,true);
  return x.segment(0,_body.nrDOF());
#else
  return Vec::Zero(0);
#endif
#else
  FUNCTION_NOT_IMPLEMENTED
  return Vec::Zero(0);
#endif
}
template <typename T>
typename PBDSimulator<T>::Vec PBDSimulator<T>::stepNMDPPGM(T dt,T lastDt,VecCM tau0,VecCM qMqMM,T* betaInit)
{
  sizeType N=_body.nrDOF();
  T betaMin=_betaMin,beta=std::max<T>(betaInit?*betaInit:betaMin,betaMin),inc=_alphaBetaInc,alpha=1;
  _qqMqMM[0].reset(_body,qMqMM.segment(0,N));
  _qqMqMM[1]=_qqMqMM[0];
  _qqMqMM[2].reset(_body,qMqMM.segment(N,N));
  MatT DDKDDTheta=MatT::Zero(N,N),DGDTheta=MatT::Zero(N,N),DGDw,H;
  Vec DKDTheta=Vec::Zero(N),q=qMqMM.segment(0,N),qNext,w=Vec::Zero(0),wNext=Vec::Zero(0),g;
  bool updateQP=true;
  //QP constraint matrix
  Eigen::Matrix<scalarD,-1,-1,Eigen::RowMajor> Ad;
  //initial projection
  q=manifoldProjection(dt,lastDt,tau0,mapCV(q),&w,&alpha);
  if(q.size()==0)
    return q;
  DGDw.resize(N,nW());
  T KVal=K(dt,mapV((Vec*)NULL),mapM((MatT*)NULL)),KValNext;
  qpOASES::SQProblem prob(w.size(),_externalWrench.size());
  for(; _iter<_maxIter; _iter++) {
    //buildQP
    if(updateQP) {
      K(dt,mapV(DKDTheta),mapM(DDKDDTheta));
      G(dt,lastDt,tau0,mapM(DGDTheta),&w,mapM(DGDw));
      SolveNewton<T>::solveLU(DGDTheta,DGDw,H,true);
      g=H.transpose()*DKDTheta;
      H=H.transpose()*(DDKDDTheta*H);
      if(Ad.size()==0) {
        Ad.setZero(_externalWrench.size(),nW());
        for(sizeType i=0,offW=0; i<(sizeType)_externalWrench.size(); offW+=_externalWrench[i]._B.cols(),i++)
          Ad.block(i,offW,1,_externalWrench[i]._B.cols()).setOnes();
      }
    }
    //solveQP
    Vec dwd=PBDSimulatorTraits<T>::solveQP(prob,Ad,H,g,w,beta);
    if(std::sqrt(dwd.squaredNorm())<_tolG) {
      if(_callback) {
        INFOV("Beta=%f!",std::to_double(beta))
      }
      break;
    }
    if(dwd.size()==0) {
      beta*=inc;
      updateQP=false;
      if(beta>1/_tolLS) {
        if(_callback) {
          INFOV("Line-Search failed (iter=%d)!",_iter+1)
        }
        q.resize(0);
        break;
      } else continue;
    } else wNext=w+dwd;
    //manifold projection
    qNext=manifoldProjection(dt,lastDt,tau0,mapCV(q),&wNext,&alpha);
    if(qNext.size()==0)
      return qNext;
    else {
      _qqMqMM[0].reset(_body,qNext);
      KValNext=K(dt,mapV((Vec*)NULL),mapM((MatT*)NULL));
    }
    //update line search info
    if(KValNext<KVal) {
      w=wNext;
      q=qNext;
      if(_callback) {
        INFOV("NPBD-Iter%d: K=%f, beta=%f!",_iter,std::to_double(KVal),std::to_double(beta))
      }
      if(KVal-KValNext<_tolK)
        break;
      KVal=KValNext;
      updateQP=true;
      beta=std::max<T>(beta/inc,betaMin);
    } else {
      updateQP=false;
      beta*=inc;
      if(beta>1/_tolLS) {
        if(_callback) {
          INFOV("Line-Search failed (iter=%d)!",_iter+1)
        }
        q.resize(0);
        break;
      }
    }
  }
  //debug
  /*for(sizeType i=0,off=0; i<(sizeType)_externalWrench.size(); i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    Vec3T f=(wi._B*w.segment(off,wi._B.cols())). template segment<3>(3);
    INFOV("f=(%f,%f,%f), w=%f",std::to_double(f[0]),std::to_double(f[1]),std::to_double(f[2]),std::to_double(w.segment(off,wi._B.cols()).sum()))
    off+=wi._B.cols();
  }*/
  if(q.size()>0 && betaInit)
    *betaInit=beta;
  return q;
}
template <typename T>
bool PBDSimulator<T>::stepNMDPPGMAdaptive(T dt,T& lastDt,VecCM tau0,VecM qMqMM,T* betaInit)
{
  _iter=0,_iterG=0;
  if(dt<_minDt)  //this algorithm has failed
    return false;
  sizeType N=qMqMM.size()/2;
  T betaInitTmp=betaInit?*betaInit:_betaMin;
  Vec ret=stepNMDPPGM(dt,lastDt,tau0,mapV2CV(qMqMM),&betaInitTmp);
  if(ret.size()==0) {
    //subdivide: piece-1
    betaInitTmp=betaInit?*betaInit:_betaMin;
    if(!stepNMDPPGMAdaptive(dt/2,lastDt,tau0,qMqMM,&betaInitTmp))
      return false;
    //subdivide: piece-2
    if(!stepNMDPPGMAdaptive(dt/2,lastDt,tau0,qMqMM,&betaInitTmp))
      return false;
    if(betaInit)
      *betaInit=betaInitTmp;
  } else {
    if(betaInit)
      *betaInit=betaInitTmp;
    qMqMM.segment(N,N)=qMqMM.segment(0,N);
    qMqMM.segment(0,N)=ret;
    lastDt=dt;
  }
  return true;
}
template <typename T>
typename PBDSimulator<T>::Vec PBDSimulator<T>::manifoldProjection(T dt,T lastDt,VecCM tau0,VecCM init,Vec* w,T* alphaInit)
{
  sizeType N=init.size();
  Vec q=init,qNext,Gi,dir;
  MatT DGDTheta=MatT::Zero(N,N);
  T alpha=alphaInit?*alphaInit:1,inc=_alphaBetaInc,V;
  T tolLS=1;//_tolLS;
  _qqMqMM[0].reset(_body,q);
  Gi=G(dt,lastDt,tau0,mapM(DGDTheta),w,mapM((MatT*)NULL));
  for(sizeType i=0; i<_maxIter; i++) {
    //termination
    if(std::sqrt(Gi.squaredNorm())<_tolG) {
      if(_callback) {
        INFOV("Manifold-Projection succeed (it=%ld)f!",i)
      }
      break;
    }
    V=Gi.squaredNorm();
    SolveNewton<T>::template solveLU<Vec>(DGDTheta,Gi,dir,true);
    //line search
    while(alpha>=tolLS) {
      qNext=q+dir*alpha;
      _qqMqMM[0].reset(_body,qNext);
      Gi=G(dt,lastDt,tau0,mapM(DGDTheta),w,mapM((MatT*)NULL));
      if(Gi.squaredNorm()<V) {
        q=qNext;
        break;
      } else alpha/=inc;
    }
    if(alpha<tolLS) {
      if(_callback) {
        INFOV("Manifold-Projection failed (iter=%d,dt=%f,lastDt=%f,GNorm=%f)!",
              i+1,std::to_double(dt),std::to_double(lastDt),std::to_double(std::sqrt(Gi.squaredNorm())))
      }
      alpha=alphaInit?*alphaInit:1;
      q.resize(0);
      break;
    }
    alpha=std::min<T>(alpha*inc,1.0f);
  }
  if(alphaInit)
    *alphaInit=alpha;
  return q;
}
//helper
template <typename T>
T PBDSimulator<T>::E(T dt,T lastDt,VecCM tau0) const
{
  return E(dt,lastDt,tau0,_qqMqMM[0],_qqMqMM[1],_qqMqMM[2]);
}
template <typename T>
T PBDSimulator<T>::E
(T dt,T lastDt,VecCM tau0,
 const PBDArticulatedGradientInfo<T>& q,
 const PBDArticulatedGradientInfo<T>& qM,
 const PBDArticulatedGradientInfo<T>& qMM) const
{
  T ret=0.0;
  Mat3X4T A;
  T alpha=dt/lastDt;
  T coef=1.0/(lastDt*lastDt*alpha*alpha);
  sizeType nrJ=_body.nrJ();
  for(sizeType k=0; k<nrJ; k++) {
    const Joint& J=_body.joint(k);
    const Mat3T PPT=J._MCCT.template cast<T>();
    const Vec3T P=J._MC.template cast<T>();
    A=TRANSI(q._TM,k)-(1+alpha)*TRANSI(qM._TM,k)+alpha*TRANSI(qMM._TM,k);
    ret+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace()*coef/2;
    ret+=(TRANSI(q._TM,k)*TRANSI(_JRCF,k).transpose()).trace();
  }
  if(tau0.data())
    ret-=q._xM.dot(tau0);
  if(_PDTarget) {
    const Vec s=_PDTarget->s().template cast<T>();
    sizeType N=s.size()/2;
    Vec tmp=q._xM-s.segment(0,N);
    ret+=(tmp.cwiseProduct(_PDTarget->_PCoef.template cast<T>())).dot(tmp)/2;
    tmp=(q._xM-qM._xM)/dt-s.segment(N,N);
    ret+=(tmp.cwiseProduct(_PDTarget->_DCoef.template cast<T>())).dot(tmp)*dt/2;
  }
  return ret;
}
template <typename T>
T PBDSimulator<T>::K(T dt,VecM DKDTheta,MatTM DDKDDTheta) const
{
  T ret=0.0;
  Mat3X4T A;
  Mat3XT G,GB,MRR,MRt,MtR,Mtt;
  T coef=1.0/(dt*dt);
  sizeType nrJ=_body.nrJ();
  sizeType nrE=(sizeType)_externalWrench.size();
  if(DKDTheta.data()) {
    DKDTheta.setZero();
    G.setZero(3,4*nrJ);
    GB.setZero(3,4*nrJ);
  }
  if(DDKDDTheta.data()) {
    DDKDDTheta.setZero();
    MRR.setZero(3,nrJ*3);
    MRt.setZero(3,nrJ*3);
    MtR.setZero(3,nrJ*3);
    Mtt.setZero(3,nrJ*3);
  }
  if(_useKEE) {
    for(sizeType i=0; i<nrE; i++) {
      const EndEffectorBounds& ee=_wrench->_externalForces[i];
      const Mat3T PPT=(ee._localPos*ee._localPos.transpose()).template cast<T>();
      const Vec3T P=ee._localPos.template cast<T>();
      sizeType k=ee.jointId();
      //if(_wrench->_env->phi(ROTI(_qqMqMM[0]._TM,k)*P+CTRI(_qqMqMM[0]._TM,k))>ee._phi0)
      //  continue;
      A=TRANSI(_qqMqMM[0]._TM,k)-TRANSI(_qqMqMM[1]._TM,k);
      ret+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()).trace();
      if(DKDTheta.data() || DDKDDTheta.data()) {
        ROTI(G,k)+=ROT(A)*PPT+CTR(A)*P.transpose();
        CTRI(G,k)+=CTR(A)+ROT(A)*P;
      }
    }
    if(DKDTheta.data()) {
      _qqMqMM[0].DTG(_body,mapM(GB=G),DKDTheta);
      DKDTheta*=coef;
    }
    if(DDKDDTheta.data()) {
      for(sizeType i=0; i<nrE; i++) {
        const EndEffectorBounds& ee=_wrench->_externalForces[i];
        const Mat3T PPT=(ee._localPos*ee._localPos.transpose()).template cast<T>();
        const Vec3T P=ee._localPos.template cast<T>();
        sizeType k=ee.jointId();
        //if(_wrench->_env->phi(ROTI(_qqMqMM[0]._TM,k)*P+CTRI(_qqMqMM[0]._TM,k))>ee._phi0)
        //  continue;
        MRR.template block<3,3>(0,k*3)-=invDoubleCrossMatTrace<T>(ROTI(_qqMqMM[0]._TM,k)*PPT*ROTI(_qqMqMM[0]._TM,k).transpose());
        MRt.template block<3,3>(0,k*3)+=cross<T>(ROTI(_qqMqMM[0]._TM,k)*P);
        MtR.template block<3,3>(0,k*3)-=cross<T>(ROTI(_qqMqMM[0]._TM,k)*P);
        Mtt.template block<3,3>(0,k*3)+=Mat3T::Identity();
      }
      _qqMqMM[0].toolAB(_body,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),mapM(GB=G),[&](sizeType row,sizeType col,T val) {
        DDKDDTheta(row,col)+=val;
      });
      DDKDDTheta*=coef;
    }
  } else {
    for(sizeType k=0; k<nrJ; k++) {
      const Joint& J=_body.joint(k);
      const Mat3T PPT=J._MCCT.template cast<T>();
      const Vec3T P=J._MC.template cast<T>();
      A=TRANSI(_qqMqMM[0]._TM,k)-TRANSI(_qqMqMM[1]._TM,k);
      ret+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace();
      if(DKDTheta.data() || DDKDDTheta.data()) {
        ROTI(G,k)+=ROT(A)*PPT+CTR(A)*P.transpose();
        CTRI(G,k)+=CTR(A)*J._M+ROT(A)*P;
      }
    }
    if(DKDTheta.data()) {
      _qqMqMM[0].DTG(_body,mapM(GB=G),DKDTheta);
      DKDTheta*=coef;
    }
    if(DDKDDTheta.data()) {
      for(sizeType k=0; k<nrJ; k++) {
        const Joint& J=_body.joint(k);
        MRR.template block<3,3>(0,k*3)-=invDoubleCrossMatTrace<T>(ROTI(_qqMqMM[0]._TM,k)*J._MCCT.template cast<T>()*ROTI(_qqMqMM[0]._TM,k).transpose());
        MRt.template block<3,3>(0,k*3)+=cross<T>(ROTI(_qqMqMM[0]._TM,k)*J._MC.template cast<T>());
        MtR.template block<3,3>(0,k*3)-=cross<T>(ROTI(_qqMqMM[0]._TM,k)*J._MC.template cast<T>());
        Mtt.template block<3,3>(0,k*3)+=Mat3T::Identity()*J._M;
      }
      _qqMqMM[0].toolAB(_body,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),mapM(GB=G),[&](sizeType row,sizeType col,T val) {
        DDKDDTheta(row,col)+=val;
      });
      DDKDDTheta*=coef;
    }
  }
  return ret*coef/2;
}
template <typename T>
typename PBDSimulator<T>::Vec PBDSimulator<T>::G(T dt,T lastDt,VecCM tau0,MatTM DGDTheta,Vec* w,MatTM DGDw)
{
  return G(dt,lastDt,tau0,DGDTheta,w,DGDw,_qqMqMM[0],_qqMqMM[1],_qqMqMM[2],_externalWrench,_iterG);
}
template <typename T>
typename PBDSimulator<T>::Vec PBDSimulator<T>::G
(T dt,T lastDt,VecCM tau0,MatTM DGDTheta,Vec* w,MatTM DGDw,
 const PBDArticulatedGradientInfo<T>& q,
 const PBDArticulatedGradientInfo<T>& qM,
 const PBDArticulatedGradientInfo<T>& qMM,
 std::vector<ExternalWrench<T>>& externalWrench,sizeType& iterG) const
{
  iterG++;
  Mat3X4T A;
  Mat3T DFDPos;
  Mat3XT G,GB,MRR,MRt,MtR,Mtt;
  T alpha=dt/lastDt;
  T coef=1.0/(lastDt*lastDt*alpha*alpha);
  sizeType nrJ=_body.nrJ();
  G.setZero(3,4*nrJ);
  GB.setZero(3,4*nrJ);
  if(DGDTheta.data()) {
    MRR.setZero(3,nrJ*3);
    MRt.setZero(3,nrJ*3);
    MtR.setZero(3,nrJ*3);
    Mtt.setZero(3,nrJ*3);
    DGDTheta.setZero();
  }
  //G
  for(sizeType k=0; k<nrJ; k++) {
    const Joint& J=_body.joint(k);
    const Mat3T PPT=J._MCCT.template cast<T>();
    const Vec3T P=J._MC.template cast<T>();
    A=(TRANSI(q._TM,k)-(1+alpha)*TRANSI(qM._TM,k)+alpha*TRANSI(qMM._TM,k))*coef;
    ROTI(G,k)+=ROT(A)*PPT+CTR(A)*P.transpose();
    CTRI(G,k)+=CTR(A)*J._M+ROT(A)*P;
    TRANSI(G,k)+=TRANSI(_JRCF,k);
  }
  if(_mode==PBTO && _wrench && _wrench->_env) {
    for(sizeType i=0; i<(sizeType)_wrench->_externalForces.size(); i++) {
      const EndEffectorBounds& ee=_wrench->_externalForces[i];
      Vec3T fPBTO=forcePBTO(dt,ee,DGDTheta.data()?&DFDPos:NULL,NULL);
      if(std::abs(fPBTO[0])>0 || std::abs(fPBTO[1])>0 || std::abs(fPBTO[2])>0) {
        if(DGDTheta.data()) {
          Mat3T CRp=-cross<T>(ROTI(q._TM,ee.jointId())*ee._localPos.template cast<T>());
          MRR.template block<3,3>(0,ee.jointId()*3)-=CRp*DFDPos*CRp.transpose();
          MtR.template block<3,3>(0,ee.jointId()*3)+=DFDPos*CRp.transpose();
          MRt.template block<3,3>(0,ee.jointId()*3)+=CRp*DFDPos;
          Mtt.template block<3,3>(0,ee.jointId()*3)-=DFDPos;
        }
        ROTI(G,ee.jointId())-=fPBTO*ee._localPos.transpose().template cast<T>();
        CTRI(G,ee.jointId())-=fPBTO;
      }
    }
  }
  if((_mode==NPBD_SQP || _mode==NPBD_PGM || _mode==NPBD_ZOPGM) && _wrench) {
    (*_wrench)(externalWrench,q,_mode!=NPBD_ZOPGM);
    ASSERT_MSGV(w && externalWrench.size()==_wrench->_externalForces.size(),"size of externalWrench (%d) mismatch size of externalForce (%d)",externalWrench.size(),_wrench->_externalForces.size())
    if(w->size()==0)
      w->setZero(nW());
    for(sizeType i=0,offW=0; i<(sizeType)externalWrench.size(); i++) {
      ExternalWrench<T>& wi=externalWrench[i];
      //wi*=dt;
      sizeType nW=wi._B.cols();
      wi._w=w->segment(offW,nW);
      const EndEffectorBounds& ee=_wrench->_externalForces[i];
      Vec3T f=wi._B.block(3,0,3,nW)*w->segment(offW,nW);
      if(DGDw.data() || DGDTheta.data()) {
        SMat dPos=DPos(ee);
        if(DGDw.data())
          DGDw.block(0,offW,q._xM.size(),nW)=-dPos*wi._B.block(3,0,3,nW);
        if(DGDTheta.data() && _mode!=NPBD_ZOPGM)
          DFDTheta(DGDTheta,wi,ee,dPos,w->segment(offW,nW));
      }
      ROTI(G,ee.jointId())-=f*ee._localPos.transpose().template cast<T>();
      CTRI(G,ee.jointId())-=f;
      offW+=nW;
    }
  }
  //FD
  Vec FD=Vec::Zero(_body.nrDOF());
  q.DTG(_body,mapM(GB=G),mapV(FD));
  if(tau0.data())
    FD-=tau0;
  if(_PDTarget) {
    const Vec s=_PDTarget->s().template cast<T>();
    sizeType N=s.size()/2;
    FD+=(q._xM-s.segment(0,N)).cwiseProduct(_PDTarget->_PCoef.template cast<T>());
    FD+=((q._xM-qM._xM)/dt-s.segment(N,N)).cwiseProduct(_PDTarget->_DCoef.template cast<T>());
    if(DGDTheta.data()) {
      DGDTheta.diagonal()+=_PDTarget->_PCoef.template cast<T>();
      DGDTheta.diagonal()+=_PDTarget->_DCoef.template cast<T>()/dt;
    }
  }
  //DGDTheta
  if(DGDTheta.data()) {
    for(sizeType k=0; k<nrJ; k++) {
      const Joint& J=_body.joint(k);
      MRR.template block<3,3>(0,k*3)-=invDoubleCrossMatTrace<T>(ROTI(q._TM,k)*J._MCCT.template cast<T>()*ROTI(q._TM,k).transpose())*coef;
      MRt.template block<3,3>(0,k*3)+=cross<T>(ROTI(q._TM,k)*J._MC.template cast<T>())*coef;
      MtR.template block<3,3>(0,k*3)-=cross<T>(ROTI(q._TM,k)*J._MC.template cast<T>())*coef;
      Mtt.template block<3,3>(0,k*3)+=Mat3T::Identity()*J._M*coef;
    }
    q.toolAB(_body,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),mapM(GB=G),[&](sizeType row,sizeType col,T val) {
      DGDTheta(row,col)+=val;
    });
  }
  return FD;
}
template <typename T>
typename PBDSimulator<T>::Vec3T PBDSimulator<T>::forcePBTO(T dt,T phi0,const Vec3T qqM[2],Mat3T* DFDPos,Mat3T* DFDPosN) const
{
  const Vec3T& pos=qqM[0];
  const Vec3T& posN=qqM[1];
  Vec3T vel=(pos-posN)/dt,velT,dPhi,n,force=Vec3T::Zero();
  Mat3T h;
  if(DFDPos)
    DFDPos->setZero();
  if(DFDPosN)
    DFDPosN->setZero();
  T phi=_wrench->_env->phi(pos,DFDPos?&dPhi:NULL)+phi0,velN;
  if(phi<0) {
    n=_wrench->_env->phiGrad(pos,DFDPos?&h:NULL);
    //normal
    force-=_alphaPBTO*phi*n;
    if(DFDPos)
      *DFDPos-=_alphaPBTO*(n*dPhi.transpose()+phi*h);
    //tangent
    velN=vel.dot(n);
    velT=vel-velN*n;
    force+=_betaPBTO*phi*velT;
    if(DFDPos)
      *DFDPos+=_betaPBTO*(velT*dPhi.transpose()+phi*(Mat3T::Identity()-n*n.transpose())/dt-phi*(velN*h+n*vel.transpose()*h));
    if(DFDPosN)
      *DFDPosN-=_betaPBTO*phi*(Mat3T::Identity()-n*n.transpose())/dt;
  }
  return force;
}
template <typename T>
typename PBDSimulator<T>::Vec3T PBDSimulator<T>::forcePBTO(T dt,const EndEffectorBounds& EE,Mat3T* DFDPos,Mat3T* DFDPosN) const
{
  Vec3T pos[2];
  pos[0]=ROTI(_qqMqMM[0]._TM,EE.jointId())*EE._localPos.template cast<T>()+CTRI(_qqMqMM[0]._TM,EE.jointId());
  pos[1]=ROTI(_qqMqMM[1]._TM,EE.jointId())*EE._localPos.template cast<T>()+CTRI(_qqMqMM[1]._TM,EE.jointId());
  return forcePBTO(dt,EE._phi0,pos,DFDPos,DFDPosN);
}
template <typename T>
void PBDSimulator<T>::DFDTheta(MatTM DGDTheta,const ExternalWrench<T>& ew,const EndEffectorBounds& ee,const SMat& dPos,const Vec& w) const
{
  SMat DFDq;
  STrips trips;
  if(!ew._DBDq.empty()) {
    for(sizeType i=0; i<(sizeType)ew._DBDq.size(); i++)
      addBlock(trips,0,ew._DBDq[i].first,(ew._DBDq[i].second*w).template segment<3>(3));
  } else {
    Mat3X4T A;
    for(sizeType d=0; d<3; d++) {
      for(sizeType r=0; r<3; r++)
        for(sizeType c=0; c<4; c++)
          A(r,c)=(ew._DBDX[r][c]*w)[3+d];
      _qqMqMM[0].DTG(ee.jointId(),_body,A,[&](sizeType jid,T coef) {
        trips.push_back(STrip(d,jid,coef));
      });
    }
  }
  DFDq.resize(3,DGDTheta.rows());
  DFDq.setFromTriplets(trips.begin(),trips.end());
  DGDTheta-=dPos*DFDq;
}
template <typename T>
typename PBDSimulator<T>::SMat PBDSimulator<T>::DPos(const EndEffectorBounds& EE) const
{
  Mat3X4T A;
  SMat DPos;
  STrips trips;
  for(sizeType d=0; d<3; d++) {
    A.setZero();
    A.template block<1,3>(d,0)=EE._localPos.transpose().template cast<T>();
    A(d,3)=1;
    _qqMqMM[0].DTG(EE.jointId(),_body,A,[&](sizeType jid,T coef) {
      trips.push_back(STrip(jid,d,coef));
    });
  }
  DPos.resize(_qqMqMM[0]._xM.size(),3);
  DPos.setFromTriplets(trips.begin(),trips.end());
  return DPos;
}
template <typename T>
sizeType PBDSimulator<T>::nW() const
{
  sizeType nW=0;
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++)
    nW+=_externalWrench[i]._B.cols();
  return nW;
}
//objective
#ifdef OPTIMIZER_SUPPORT
template <typename T>
PBDNMDPObjective<T>::PBDNMDPObjective(T dt,T lastDt,VecCM tau0,PBDSimulator<T>& sim):_dt(dt),_lastDt(lastDt),_tau0(tau0),_sim(sim) {}
template <typename T>
typename PBDNMDPObjective<T>::Vec PBDNMDPObjective<T>::lb() const
{
  return concat<Vec>(Vec::Constant(_sim._body.nrDOF(),-DSSQPObjective<T>::infty()),Vec::Zero(_sim.nW()));
}
template <typename T>
typename PBDNMDPObjective<T>::Vec PBDNMDPObjective<T>::ub() const
{
  return concat<Vec>(Vec::Constant(_sim._body.nrDOF(), DSSQPObjective<T>::infty()),Vec::Ones(_sim.nW()));
}
template <typename T>
typename PBDNMDPObjective<T>::Vec PBDNMDPObjective<T>::gl() const
{
  return concat<Vec>(Vec::Zero(_sim._body.nrDOF()),Vec::Constant((sizeType)_sim._externalWrench.size(),0));
}
template <typename T>
typename PBDNMDPObjective<T>::Vec PBDNMDPObjective<T>::gu() const
{
  return concat<Vec>(Vec::Zero(_sim._body.nrDOF()),Vec::Constant((sizeType)_sim._externalWrench.size(),1));
}
template <typename T>
typename PBDNMDPObjective<T>::Vec PBDNMDPObjective<T>::init() const
{
  return concat<Vec>(_sim._qqMqMM[1]._xM,Vec::Zero(_sim.nW()));
}
//constraint
template <typename T>
int PBDNMDPObjective<T>::operator()(const Vec& x,Vec& fvec,STrips* fjac)
{
  sizeType N=_sim._qqMqMM[1]._xM.size();
  _sim._qqMqMM[0].reset(_sim._body,x.segment(0,N));
  (*_sim._wrench)(_sim._externalWrench,_sim._qqMqMM[0]);
  Vec w=Eigen::Map<const Vec>(x.data()+N,_sim.nW());
  if(fjac) {
    MatT DGDTheta,DGDw;
    DGDTheta.resize(N,N);
    DGDw.resize(N,_sim.nW());
    fvec.setZero(N+_sim._externalWrench.size());
    fvec.segment(0,N)=_sim.G(_dt,_lastDt,_tau0,mapM(DGDTheta),&w,mapM(DGDw));
    addBlock(*fjac,0,0,DGDTheta);
    addBlock(*fjac,0,N,DGDw);
    for(sizeType i=0,offW=N; i<(sizeType)_sim._externalWrench.size(); i++) {
      fvec[N+i]=x.segment(offW,_sim._externalWrench[i]._B.cols()).sum();
      addBlock(*fjac,N+i,offW,MatT::Constant(1,_sim._externalWrench[i]._B.cols(),1));
      offW+=_sim._externalWrench[i]._B.cols();
    }
  } else {
    fvec.setZero(N+_sim._externalWrench.size());
    fvec.segment(0,N)=_sim.G(_dt,_lastDt,_tau0,mapM((MatT*)NULL),&w,mapM((MatT*)NULL));
    for(sizeType i=0,offW=N; i<(sizeType)_sim._externalWrench.size(); i++) {
      fvec[N+i]=x.segment(offW,_sim._externalWrench[i]._B.cols()).sum();
      offW+=_sim._externalWrench[i]._B.cols();
    }
  }
  return 0;
}
//objective
template <typename T>
T PBDNMDPObjective<T>::operator()(const Vec& x,Vec* fgrad)
{
  sizeType N=_sim._qqMqMM[1]._xM.size();
  _sim._qqMqMM[0].reset(_sim._body,x.segment(0,N));
  if(fgrad) {
    Vec DKDTheta=Vec::Zero(N);
    T ret=_sim.K(_dt,mapV(DKDTheta),mapM((MatT*)NULL));
    fgrad->setZero(N+_sim.nW());
    fgrad->segment(0,N)=DKDTheta;
    return ret;
  } else return _sim.K(_dt,mapV((Vec*)NULL),mapM((MatT*)NULL));
}
//problem size
template <typename T>
int PBDNMDPObjective<T>::inputs() const
{
  return _sim._body.nrDOF()+_sim.nW();
}
template <typename T>
int PBDNMDPObjective<T>::values() const
{
  return _sim._body.nrDOF()+(sizeType)_sim._externalWrench.size();
}
#endif
//instance
template class PBDSimulator<double>;
#ifdef ALL_TYPES
template class PBDSimulator<__float128>;
template class PBDSimulator<mpfr::mpreal>;
#endif
PRJ_END
#endif
