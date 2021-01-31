#ifdef ENVIRONMENT_SUPPORT
#include "MDPSimulator.h"
#include "PBDSimulator.h"
#include "PDTarget.h"
#include <CommonFile/Timing.h>
#include <Optimizer/QCQPSolverMosek.h>
#include <Environment/GranularWrenchConstructor.h>
#include <Utils/SpatialRotationUtil.h>
#include <Utils/DebugGradient.h>
#include <Utils/Utils.h>
#include <omp.h>

PRJ_BEGIN

template <>
struct MDPSimulatorTraits<double>
{
  typedef double T;
  DECL_MAP_TYPES_T
  static void registerOptions(Options& ops)
  {
    REGISTER_FLOAT_TYPE("tolK",MDPSimulator<double>,double,t._tolK)
    REGISTER_FLOAT_TYPE("minDt",MDPSimulator<double>,double,t._minDt)
    REGISTER_FLOAT_TYPE("betaMin",MDPSimulator<double>,double,t._betaMin)
    REGISTER_FLOAT_TYPE("alphaBetaInc",MDPSimulator<double>,double,t._alphaBetaInc)
    REGISTER_BOOL_TYPE("useKEE",MDPSimulator<double>,bool,t._useKEE)
  }
};
template <>
struct MDPSimulatorTraits<__float128>
{
  typedef __float128 T;
  DECL_MAP_TYPES_T
  static void registerOptions(Options& ops)
  {
    REGISTER_FLOAT128_TYPE("tolK",MDPSimulator<__float128>,__float128,t._tolK)
    REGISTER_FLOAT128_TYPE("minDt",MDPSimulator<__float128>,__float128,t._minDt)
    REGISTER_FLOAT128_TYPE("betaMin",MDPSimulator<__float128>,__float128,t._betaMin)
    REGISTER_FLOAT128_TYPE("alphaBetaInc",MDPSimulator<__float128>,__float128,t._alphaBetaInc)
    REGISTER_FLOAT128_TYPE("useKEE",MDPSimulator<__float128>,bool,t._useKEE)
  }
};
template <>
struct MDPSimulatorTraits<mpfr::mpreal>
{
  typedef mpfr::mpreal T;
  DECL_MAP_TYPES_T
  static void registerOptions(Options& ops)
  {
    REGISTER_MPFR_TYPE("tolK",MDPSimulator<mpfr::mpreal>,mpfr::mpreal,t._tolK)
    REGISTER_MPFR_TYPE("minDt",MDPSimulator<mpfr::mpreal>,mpfr::mpreal,t._minDt)
    REGISTER_MPFR_TYPE("betaMin",MDPSimulator<mpfr::mpreal>,mpfr::mpreal,t._betaMin)
    REGISTER_MPFR_TYPE("alphaBetaInc",MDPSimulator<mpfr::mpreal>,mpfr::mpreal,t._alphaBetaInc)
    REGISTER_MPFR_TYPE("useKEE",MDPSimulator<mpfr::mpreal>,bool,t._useKEE)
  }
};
//MDPSimulator
template <typename T>
MDPSimulator<T>::MDPSimulator(const ArticulatedBody& body,Options& ops,const Vec3T& g,MDP_SIMULATOR_MODE mode):MDP<T>(body,ops),_mode(mode),_warmStart(false)
{
  _a0=concat<Vec>(Vec3T::Zero(),-g);
  sizeType N=body.nrDOF();
  _DdqHatDq.resize(N,N);
  _DdqHatDdq.resize(N,N);
  _Dqdq2.resize(N*2,N*2);
  _Dqdq3.resize(N*2,N*2);
  _Dqdq4.resize(N*2,N*2);
  _Dtau2.resize(N*2,N);
  _Dtau3.resize(N*2,N);
  _Dtau4.resize(N*2,N);
  _info.reset(body,Vec::Zero(N),Vec::Zero(N),Vec::Zero(N));
  _IMCustom.setZero(6,body.nrJ()*6);

  if(!ops.hasType<MDPSimulator<T>>())
    MDPSimulatorTraits<T>::registerOptions(ops);
  reset(ops);
}
template <typename T>
MDPSimulator<T>::MDPSimulator(const MDPSimulator<T>& other):MDP<T>(other)
{
  operator=(other);
}
template <typename T>
MDPSimulator<T>& MDPSimulator<T>::operator=(const MDPSimulator<T>& other)
{
  MDP<T>::operator=(other);
  _IMCustom=other._IMCustom;
  _wrench=other._wrench;
  _PDTarget=other._PDTarget;
  _mode=other._mode;
  _warmStart=other._warmStart;
  _a0=other._a0;

  sizeType N=other._body.nrDOF();
  _DdqHatDq.resize(N,N);
  _DdqHatDdq.resize(N,N);
  _Dqdq2.resize(N*2,N*2);
  _Dqdq3.resize(N*2,N*2);
  _Dqdq4.resize(N*2,N*2);
  _Dtau2.resize(N*2,N);
  _Dtau3.resize(N*2,N);
  _Dtau4.resize(N*2,N);

  Options ops;
  ops.setOptions(const_cast<MDPSimulator<T>*>(&other));
  reset(ops);
  return *this;
}
template <typename T>
void MDPSimulator<T>::reset(Options& ops)
{
  MDP<T>::reset(ops);
  _tolK=1E-6f;
  _minDt=1E-5f;
  _betaMin=1e-4f;
  _alphaBetaInc=1.5f;
  _useKEE=true;
  ops.setOptions(this);
}
template <typename T>
typename MDPSimulator<T>::WrenchConstructor MDPSimulator<T>::getWrenchConstructor() const
{
  return _wrench;
}
template <typename T>
void MDPSimulator<T>::setWrenchConstructor(WrenchConstructor wrench)
{
  _wrench=wrench;
  if(std::dynamic_pointer_cast<C2GranularWrenchConstructor<T>>(wrench))
    calcIMCustom(std::dynamic_pointer_cast<C2GranularWrenchConstructor<T>>(wrench)->_externalForces);
  else if(std::dynamic_pointer_cast<C2EnvWrenchConstructor<T>>(wrench))
    calcIMCustom(std::dynamic_pointer_cast<C2EnvWrenchConstructor<T>>(wrench)->_externalForces);
  else _IMCustom.setZero(6,_body.nrJ()*6);
}
template <typename T>
void MDPSimulator<T>::setPDTarget(std::shared_ptr<PDTarget> PDTarget)
{
  _PDTarget=PDTarget;
  if(_PDTarget)
    _PDTarget->reset();
}
template <typename T>
void MDPSimulator<T>::writeVTKSeq(const std::string& path,T horizon,T interval,T dt,std::shared_ptr<PDTarget> PDTarget,const std::map<std::string,std::set<sizeType>>* jointMask)
{
  recreate(path);
  _PDTarget=PDTarget;
  T betaInit=_betaMin,lastDt=dt;
  sizeType id=0,N=_body.nrDOF();
  sizeType intervalFrm=sizeType(interval/dt);
  std::vector<sizeType> iter,iterG;
  std::vector<scalarD> frametime,kinetic;
  std::vector<Vec,Eigen::aligned_allocator<Vec>> qss;
  ASSERT_MSGV(intervalFrm>0,"IntervalFrm(%ld)<=0 is not allowed!",intervalFrm)
  Vec s=concat<Vec>(_info._qM,_info._dqM);
  for(T t=0; t<horizon; t+=dt) {
    TBEG("");
    s=step(dt,NULL,s,NULL,NULL,&betaInit);
    iter.push_back(_iter);
    iterG.push_back(_iterG);
    frametime.push_back(TENDV());
    _info.reset(_body,s.segment(0,N));
    kinetic.push_back(std::to_double(s.segment(N,N).dot(_info.getH(_body)*s.segment(N,N))/2));
    qss.push_back(s.segment(0,N));
    if(s.size()==0) {
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
      return (scalarD)std::to_double(in);
    }) << std::endl;
  }
}
template <typename T>
void MDPSimulator<T>::writeVTK(const Vec& q,const std::string& path,const std::map<std::string,std::set<sizeType>>* jointMask) const
{
  ASSERT_MSGV(endsWith(path,".vtk"),"writeVTK path(%s) must ends with .vtk!",path.c_str())
  NEArticulatedGradientInfo<T> info(_body,q);
  //std::cout << q.transpose().size() << " " << _body.nrDOF() << std::endl;
  Mat3Xd trans=info.getTrans().template cast<scalarD>();
  if(!jointMask || jointMask->size()==0)
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
void MDPSimulator<T>::setState(const Vec& q,const Vec& dq)
{
  _info.reset(_body,q,dq);
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::step(T dt,const Vec* tau0,const Vec& qdq,DMat* Dqdq,DMat* Dtau,T* betaInit)
{
  sizeType N=_Dtau2.cols();
  if(tau0) {
    ASSERT_MSG(tau0->size()==N,"tau0 size incorrect, should be #DOF!")
  }
  if(_mode==INVERSE_I || _mode==INVERSE_LF || _mode==INVERSE_BACKWARD_RK1F || _mode==BACKWARD_RK1F) {
    ASSERT_MSG(qdq.size()==N*3,"qdq size incorrect, should be #DOF*3!")
  } else {
    ASSERT_MSG(qdq.size()==N*2,"qdq size incorrect, should be #DOF*2!")
  }
  //ASSERT_MSG(Dqdq->rows()==N*2 && Dqdq->cols()==N*2,"Dqdq size incorrect, should be (#DOF*2)X(#DOF*2)")
  //ASSERT_MSG(Dtau->rows()==N*2 && Dtau->cols()==N,"Dtau size incorrect, should be (#DOF*2)X(#DOF*2)")
  if(_mode==INVERSE_I || _mode==INVERSE_LF) {
    if(Dqdq)
      Dqdq->resize(N,N*3);
    if(Dtau)
      Dtau->resize(0,0);
  } else if(_mode==INVERSE_BACKWARD_RK1F || _mode==BACKWARD_RK1F) {
    if(Dqdq)
      Dqdq->resize(N,N*3);
    if(Dtau)
      Dtau->resize(N,N);
  } else if(_mode==FORWARD_LF) {
    if(Dqdq)
      Dqdq->resize(N,N*2);
    if(Dtau)
      Dtau->resize(N,N);
  } else {
    if(Dqdq)
      Dqdq->resize(N*2,N*2);
    if(Dtau)
      Dtau->resize(N*2,N);
  }
  return step(dt,this->template mapCV<Vec>(tau0),this->template mapCV<Vec>(qdq),this->template mapM<DMat>(Dqdq),this->template mapM<DMat>(Dtau),betaInit);
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::step(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau,T* betaInit)
{
  Vec ret=Vec::Zero(0);
  if(_mode==INVERSE_I)
    ret=inverseI(qdq,Dqdq);
  else if(_mode==FORWARD_I)
    ret=forwardI(tau0,qdq,Dqdq,Dtau);
  else if(_mode==INVERSE_LF)
    ret=inverseLF(qdq,Dqdq);
  else if(_mode==FORWARD_LF)
    ret=forwardLF(tau0,qdq,Dqdq,Dtau);
  else if(_mode==INVERSE_BACKWARD_RK1F)
    ret=inverseBackwardRK1F(dt,tau0,qdq,Dqdq,Dtau);
  else if(_mode==BACKWARD_RK1F)
    ret=backwardRK1F(dt,tau0,qdq,Dqdq,Dtau);
  else if(_mode==FORWARD_RK1F)
    ret=forwardRK1F(dt,tau0,qdq,Dqdq,Dtau);
  else if(_mode==FORWARD_RK1I)
    ret=forwardRK1I(dt,tau0,qdq,Dqdq,Dtau);
  else if(_mode==FORWARD_RK2I)
    ret=forwardRK2I(dt,tau0,qdq,Dqdq,Dtau);
  else if(_mode==FORWARD_RK4I)
    ret=forwardRK4I(dt,tau0,qdq,Dqdq,Dtau);
  else if(_mode==NMDP_PGM || _mode==NMDP_ZOPGM) {
    ret=qdq;
    if(!stepNMDPPGMAdaptive(dt,tau0,mapV(ret),betaInit))
      ret.resize(0);
  } else {
    ASSERT(false)
  }
  if(_PDTarget)
    _PDTarget->advance(std::to_double(dt));
  return ret;
}
template <typename T>
typename MDPSimulator<T>::Vecss MDPSimulator<T>::stepBatched(T dt,const Vecss* tau0,const Vecss& qdq,DMatss* Dqdq,DMatss* Dtau)
{
  if(tau0) {
    ASSERT_MSG(tau0->size()==qdq.size(),"tau0.size()!=qdq.size()!")
  }
  if(Dqdq && Dqdq->size()!=qdq.size())
    Dqdq=NULL;
  if(Dtau && Dtau->size()!=qdq.size())
    Dtau=NULL;
  Vecss o(qdq.size());
  std::vector<std::shared_ptr<MDPSimulator<T>>> workers(OmpSettings::getOmpSettings().nrThreads());
  for(sizeType i=0; i<(sizeType)workers.size(); i++) {
    workers[i].reset(new MDPSimulator<T>(*this));
    workers[i]->_warmStart=false;
  }
  ASSERT_MSG(_mode!=NMDP_PGM && _mode!=NMDP_ZOPGM,"stepBatched does not support NMDP_PGM or NMDP_ZOPGM!")
  mpfr_prec_t prec=mpfr_get_default_prec();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)qdq.size(); i++) {
    mpfr_set_default_prec(prec);
    if(_mode==INVERSE_I)
      o[i]=workers[omp_get_thread_num()]->step(dt,tau0?&(tau0->at(i)):NULL,qdq[i],Dqdq?&(Dqdq->at(i)):NULL,Dtau?&(Dtau->at(i)):NULL,NULL);
    else
      o[i]=workers[omp_get_thread_num()]->step(dt,tau0?&(tau0->at(i)):NULL,qdq[i],Dqdq?&(Dqdq->at(i)):NULL,Dtau?&(Dtau->at(i)):NULL,NULL);
  }
  return o;
}
template <typename T>
void MDPSimulator<T>::debug(T dt,sizeType nrIter,T scale)
{
  T D;
  MatT Dqdq,Dtau;
  Vec s,s2,qdq,qdq2;
  Vec deltas,deltatau,tau0,tau02;
  DEFINE_NUMERIC_DELTA_T(T)
  std::vector<std::string> modes;
  modes.push_back("INVERSE_I");
  modes.push_back("FORWARD_I");
  modes.push_back("INVERSE_LF");
  modes.push_back("FORWARD_LF");
  modes.push_back("INVERSE_BACKWARD_RK1F");
  modes.push_back("BACKWARD_RK1F");
  modes.push_back("FORWARD_RK1F");
  modes.push_back("FORWARD_RK1I");
  modes.push_back("FORWARD_RK2I");
  modes.push_back("FORWARD_RK4I");
  for(sizeType i=0; i<nrIter; i++) {
    INFOV("-------------------------------------------------------------DebugMDPSimulator: scale=%s!",std::to_string(scale).c_str())
    //test
    for(sizeType m=0; m<(sizeType)modes.size(); m++) {
      //randomize MDP
      if(scale>0)
        MDP<T>::randomize(scale,m!=FORWARD_LF);
      else _externalWrench.clear();
      std::shared_ptr<DebugWrenchConstructor<T>> debugWrench(new DebugWrenchConstructor<T>(_externalWrench));
      setWrenchConstructor(debugWrench);
      setMode((MDP_SIMULATOR_MODE)m);
      if(m==INVERSE_I || m==INVERSE_LF || m==INVERSE_BACKWARD_RK1F || m==BACKWARD_RK1F) {
        s.setRandom(_info._qM.size()*3);
        deltas.setRandom(_info._qM.size()*3);
        Dqdq.setRandom(_info._qM.size()*1,_info._qM.size()*3);
      } else if(m==FORWARD_LF) {
        s.setRandom(_info._qM.size()*2);
        deltas.setRandom(_info._qM.size()*2);
        Dqdq.setRandom(_info._qM.size()*1,_info._qM.size()*2);
      } else {
        s.setRandom(_info._qM.size()*2);
        deltas.setRandom(_info._qM.size()*2);
        Dqdq.setRandom(_info._qM.size()*2,_info._qM.size()*2);
      }
      D=std::max<T>(DELTA,1e-20f);
      s2=s+deltas*D;
      tau0.setRandom(_info._qM.size());
      deltatau.setRandom(_info._qM.size());
      Dtau.setRandom(Dqdq.rows(),_info._qM.size());
      tau02=tau0+deltatau*D;

      _warmStart=true;
      _lastW.setZero(0);
      qdq=step(dt,NULL,s,&Dqdq,NULL);
      qdq2=step(dt,NULL,s2,NULL,NULL);
      DEBUG_GRADIENT(modes[m]+"-Dqdq",std::sqrt((Dqdq*deltas).squaredNorm()),std::sqrt((Dqdq*deltas-(qdq2-qdq)/D).squaredNorm()))

      _warmStart=true;
      _lastW.setZero(0);
      qdq=step(dt,&tau0,s,&Dqdq,NULL);
      qdq2=step(dt,&tau0,s2,NULL,NULL);
      DEBUG_GRADIENT(modes[m]+"-tau0-Dqdq",std::sqrt((Dqdq*deltas).squaredNorm()),std::sqrt((Dqdq*deltas-(qdq2-qdq)/D).squaredNorm()))

      if(m!=INVERSE_I && m!=INVERSE_LF) {
        _warmStart=true;
        _lastW.setZero(0);
        qdq=step(dt,&tau0,s,NULL,&Dtau);
        qdq2=step(dt,&tau02,s,NULL,NULL);
        DEBUG_GRADIENT(modes[m]+"-Dtau",std::sqrt((Dtau*deltatau).squaredNorm()),std::sqrt((Dtau*deltatau-(qdq2-qdq)/D).squaredNorm()))
      }

      if(m==INVERSE_I) {
        setMode(FORWARD_I);
        Vec ddq=step(dt,&tau0,s.segment(0,_info._qM.size()*2),NULL,NULL).segment(_info._qM.size(),_info._qM.size());
        setMode(INVERSE_I);
        Vec tau0Ref=step(dt,NULL,concat<Vec>(s.segment(0,_info._qM.size()*2),ddq),NULL,NULL);
        DEBUG_GRADIENT("INVERSE/FORWARD_I",std::sqrt(tau0.squaredNorm()),std::sqrt((tau0-tau0Ref).squaredNorm()))
      }

      if(m==INVERSE_LF) {
        setMode(FORWARD_LF);
        Vec ddq=step(dt,&tau0,s.segment(0,_info._qM.size()*2),NULL,NULL);
        setMode(INVERSE_LF);
        Vec tau0Ref=step(dt,NULL,concat<Vec>(s.segment(0,_info._qM.size()*2),ddq),NULL,NULL);
        DEBUG_GRADIENT("INVERSE/FORWARD_LF",std::sqrt(tau0.squaredNorm()),std::sqrt((tau0-tau0Ref).squaredNorm()))
      }

      if(m==INVERSE_BACKWARD_RK1F) {
        setMode(BACKWARD_RK1F);
        Vec CF=step(dt,&tau0,s,NULL,NULL);
        setMode(INVERSE_BACKWARD_RK1F);
        Vec CI=step(dt,&tau0,s,NULL,NULL);
        DEBUG_GRADIENT("INVERSE/FORWARD_BACKWARD_RK1F",std::sqrt(CF.squaredNorm()),std::sqrt((-_info._HM*CF/dt-CI).squaredNorm()))
      }

      //batched
      DMatss Dqdqss,Dtauss;
      Vecss sss,s2ss,tau0ss,tau02ss,qdqss,qdq2ss;
      for(sizeType b=0,bs=RandEngine::randI(5,10); b<bs; b++) {
        sss.push_back(Vec::Random(s.size()));
        s2ss.push_back(sss.back()+deltas*D);
        tau0ss.push_back(Vec::Random(tau0.size()));
        tau02ss.push_back(tau0ss.back()+deltatau*D);
        Dqdqss.push_back(DMat());
        Dtauss.push_back(DMat());
      }

      _warmStart=true;
      qdqss=stepBatched(dt,NULL,sss,&Dqdqss,NULL);
      qdq2ss=stepBatched(dt,NULL,s2ss,NULL,NULL);
      for(sizeType b=0; b<(sizeType)sss.size(); b++) {
        DEBUG_GRADIENT(modes[m]+"-Batched-Dqdq",std::sqrt((Dqdqss[b]*deltas).squaredNorm()),std::sqrt((Dqdqss[b]*deltas-(qdq2ss[b]-qdqss[b])/D).squaredNorm()))
      }

      _warmStart=true;
      qdqss=stepBatched(dt,&tau0ss,sss,&Dqdqss,NULL);
      qdq2ss=stepBatched(dt,&tau0ss,s2ss,NULL,NULL);
      for(sizeType b=0; b<(sizeType)sss.size(); b++) {
        DEBUG_GRADIENT(modes[m]+"-Batched-tau0-Dqdq",std::sqrt((Dqdqss[b]*deltas).squaredNorm()),std::sqrt((Dqdqss[b]*deltas-(qdq2ss[b]-qdqss[b])/D).squaredNorm()))
      }

      if(m!=INVERSE_I && m!=INVERSE_LF) {
        _warmStart=true;
        qdqss=stepBatched(dt,&tau0ss,sss,NULL,&Dtauss);
        qdq2ss=stepBatched(dt,&tau02ss,sss,NULL,NULL);
        for(sizeType b=0; b<(sizeType)sss.size(); b++) {
          DEBUG_GRADIENT(modes[m]+"-Batched-Dtau",std::sqrt((Dtauss[b]*deltatau).squaredNorm()),std::sqrt((Dtauss[b]*deltatau-(qdq2ss[b]-qdqss[b])/D).squaredNorm()))
        }
      }
    }
  }
}
template <typename T>
void MDPSimulator<T>::debugNMDP(T dt,sizeType nrIter,T scale)
{
  T KVal,KVal2;
  bool useKEE=_useKEE;
  sizeType N=_body.nrDOF();
  MatT DGDTheta=MatT::Zero(N,N),DDKDDTheta=MatT::Zero(N,N),DGDw;
  Vec GVal,GVal2,DKDTheta=Vec::Zero(N),DKDTheta2=Vec::Zero(N),delta,qMqMM,w,w2,tau0;
  MDP_SIMULATOR_MODE modeTmp=_mode;
  _mode=NMDP_PGM;
  DEFINE_NUMERIC_DELTA_T(T)
  INFOV("-------------------------------------------------------------DebugMDPSimulator-NMDP: scale=%s!",std::to_string(scale).c_str())
  for(sizeType i=0; i<nrIter; i++) {
    qMqMM=Vec::Random(N*2);
    tau0=Vec::Random(N);
    _info.reset(_body,_info._qM);
    (*_wrench)(_externalWrench,_info);
    _PDTarget.reset(new PDTarget(PDTarget::Vec::Random(N),PDTarget::Vec::Random(N),PDTarget::Vec::Random(2*N)));
    for(sizeType KEEVal=0; KEEVal<2; KEEVal++) {
      _useKEE=KEEVal;
      w.setRandom(MDP<T>::nW());
      DGDw.resize(N,MDP<T>::nW());
      KVal=K(dt,mapCV(qMqMM),mapV(DKDTheta),mapM(DDKDDTheta));
      GVal=G(dt,mapCV(tau0),mapCV(qMqMM),mapM(DGDTheta),&w,mapM(DGDw));
      //DGDw
      delta=Vec::Random(MDP<T>::nW());
      w2=w+delta*DELTA;
      GVal2=G(dt,mapCV(tau0),mapCV(qMqMM),mapM((MatT*)NULL),&w2,mapM((MatT*)NULL));
      DEBUG_GRADIENT("NMDP-DGDw",std::sqrt((DGDw*delta).squaredNorm()),std::sqrt((DGDw*delta-(GVal2-GVal)/DELTA).squaredNorm()))
      //K
      delta=Vec::Random(N);
      _info.reset(_body,_info._qM+delta*DELTA);
      (*_wrench)(_externalWrench,_info);
      KVal2=K(dt,mapCV(qMqMM),mapV(DKDTheta2),mapM((MatT*)NULL));
      DEBUG_GRADIENT("NMDP-DKDTheta-useKEE="+std::to_string(_useKEE?"1":"0"),DKDTheta.dot(delta),DKDTheta.dot(delta)-(KVal2-KVal)/DELTA)
      //DGDTheta
      GVal2=G(dt,mapCV(tau0),mapCV(qMqMM),mapM((MatT*)NULL),&w,mapM((MatT*)NULL));
      DEBUG_GRADIENT("NMDP-DGDTheta",std::sqrt((DGDTheta*delta).squaredNorm()),std::sqrt((DGDTheta*delta-(GVal2-GVal)/DELTA).squaredNorm()))
    }
  }
  _mode=modeTmp;
  _useKEE=useKEE;
}
template <typename T>
void MDPSimulator<T>::setMode(MDP_SIMULATOR_MODE m)
{
  _mode=m;
}
template <typename T>
MDP_SIMULATOR_MODE MDPSimulator<T>::getMode() const
{
  return _mode;
}
template <typename T>
void MDPSimulator<T>::setWarmStart(bool warm)
{
  _warmStart=warm;
}
template <typename T>
void MDPSimulator<T>::begLog(const std::string& path)
{
  _log.reset(new std::ofstream(path));
}
template <typename T>
void MDPSimulator<T>::endLog()
{
  _log=NULL;
}
template <typename T>
void MDPSimulator<T>::log(T tw,T ts,T tne)
{
  if(_log)
    *_log << "wrench=" << tw << ",solve=" << ts << ",newtonEuler=" << tne << std::endl;
}
template <typename T>
typename MDPSimulator<T>::Vec3T MDPSimulator<T>::g() const
{
  return -_a0.template segment<3>(3);
}
//helper
template <typename T>
void MDPSimulator<T>::calcIMCustom(const std::vector<EndEffectorBounds>& EE)
{
  Mat6X3T L;
  _IMCustom.setZero(6,6*_body.nrJ());
  Mat6T P;
  P.setZero();
  P.template block<3,3>(0,3).setIdentity();
  P.template block<3,3>(3,0).setIdentity();
  for(sizeType i=0; i<(sizeType)EE.size(); i++) {
    const EndEffectorBounds& ee=EE[i];
    L.template block<3,3>(0,0).setIdentity();
    L.template block<3,3>(3,0)=cross<T>(ee._localPos.template cast<T>());
    _IMCustom.template block<6,6>(0,6*ee.jointId())+=P*L*L.transpose()*P;
  }
}
//forward/backward impulse
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::inverseI(VecCM qdqddq,MatTM jac)
{
  sizeType N=_info._qM.size();
  //MDP
  if(!_warmStart)
    _lastW.resize(0);
  _info.NEArticulatedGradientInfoMap<T>::reset(_body,VecCM(&(qdqddq.coeff(0)),N),VecCM(&(qdqddq.coeff(N)),N),VecCM(&(qdqddq.coeff(N*2)),N));
  TBEG();
  (*_wrench)(_externalWrench,_info);
  T tw=TENDV();
  TBEG();
  bool modified=MDP<T>::solveMDPQPF(jac.data()?this->template mapM<MatT>(_DdqHatDq):this->template mapM<MatT>((MatT*)NULL),
                                    jac.data()?this->template mapM<MatT>(_DdqHatDdq):this->template mapM<MatT>((MatT*)NULL));
  T ts=TENDV();

  //Newton-Euler
  TBEG();
  if(modified)
    _info._dqM=_dqMDP;
  if(jac.data()) {
    _info.RNEADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);
    _info.calcH(_body);
    jac.block(0,0,N,N)=_info._DtauDqM+_info._DtauDdqM*_DdqHatDq;
    jac.block(0,N,N,N)=_info._DtauDdqM*_DdqHatDdq;
    jac.block(0,N*2,N,N)=_info._HM;
  } else _info.RNEA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);
  T tne=TENDV();
  log(tw,ts,tne);
  return _info._tauM;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::forwardI(VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau)
{
  sizeType N=_info._qM.size();
  //MDP
  if(!_warmStart)
    _lastW.resize(0);
  _info.NEArticulatedGradientInfoMap<T>::reset(_body,VecCM(&(qdq.coeff(0)),N),VecCM(&(qdq.coeff(N)),N));
  TBEG();
  (*_wrench)(_externalWrench,_info);
  T tw=TENDV();
  TBEG();
  MDP<T>::solveMDPQPF(Dqdq.data()?this->template mapM<MatT>(_DdqHatDq):this->template mapM<MatT>((MatT*)NULL),
                      Dqdq.data()?this->template mapM<MatT>(_DdqHatDdq):this->template mapM<MatT>((MatT*)NULL));
  T ts=TENDV();

  //Newton-Euler
  TBEG();
  _info._dqM=_dqMDP;
  Vec ctrl=getCtrl(tau0,qdq);
  if(Dqdq.data()) {
    _info.ABADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),mapCV(ctrl),true);
    Dqdq.block(0,0,N,N)=_DdqHatDq;
    Dqdq.block(0,N,N,N)=_DdqHatDdq;
    Dqdq.block(N,0,N,N)=_info.getDddqDq()+_info.getDddqDdq()*_DdqHatDq;
    Dqdq.block(N,N,N,N)=_info.getDddqDdq()*_DdqHatDdq;
  } else _info.ABA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),mapCV(ctrl));
  if(Dtau.data()) {
    _info.calcHInvH(_body);
    Dtau.block(0,0,N,N).setZero();
    Dtau.block(N,0,N,N)=_info._invHM;
  }
  T tne=TENDV();
  log(tw,ts,tne);
  return concat<Vec>(_info._dqM,_info._ddqM);
}
//forward/backward limiting force
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::inverseLF(VecCM qdqddq,MatTM jac)
{
  sizeType N=_info._qM.size();
  //MDP
  if(!_warmStart)
    _lastW.resize(0);
  _info.NEArticulatedGradientInfoMap<T>::reset(_body,VecCM(&(qdqddq.coeff(0)),N),VecCM(&(qdqddq.coeff(N)),N),VecCM(&(qdqddq.coeff(N*2)),N));
  TBEG();
  (*_wrench)(_externalWrench,_info);
  T tw=TENDV();
  TBEG();
  bool modified=MDP<T>::solveMDPLPI(jac.data()?this->template mapM<MatT>(_DdqHatDq):this->template mapM<MatT>((MatT*)NULL),
                                    jac.data()?this->template mapM<MatT>(_DdqHatDdq):this->template mapM<MatT>((MatT*)NULL));
  T ts=TENDV();

  //Newton-Euler
  TBEG();
  if(jac.data()) {
    _info.RNEADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);
    _info.calcH(_body);
    jac.block(0,0,N,N)=_info._DtauDqM-_DdqHatDq;
    jac.block(0,N,N,N)=_info._DtauDdqM-_DdqHatDdq;
    jac.block(0,N*2,N,N)=_info._HM;
  } else _info.RNEA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);
  if(modified)
    _info._tauM-=_deltaDq;
  T tne=TENDV();
  log(tw,ts,tne);
  return _info._tauM;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::forwardLF(VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau)
{
  sizeType N=_info._qM.size();
  //MDP
  if(!_warmStart)
    _lastW.resize(0);
  _info.NEArticulatedGradientInfoMap<T>::reset(_body,VecCM(&(qdq.coeff(0)),N),VecCM(&(qdq.coeff(N)),N));
  TBEG();
  (*_wrench)(_externalWrench,_info);
  T tw=TENDV();
  TBEG();
  bool modified=MDP<T>::solveMDPLPF(Dqdq.data()?this->template mapM<MatT>(_DdqHatDq):this->template mapM<MatT>((MatT*)NULL),
                                    Dqdq.data()?this->template mapM<MatT>(_DdqHatDdq):this->template mapM<MatT>((MatT*)NULL));
  T ts=TENDV();

  //Newton-Euler
  TBEG();
  Vec ctrl=getCtrl(tau0,qdq);
  if(Dqdq.data()) {
    _info.ABADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),mapCV(ctrl),true);
    Dqdq.block(0,0,N,N)=_info.getDddqDq()+_DdqHatDq;
    Dqdq.block(0,N,N,N)=_info.getDddqDdq()+_DdqHatDdq;
  } else _info.ABA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),mapCV(ctrl));
  if(Dtau.data()) {
    _info.calcHInvH(_body);
    Dtau.block(0,0,N,N)=_info._invHM;
  }
  if(modified)
    _info._ddqM+=_deltaDq;
  T tne=TENDV();
  log(tw,ts,tne);
  return _info._ddqM;
}
//forward/backward force+RK1
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::inverseBackwardRK1F(T dt,VecCM tau0,VecCM qdqNextdq,MatTM DqdqNextdq,MatTM Dtau)
{
  sizeType N=_info._qM.size();
  //Newton-Euler
  _info.NEArticulatedGradientInfoMap<T>::reset(_body,VecCM(&(qdqNextdq.coeff(N*0)),N),VecCM(&(qdqNextdq.coeff(N*1)),N));
  TBEG();
  (*_wrench)(_externalWrench,_info);
  T tw=TENDV();
  TBEG();
  for(ExternalWrench<T>& e:_externalWrench)
    e*=dt;
  if(DqdqNextdq.data())
    _info.ABADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),tau0,true);
  else _info.ABA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),tau0);
  _info._dqM=qdqNextdq.segment(N*2,N)+_info._ddqM*dt;
  T tne=TENDV();

  //MDP
  TBEG();
  if(!_warmStart)
    _lastW.resize(0);
  MDP<T>::solveMDPQPI(DqdqNextdq.data()||Dtau.data()?this->template mapM<MatT>(_DdqHatDq):this->template mapM<MatT>((MatT*)NULL),
                      DqdqNextdq.data()||Dtau.data()?this->template mapM<MatT>(_DdqHatDdq):this->template mapM<MatT>((MatT*)NULL));
  if(DqdqNextdq.data()) {
    DqdqNextdq.block(0,N*0,N,N)=-_DdqHatDq/dt-_DdqHatDdq*_info.getDddqDq();
    DqdqNextdq.block(0,N*1,N,N)=-_DdqHatDdq*_info.getDddqDdq();
    DqdqNextdq.block(0,N*2,N,N)=-_DdqHatDdq/dt;
  }
  T ts=TENDV();
  log(tw,ts,tne);

  //inverse
  _info._dqM=qdqNextdq.segment(N*1,N);
  _info._ddqM=(qdqNextdq.segment(N*1,N)-qdqNextdq.segment(N*2,N))/dt;
  if(DqdqNextdq.data()) {
    _info.RNEADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);
    _info.calcH(_body);
    DqdqNextdq.block(0,N*0,N,N)+=_info._DtauDqM;
    DqdqNextdq.block(0,N*1,N,N)+=_info._DtauDdqM+_info._HM/dt;
    DqdqNextdq.block(0,N*2,N,N)-=_info._HM/dt;
  } else _info.RNEA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);
  if(Dtau.data()) {
    _info.calcHInvH(_body);
    Dtau=-_DdqHatDdq*_info._invHM-MatT::Identity(N,N);
  }
  if(tau0.data())
    return _info._tauM-_deltaDq/dt-tau0;
  else return _info._tauM-_deltaDq/dt;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::backwardRK1F(T dt,VecCM tau0,VecCM qdqNextdq,MatTM DqdqNextdq,MatTM Dtau)
{
  sizeType N=_info._qM.size();
  //Newton-Euler
  _info.NEArticulatedGradientInfoMap<T>::reset(_body,VecCM(&(qdqNextdq.coeff(N*0)),N),VecCM(&(qdqNextdq.coeff(N*1)),N));
  TBEG();
  (*_wrench)(_externalWrench,_info);
  T tw=TENDV();
  TBEG();
  for(ExternalWrench<T>& e:_externalWrench)
    e*=dt;
  if(DqdqNextdq.data())
    _info.ABADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),tau0,true);
  else _info.ABA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),tau0);
  _info._dqM=qdqNextdq.segment(N*2,N)+_info._ddqM*dt;
  T tne=TENDV();

  //MDP
  TBEG();
  if(!_warmStart)
    _lastW.resize(0);
  MDP<T>::solveMDPQPF(DqdqNextdq.data()||Dtau.data()?this->template mapM<MatT>(_DdqHatDq):this->template mapM<MatT>((MatT*)NULL),
                      DqdqNextdq.data()||Dtau.data()?this->template mapM<MatT>(_DdqHatDdq):this->template mapM<MatT>((MatT*)NULL));
  if(DqdqNextdq.data()) {
    DqdqNextdq.block(0,N*0,N,N)=_DdqHatDq+_DdqHatDdq*_info.getDddqDq()*dt;
    DqdqNextdq.block(0,N*1,N,N)=_DdqHatDdq*_info.getDddqDdq()*dt-MatT::Identity(N,N);
    DqdqNextdq.block(0,N*2,N,N)=_DdqHatDdq;
  }
  if(Dtau.data()) {
    _info.calcHInvH(_body);
    Dtau=_DdqHatDdq*_info._invHM*dt;
  }
  T ts=TENDV();
  log(tw,ts,tne);
  return _dqMDP-qdqNextdq.segment(N,N);
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::forwardRK1F(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau)
{
  sizeType N=_info._qM.size();
  //Newton-Euler
  _info.NEArticulatedGradientInfoMap<T>::reset(_body,VecCM(&(qdq.coeff(0)),N),VecCM(&(qdq.coeff(N)),N));
  TBEG();
  (*_wrench)(_externalWrench,_info);
  T tw=TENDV();
  TBEG();
  for(ExternalWrench<T>& e:_externalWrench)
    e*=dt;
  Vec ctrl=getCtrl(tau0,qdq);
  if(Dqdq.data())
    _info.ABADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),mapCV(ctrl),true);
  else _info.ABA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),mapCV(ctrl));
  _info._dqM+=_info._ddqM*dt;
  T tne=TENDV();

  //MDP
  TBEG();
  if(!_warmStart)
    _lastW.resize(0);
  MDP<T>::solveMDPQPF(Dqdq.data()||Dtau.data()?this->template mapM<MatT>(_DdqHatDq):this->template mapM<MatT>((MatT*)NULL),
                      Dqdq.data()||Dtau.data()?this->template mapM<MatT>(_DdqHatDdq):this->template mapM<MatT>((MatT*)NULL));
  if(Dqdq.data()) {
    _DdqHatDq+=_DdqHatDdq*_info.getDddqDq()*dt;
    _DdqHatDdq+=_DdqHatDdq*_info.getDddqDdq()*dt;
    Dqdq.block(0,0,N,N)=_DdqHatDq*dt+MatT::Identity(N,N);
    Dqdq.block(0,N,N,N)=_DdqHatDdq*dt;
    Dqdq.block(N,0,N,N)=_DdqHatDq;
    Dqdq.block(N,N,N,N)=_DdqHatDdq;
  }
  if(Dtau.data()) {
    _info.calcHInvH(_body);
    Dtau.block(N,0,N,N)=_DdqHatDdq*_info._invHM*dt;
    Dtau.block(0,0,N,N)=Dtau.block(N,0,N,N)*dt;
  }
  _info._dqM=_dqMDP;
  _info._qM+=_info._dqM*dt;
  T ts=TENDV();
  log(tw,ts,tne);
  return concat<Vec>(_info._qM,_info._dqM);
}
//forward impulse+RK1/2/4
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::forwardRK1I(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau)
{
  Vec k1=forwardI(tau0,qdq,Dqdq,Dtau);
  if(Dtau.data())
    Dtau*=dt;
  if(Dqdq.data())
    Dqdq=Dqdq*dt+MatT::Identity(Dqdq.rows(),Dqdq.cols());
  return qdq+k1*dt;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::forwardRK2I(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau)
{
  Vec tmp;
  Vec k1=forwardI(tau0,qdq,Dqdq,Dtau);
  Vec k2=forwardI(tau0,this->template mapCV<Vec>(tmp=qdq+k1*dt/2),(Dqdq.data()||Dtau.data())?this->template mapM<MatT>(_Dqdq2):this->template mapM<MatT>((MatT*)NULL),this->template mapM<MatT>(_Dtau2));
  if(Dtau.data()) {
    _Dtau2=_Dtau2+_Dqdq2*Dtau*(dt/2);
    Dtau=_Dtau2*dt;
  }
  if(Dqdq.data()) {
    _Dqdq2=_Dqdq2+_Dqdq2*Dqdq*(dt/2);
    Dqdq=_Dqdq2*dt+MatT::Identity(Dqdq.rows(),Dqdq.cols());
  }
  return qdq+k2*dt;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::forwardRK4I(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau)
{
  Vec tmp;
  Vec k1=forwardI(tau0,qdq,Dqdq,Dtau);
  Vec k2=forwardI(tau0,this->template mapCV<Vec>(tmp=qdq+k1*dt/2),(Dqdq.data()||Dtau.data())?this->template mapM<MatT>(_Dqdq2):this->template mapM<MatT>((MatT*)NULL),this->template mapM<MatT>(_Dtau2));
  Vec k3=forwardI(tau0,this->template mapCV<Vec>(tmp=qdq+k2*dt/2),(Dqdq.data()||Dtau.data())?this->template mapM<MatT>(_Dqdq3):this->template mapM<MatT>((MatT*)NULL),this->template mapM<MatT>(_Dtau3));
  Vec k4=forwardI(tau0,this->template mapCV<Vec>(tmp=qdq+k3*dt)  ,(Dqdq.data()||Dtau.data())?this->template mapM<MatT>(_Dqdq4):this->template mapM<MatT>((MatT*)NULL),this->template mapM<MatT>(_Dtau4));
  if(Dtau.data()) {
    _Dtau2=_Dtau2+_Dqdq2*Dtau*(dt/2);
    _Dtau3=_Dtau3+_Dqdq3*_Dtau2*(dt/2);
    _Dtau4=_Dtau4+_Dqdq4*_Dtau3*dt;
    Dtau=(Dtau+_Dtau2*2+_Dtau3*2+_Dtau4)*dt/6;
  }
  if(Dqdq.data()) {
    _Dqdq2=_Dqdq2+_Dqdq2*Dqdq*(dt/2);
    _Dqdq3=_Dqdq3+_Dqdq3*_Dqdq2*(dt/2);
    _Dqdq4=_Dqdq4+_Dqdq4*_Dqdq3*dt;
    Dqdq=(Dqdq+2*_Dqdq2+2*_Dqdq3+_Dqdq4)*dt/6+MatT::Identity(Dqdq.rows(),Dqdq.cols());
  }
  return qdq+(k1+2*k2+2*k3+k4)*dt/6;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::getCtrl(VecCM tau0,VecCM qdq) const
{
  if(!tau0.data() && !_PDTarget)
    return Vec();
  sizeType N=qdq.size()/2;
  Vec ret=Vec::Zero(N);
  if(tau0.data())
    ret=tau0;
  if(_PDTarget) {
    const Vec s=_PDTarget->s().template cast<T>();
    ret+=(s-qdq).segment(0,N).cwiseProduct(_PDTarget->_PCoef.template cast<T>());
    ret+=(s-qdq).segment(N,N).cwiseProduct(_PDTarget->_DCoef.template cast<T>());
  }
  return ret;
}
//NMDP
template <typename T>
T MDPSimulator<T>::K(T dt,VecCM qdq,VecM DKDTheta,MatTM DDKDDTheta)
{
  MatT DHDq,H;
  sizeType N=qdq.size()/2;
  Vec dqdt=(_info._qM-qdq.segment(0,N))/dt;
  //calcH
  NEArticulatedGradientInfo<T> info=_info;
  if(_useKEE) {
    H.setZero(N,N);
    _info.calcHInner(_body,mapM(H),mapCM(_IMCustom));
  } else H=_info.getH(_body);
  //calcDHDq
  if(DKDTheta.data() || DDKDDTheta.data()) {
    DHDq.resize(N,N);
    if(_useKEE)
      _info.DHDq(_body,mapM(DHDq),mapCV(dqdt),mapCM(_IMCustom));
    else _info.DHDq(_body,mapM(DHDq),mapCV(dqdt));
  }
  //assemble
  if(DKDTheta.data())
    DKDTheta=(H*dqdt)/dt+DHDq.transpose()*dqdt/2;
  if(DDKDDTheta.data())
    DDKDDTheta=H/(dt*dt);
  _info=info;
  return dqdt.dot(H*dqdt)/2;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::G(T dt,VecCM tau0,VecCM qdq,MatTM DGDTheta,Vec* w,MatTM DGDw)
{
  _iterG++;
  sizeType N=qdq.size()/2;
  _info._dqM=(_info._qM-qdq.segment(0,N))/dt;
  _info._ddqM=(_info._dqM-qdq.segment(N,N))/dt;
  if(DGDTheta.data())
    _info.RNEADerivatives(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);
  else _info.RNEA(_body,this->template mapCV<Vec6T>(_a0),this->template mapCM<Mat6XT>((Mat6XT*)NULL),true);

  (*_wrench)(_externalWrench,_info,_mode!=NMDP_ZOPGM);
  if(w->size()==0)
    w->setZero(MDP<T>::nW());
  _joints.resize(_externalWrench.size());
  ColiCM jointsM(&_joints[0],(sizeType)_externalWrench.size());
  _Sc.resize(6*(sizeType)_externalWrench.size(),N);
  MatT B=MatT::Zero(6*(sizeType)_externalWrench.size(),w->size()),DBDqw;
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); i++) {
    ExternalWrench<T>& wi=_externalWrench[i];
    //wi*=dt;
    wi._w=w->segment(offi,wi._B.cols());
    B.block(i*6,offi,wi._B.rows(),wi._B.cols())=wi._B;
    _joints[i]=wi._jointId;
    offi+=wi._B.cols();
  }
  _info.Sc(_body,jointsM,mapM(_Sc));
  Vec negBw=-B**w;
  if(DGDTheta.data()) {
    _info.DScDqT(_body,jointsM,DGDTheta,mapCV(negBw));
    DGDTheta+=_info._DtauDqM+_info._DtauDdqM/dt+_info.getH(_body)/(dt*dt);
    if(_mode!=NMDP_ZOPGM) {
      MDP<T>::DBDq(*w,DBDqw);
      DGDTheta-=_Sc.transpose()*DBDqw;
    }
  }
  if(DGDw.data())
    DGDw=-_Sc.transpose()*B;
  Vec ret=_info._tauM+_Sc.transpose()*negBw;
  if(tau0.data())
    ret-=tau0;
  if(_PDTarget) {
    const Vec s=_PDTarget->s().template cast<T>();
    ret+=(_info._qM-s.segment(0,N)).cwiseProduct(_PDTarget->_PCoef.template cast<T>());
    ret+=((_info._qM-qdq.segment(0,N))/dt-s.segment(N,N)).cwiseProduct(_PDTarget->_DCoef.template cast<T>());
    if(DGDTheta.data()) {
      DGDTheta.diagonal()+=_PDTarget->_PCoef.template cast<T>();
      DGDTheta.diagonal()+=_PDTarget->_DCoef.template cast<T>()/dt;
    }
  }
  return ret;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::stepNMDPPGM(T dt,VecCM tau0,VecCM qdq,T* betaInit)
{
  sizeType N=_body.nrDOF();
  T betaMin=_betaMin,beta=std::max<T>(betaInit?*betaInit:betaMin,betaMin),inc=_alphaBetaInc,alpha=1;
  _info.reset(_body,qdq.segment(0,N));
  MatT DDKDDTheta=MatT::Zero(N,N),DGDTheta=MatT::Zero(N,N),DGDw,H,A;
  Vec DKDTheta=Vec::Zero(N),q=qdq.segment(0,N),qNext,w=Vec::Zero(0),wNext=Vec::Zero(0),g;
  bool updateQP=true;
  //initial projection
  q=manifoldProjection(dt,tau0,qdq,mapCV(q),&w,&alpha);
  if(q.size()==0)
    return q;
  QCQPSolverMosek<T> sol;
  DGDw.resize(N,MDP<T>::nW());
  T KVal=K(dt,qdq,mapV((Vec*)NULL),mapM((MatT*)NULL)),KValNext;
  for(; _iter<_solLQP.maxIter(); _iter++) {
    //buildQP
    if(updateQP) {
      K(dt,qdq,mapV(DKDTheta),mapM(DDKDDTheta));
      G(dt,tau0,qdq,mapM(DGDTheta),&w,mapM(DGDw));
      SolveNewton<T>::solveLU(DGDTheta,DGDw,H,true);
      g=H.transpose()*DKDTheta;
      H=H.transpose()*(DDKDDTheta*H);
      if(A.size()==0) {
        A.setZero(_externalWrench.size(),MDP<T>::nW());
        for(sizeType i=0,offW=0; i<(sizeType)_externalWrench.size(); offW+=_externalWrench[i]._B.cols(),i++)
          A.block(i,offW,1,_externalWrench[i]._B.cols()).setOnes();
      }
    }
    //solveQP
    Vec dwd;
    if(sol.QCQPSolver<T>::solveQP(dwd,H,g,A,w,beta)!=QCQPSolver<T>::SOLVED) {
      if(_solLQP.callback()) {
        INFOV("qp failed (iter=%d)!",_iter+1)
      }
      q.resize(0);
      break;
    }
    if(std::sqrt(dwd.squaredNorm())<_solLQP.tolGFinal()) {
      if(_solLQP.callback()) {
        INFOV("Beta: %f",std::to_double(beta))
      }
      break;
    }
    if(dwd.size()==0) {
      beta*=inc;
      updateQP=false;
      if(beta>1/_solLQP.tolAlpha()) {
        if(_solLQP.callback()) {
          INFOV("Line-Search failed (iter=%ld)!",_iter+1)
        }
        q.resize(0);
        break;
      } else continue;
    } else wNext=w+dwd;
    //manifold projection
    qNext=manifoldProjection(dt,tau0,qdq,mapCV(q),&wNext,&alpha);
    if(qNext.size()==0)
      return qNext;
    else {
      _info.reset(_body,qNext);
      KValNext=K(dt,qdq,mapV((Vec*)NULL),mapM((MatT*)NULL));
    }
    //update line search info
    if(KValNext<KVal) {
      w=wNext;
      q=qNext;
      if(_solLQP.callback()) {
        INFOV("NMDP-Iter%d: K=%f, beta=%f!",_iter,std::to_double(KVal),std::to_double(beta))
      }
      if(KVal-KValNext<_tolK)
        break;
      KVal=KValNext;
      updateQP=true;
      beta=std::max<T>(beta/inc,betaMin);
    } else {
      updateQP=false;
      beta*=inc;
      if(beta>1/_solLQP.tolAlpha()) {
        if(_solLQP.callback()) {
          INFOV("Line-Search failed (iter=%ld)!",_iter+1)
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
bool MDPSimulator<T>::stepNMDPPGMAdaptive(T dt,VecCM tau0,VecM qdq,T* betaInit)
{
  _iter=0,_iterG=0;
  if(dt<_minDt)  //this algorithm has failed
    return false;
  sizeType N=qdq.size()/2;
  T betaInitTmp=betaInit?*betaInit:_betaMin;
  Vec ret=stepNMDPPGM(dt,tau0,mapV2CV(qdq),&betaInitTmp);
  if(ret.size()==0) {
    //subdivide: piece-1
    betaInitTmp=betaInit?*betaInit:_betaMin;
    if(!stepNMDPPGMAdaptive(dt/2,tau0,qdq,&betaInitTmp))
      return false;
    //subdivide: piece-2
    if(!stepNMDPPGMAdaptive(dt/2,tau0,qdq,&betaInitTmp))
      return false;
    if(betaInit)
      *betaInit=betaInitTmp;
  } else {
    if(betaInit)
      *betaInit=betaInitTmp;
    qdq.segment(N,N)=(ret-qdq.segment(0,N))/dt;
    qdq.segment(0,N)=ret;
  }
  return true;
}
template <typename T>
typename MDPSimulator<T>::Vec MDPSimulator<T>::manifoldProjection(T dt,VecCM tau0,VecCM qdq,VecCM init,Vec* w,T* alphaInit)
{
  sizeType N=init.size();
  Vec q=init,qNext,Gi,dir;
  MatT DGDTheta=MatT::Zero(N,N);
  T alpha=alphaInit?*alphaInit:1,inc=_alphaBetaInc,V;
  T tolLS=1;//_solLQP.tolAlpha();
  _info.reset(_body,q);
  Gi=G(dt,tau0,qdq,mapM(DGDTheta),w,mapM((MatT*)NULL));
  for(sizeType i=0; i<_solLQP.maxIter(); i++) {
    //termination
    if(std::sqrt(Gi.squaredNorm())<_solLQP.tolGFinal()) {
      if(_solLQP.callback()) {
        INFOV("Manifold-Projection succeed (it=%ld)!",i)
      }
      break;
    }
    V=Gi.squaredNorm();
    SolveNewton<T>::template solveLU<Vec>(DGDTheta,Gi,dir,true);
    //line search
    while(alpha>=tolLS) {
      qNext=q+dir*alpha;
      _info.reset(_body,qNext);
      Gi=G(dt,tau0,qdq,mapM(DGDTheta),w,mapM((MatT*)NULL));
      if(Gi.squaredNorm()<V) {
        q=qNext;
        break;
      } else alpha/=inc;
    }
    if(alpha<tolLS) {
      if(_solLQP.callback()) {
        INFOV("Manifold-Projection failed (iter=%d,dt=%f,GNorm=%f)!",
              i+1,std::to_double(dt),std::to_double(std::sqrt(Gi.squaredNorm())))
      }
      alpha=alphaInit?*alphaInit:1;
      q.resize(0);
      break;
    }
    alpha=std::min<T>(alpha*inc,1);
  }
  if(alphaInit)
    *alphaInit=alpha;
  return q;
}
//instance
template class MDPSimulator<double>;
#ifdef ALL_TYPES
template class MDPSimulator<__float128>;
template class MDPSimulator<mpfr::mpreal>;
#endif
PRJ_END
#endif
