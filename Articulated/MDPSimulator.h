#ifdef ENVIRONMENT_SUPPORT
#ifndef MDP_SIMULATOR_H
#define MDP_SIMULATOR_H

#include "MDP.h"
#include <Environment/EnvWrenchConstructor.h>

PRJ_BEGIN

enum MDP_SIMULATOR_MODE {
  INVERSE_I,
  FORWARD_I,
  INVERSE_LF,
  FORWARD_LF,
  INVERSE_BACKWARD_RK1F,
  BACKWARD_RK1F,
  FORWARD_RK1F,
  FORWARD_RK1I,
  FORWARD_RK2I,
  FORWARD_RK4I,
  NMDP_PGM,
  NMDP_ZOPGM,
};
template <typename T>
struct MDPSimulatorTraits;
struct PDTarget;
template <typename T>
class MDPSimulator : public MDP<T>
{
public:
  DECL_MAP_TYPES_T
  using MDP<T>::mapV;
  using MDP<T>::mapCV;
  using MDP<T>::mapV2CV;
  using MDP<T>::mapM;
  using MDP<T>::mapCM;
  using MDP<T>::reset;
  using MDP<T>::_deltaDq;
  using MDP<T>::_dqMDP;
  using MDP<T>::_lastW;
  using MDP<T>::_solLQP;
  using MDP<T>::_body;
  using MDP<T>::_info;
  using MDP<T>::_externalWrench;
  using MDP<T>::_joints;
  using MDP<T>::_Sc;
  friend struct MDPSimulatorTraits<T>;
  typedef std::vector<Vec,Eigen::aligned_allocator<Vec>> Vecss;
  typedef std::vector<DMat,Eigen::aligned_allocator<DMat>> DMatss;
  typedef std::shared_ptr<EnvWrenchConstructor<T>> WrenchConstructor;
  MDPSimulator(const ArticulatedBody& body,Options& ops,const Vec3T& g,MDP_SIMULATOR_MODE integrator=FORWARD_RK1F);
  MDPSimulator(const MDPSimulator<T>& other);
  MDPSimulator<T>& operator=(const MDPSimulator<T>& other);
  virtual void reset(Options& ops) override;
  WrenchConstructor getWrenchConstructor() const;
  void setWrenchConstructor(WrenchConstructor wrench);
  void setPDTarget(std::shared_ptr<PDTarget> PDTarget=NULL);
  void writeVTKSeq(const std::string& path,T horizon,T interval,T dt,std::shared_ptr<PDTarget> PDTarget=NULL,const std::map<std::string,std::set<sizeType>>* jointMask=NULL);
  void writeVTK(const Vec& q,const std::string& path,const std::map<std::string,std::set<sizeType>>* jointMask=NULL) const;
  void setState(const Vec& q,const Vec& dq);
  Vec step(T dt,const Vec* tau0,const Vec& qdq,DMat* Dqdq,DMat* Dtau,T* betaInit=NULL);
  Vec step(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau,T* betaInit=NULL);
  Vecss stepBatched(T dt,const Vecss* tau0,const Vecss& qdq,DMatss* Dqdq,DMatss* Dtau);
  void debug(T dt,sizeType nrIter,T scale=1000);
  void debugNMDP(T dt,sizeType nrIter,T scale=1000);
  void setMode(MDP_SIMULATOR_MODE m);
  MDP_SIMULATOR_MODE getMode() const;
  void setWarmStart(bool warm);
  void begLog(const std::string& path);
  void endLog();
  Vec3T g() const;
  void log(T tw,T ts,T tne);
private:
  void calcIMCustom(const std::vector<EndEffectorBounds>& EE);
  //forward/backward impulse
  Vec inverseI(VecCM qdqddq,MatTM jac);
  Vec forwardI(VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau);
  //forward/backward limiting force
  Vec inverseLF(VecCM qdqddq,MatTM jac);
  Vec forwardLF(VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau);
  //forward/backward force+RK1
  Vec inverseBackwardRK1F(T dt,VecCM tau0,VecCM qdqNext_qdq,MatTM DqdqNext_qdq,MatTM Dtau);
  Vec backwardRK1F(T dt,VecCM tau0,VecCM qdqNext_qdq,MatTM DqdqNext_qdq,MatTM Dtau);
  Vec forwardRK1F(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau);
  //forward impulse+RK1/2/4
  Vec forwardRK1I(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau);
  Vec forwardRK2I(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau);
  Vec forwardRK4I(T dt,VecCM tau0,VecCM qdq,MatTM Dqdq,MatTM Dtau);
  Vec getCtrl(VecCM tau0,VecCM qdq) const;
  //NMDP
  T K(T dt,VecCM qdq,VecM DKDTheta,MatTM DDKDDTheta);
  Vec G(T dt,VecCM tau0,VecCM qdq,MatTM DGDTheta,Vec* w,MatTM DGDw);
  Vec stepNMDPPGM(T dt,VecCM tau0,VecCM qdq,T* betaInit=NULL);
  bool stepNMDPPGMAdaptive(T dt,VecCM tau0,VecM qdq,T* betaInit=NULL);
  Vec manifoldProjection(T dt,VecCM tau0,VecCM qdq,VecCM init,Vec* w,T* alphaInit=NULL);
  //data
  std::shared_ptr<std::ofstream> _log;
  MatT _DdqHatDq,_DdqHatDdq;
  MatT _Dqdq2,_Dqdq3,_Dqdq4;
  MatT _Dtau2,_Dtau3,_Dtau4;
  Mat6XT _IMCustom;
  WrenchConstructor _wrench;
  std::shared_ptr<PDTarget> _PDTarget;
  MDP_SIMULATOR_MODE _mode;
  bool _warmStart;
  Vec6T _a0;
  T _tolK,_minDt;
  T _betaMin,_alphaBetaInc;
  bool _useKEE;
  //transient data
  sizeType _iter,_iterG;
};

PRJ_END

#endif
#endif
