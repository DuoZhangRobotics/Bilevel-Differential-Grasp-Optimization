#ifdef ENVIRONMENT_SUPPORT
#ifndef MDP_H
#define MDP_H

#include "MultiPrecisionLQP.h"
#include <Environment/EnvWrenchConstructor.h>
#include "NEArticulatedGradientInfo.h"
#include <Utils/DebugGradient.h>
#include <Utils/Options.h>

PRJ_BEGIN

//Solve the problem of following form:
//argmin c*w+w^THw/2
//  s.t. w_i>=0
//       \sum_i w_i<=1
//
//We solve this by minimizing:
//argmin f(w)=c*w+w^THw/2-\mu\sum_i\log(w_i)-\mu\log(1-\sum_i w_i)
//
//We then take Newton's step:
//H(w)\Delta(w)+G(w)=0
//We can simply use a line-search so that w+\Delta(w) does not violate the constraint
template <typename T>
class DebugWrenchConstructor : public EnvWrenchConstructor<T>
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  DebugWrenchConstructor(const std::vector<ExternalWrench<T>>& ref);
  void operator()(std::vector<ExternalWrench<T>>& externalWrench,const NEArticulatedGradientInfo<T>& info,bool grad=true) override;
  NEArticulatedGradientInfo<T> operator()(std::vector<ExternalWrench<T>>& externalWrench,const ArticulatedBody& body,VecCM qM);
  std::vector<ExternalWrench<T>> _refWrench;
};
template <typename T>
class MDP
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  MDP(const ArticulatedBody& body,Options& ops);
  MDP(const MDP<T>& other);
  bool solveMDPQPF(MatTM DdqHatDq,MatTM DdqHatDdq,bool profileMDPError=false);
  bool solveMDPQPI(MatTM DdqHatDq,MatTM DdqHatDdq,bool profileMDPError=false);
  bool solveMDPLPF(MatTM DddqDq,MatTM DddqDdq,bool profileMDPError=false);
  bool solveMDPLPI(MatTM DddqDq,MatTM DddqDdq,bool profileMDPError=false);
  MDP<T>& operator=(const MDP<T>& other);
  virtual void reset(Options& ops);
  void randomize(T scale,bool QP);
  void readAndTestProb();
  void assembleHcQP();
  void assembleHcLP();
  void assembleHcLPI();
  sizeType nW() const;
  MatT B() const;
  void Bf(const Vec& f);
  void BT(const MatT& b,MatT& BTb) const;
  void DBDq(const Vec& b,MatT& DBDqb) const;
  void DBDqT(const Vec& b,MatT& DBDqTb) const;
  void DXDq(MatT& DBDq,sizeType r,sizeType JID,const Mat3X4T& DBDX) const;
  T kineticEnergyInfo(const Vec& w);
  T kineticEnergyAssembled(const Vec& w) const;
  static void debugWrenchConstructor(const ArticulatedBody& body,EnvWrenchConstructor<T>& debugWrench,sizeType nrIter);
  static void debug(const ArticulatedBody& body,sizeType nrIter,T scale=1000);
  //in/out
  Vec _deltaDq,_dqMDP,_lastW,_Bf;
  MultiPrecisionLQP<T> _solLQP;
  const ArticulatedBody& _body;
  NEArticulatedGradientInfo<T> _info;
  std::vector<ExternalWrench<T>> _externalWrench;
protected:
  std::vector<sizeType> _joints;
  MatT _Sc,_ScInvH,_ScInvHScT;
  MatT _BTSc,_BTScInvH,_Sigma,_invHLQP;
  MatT _DBDqT,_DScDq,_BTDScDq;
  MatT _DHDq,_DScDqT,_DBDq;
};

PRJ_END

#endif
#endif
