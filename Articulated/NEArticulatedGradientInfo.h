#ifndef NE_ARTICULATED_GRADIENT_INFO_H
#define NE_ARTICULATED_GRADIENT_INFO_H

#include "ArticulatedBody.h"
#include "LeviCivita.h"
#include <Utils/SparseUtils.h>
#include <functional>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct NEArticulatedGradientInfoMap
{
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  NEArticulatedGradientInfoMap();
  NEArticulatedGradientInfoMap(const NEArticulatedGradientInfoMap& other);
  void reset(const ArticulatedBody& body,VecCM qMap,VecCM dqMap,VecCM ddqMap);
  void reset(const ArticulatedBody& body,VecCM qMap,VecCM dqMap);
  void reset(const ArticulatedBody& body,VecCM qMap);
  //-------------------------------------------------------------RNEA
  void RNEA(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,bool useDDQ=true);
  //-------------------------------------------------------------RNEADerivatives
  void DvDdq(const ArticulatedBody& body,sizeType i);
  void DvDq(const ArticulatedBody& body,sizeType i);
  void DaDdq(const ArticulatedBody& body,sizeType i,Vec6TCM vJ);
  void DaDq(const ArticulatedBody& body,sizeType i,Vec6TCM vJ);
  void Dfx(const ArticulatedBody& body,sizeType i,Vec6TCM fx);
  void Df(const ArticulatedBody& body,sizeType i,Vec6TCM vJ,Mat6XTCM IMCustom);
  void Dtau(const ArticulatedBody& body,sizeType i,sizeType j);
  void DHDq(const ArticulatedBody& body,sizeType i,sizeType j,MatTM dHdq) const;
  void DfP(const ArticulatedBody& body,sizeType i,sizeType j,bool dq);
  void RNEADerivatives(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,bool useDDQ=true);
  void RNEADerivativesFD(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecM q1,VecM dq1,VecM ddq1,T delta,bool useDDQ=true);
  MatTCM getDtauDddq() const;
  //-------------------------------------------------------------MDP
  void Sc(const ArticulatedBody& body,ColiCM joints,MatTM sc) const;
  void DScDqT(const ArticulatedBody& body,ColiCM joints,MatTM DscDqTr,VecCM r);
  void DScDq(const ArticulatedBody& body,ColiCM joints,MatTM DscDqr,VecCM r) const;
  void ScInvH(const ArticulatedBody& body,ColiCM joints,MatTM sc,MatTM scInvH) const;
  void ScInvHScT(const ArticulatedBody& body,ColiCM joints,MatTM sc,MatTM scInvH,MatTM scInvHScT) const;
  void DTG(sizeType k,const ArticulatedBody& body,Mat3X4T GK,std::function<void(sizeType,T)> dtg) const;    //derivative with respect to T
  void DTG(sizeType k,const ArticulatedBody& body,Mat3X4T GK,VecM dtg) const;                               //derivative with respect to T
  void DHDq(const ArticulatedBody& body,MatTM dHdq,VecCM q,Mat6XTCM IMCustom=mapCM((const Mat6XT*)NULL));
  Mat3X4T getTrans(sizeType JID) const;
  //-------------------------------------------------------------CRBA
  void CRBA(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0,bool useLTDL=true);
  void calcHInvH(const ArticulatedBody& body,bool useLTDL=true);
  void calcHInner(const ArticulatedBody& body,MatTM H,Mat6XTCM IMCustom=mapCM((const Mat6XT*)NULL));
  void calcH(const ArticulatedBody& body);
  void LTDL();
  template <typename M>
  void LTDLSolve(M x);
  void LTL();
  template <typename M>
  void LTLSolve(M x);
  //-------------------------------------------------------------ABA
  void ABA(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0);
  template <typename M,typename IM>
  void invert(M in,IM out);
  //-------------------------------------------------------------ABADerivatives
  void ABADerivatives(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0,bool useLTDL);
  void ABADerivativesFD(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0,VecM q1,VecM dq1,VecM ddq1,T delta,bool useLTDL);
  MatTCM getDddqDq() const;
  MatTCM getDddqDdq() const;
  MatTCM getDddqDtau0() const;
  MatTCM getH(const ArticulatedBody& body);
  MatTCM getInvH(const ArticulatedBody& body);
  //data
  ColiM _sparsityM,_nonZeroM;
  VecM _qM,_dqM,_ddqM,_tauM;
  Mat6XTM _invTPM,_invTM,_SM,_UM,_vM,_vPM,_aM,_fM,_ICM,_JTransM,_IM,_DvJDqM,_DdvJDqM,_DdvJDdqM;
  MatTM _DtauDqM,_DtauDdqM,_DvDqM,_DvDdqM,_DaDqM,_DaDdqM,_DfDqM,_DfDdqM;
protected:
  MatTM _HM,_invHM;
  bool _updateH,_updateInvH;
};
template <typename T>
struct NEArticulatedGradientInfo : public NEArticulatedGradientInfoMap<T>
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  using NEArticulatedGradientInfoMap<T>::_sparsityM;
  using NEArticulatedGradientInfoMap<T>::_nonZeroM;
  using NEArticulatedGradientInfoMap<T>::_qM;
  using NEArticulatedGradientInfoMap<T>::_dqM;
  using NEArticulatedGradientInfoMap<T>::_ddqM;
  using NEArticulatedGradientInfoMap<T>::_tauM;

  using NEArticulatedGradientInfoMap<T>::_invTPM;
  using NEArticulatedGradientInfoMap<T>::_invTM;
  using NEArticulatedGradientInfoMap<T>::_SM;
  using NEArticulatedGradientInfoMap<T>::_UM;
  using NEArticulatedGradientInfoMap<T>::_vM;
  using NEArticulatedGradientInfoMap<T>::_vPM;
  using NEArticulatedGradientInfoMap<T>::_aM;
  using NEArticulatedGradientInfoMap<T>::_fM;
  using NEArticulatedGradientInfoMap<T>::_ICM;
  using NEArticulatedGradientInfoMap<T>::_JTransM;
  using NEArticulatedGradientInfoMap<T>::_IM;
  using NEArticulatedGradientInfoMap<T>::_DvJDqM;
  using NEArticulatedGradientInfoMap<T>::_DdvJDqM;
  using NEArticulatedGradientInfoMap<T>::_DdvJDdqM;

  using NEArticulatedGradientInfoMap<T>::_DtauDqM;
  using NEArticulatedGradientInfoMap<T>::_DtauDdqM;
  using NEArticulatedGradientInfoMap<T>::_DvDqM;
  using NEArticulatedGradientInfoMap<T>::_DvDdqM;
  using NEArticulatedGradientInfoMap<T>::_DaDqM;
  using NEArticulatedGradientInfoMap<T>::_DaDdqM;
  using NEArticulatedGradientInfoMap<T>::_DfDqM;
  using NEArticulatedGradientInfoMap<T>::_DfDdqM;
  using NEArticulatedGradientInfoMap<T>::_HM;
  using NEArticulatedGradientInfoMap<T>::_invHM;
  NEArticulatedGradientInfo();
  NEArticulatedGradientInfo(const NEArticulatedGradientInfo& other);
  NEArticulatedGradientInfo(const ArticulatedBody& body,const Vec& q,const Vec& dq,const Vec& ddq);
  NEArticulatedGradientInfo(const ArticulatedBody& body,const Vec& q,const Vec& dq);
  NEArticulatedGradientInfo(const ArticulatedBody& body,const Vec& q);
  NEArticulatedGradientInfo& operator=(const NEArticulatedGradientInfo& other);
  void reset(const ArticulatedBody& body);
  void reset(const ArticulatedBody& body,const Vec& q,const Vec& dq,const Vec& ddq);
  void reset(const ArticulatedBody& body,const Vec& q,const Vec& dq);
  void reset(const ArticulatedBody& body,const Vec& q);
  void reorthogonalize(const ArticulatedBody& body);
  Mat3XT getTrans() const;
  //debug
  static void debug(const ArticulatedBody& body);
private:
  void resetPtr();
  //data
  Coli _sparsity,_nonZero;
  Vec _q,_dq,_ddq,_tau;
  Mat6XT _invTP,_invT,_S,_U,_v,_vP,_a,_f,_IC,_DvJDq,_DdvJDq,_DdvJDdq,_JTrans,_I;
  MatT _DtauDq,_DtauDdq,_DvDq,_DvDdq,_DaDq,_DaDdq,_DfDq,_DfDdq,_H,_invH;
};

PRJ_END

#endif
