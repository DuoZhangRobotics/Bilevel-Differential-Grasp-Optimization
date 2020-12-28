#ifndef PBD_ARTICULATED_GRADIENT_INFO_H
#define PBD_ARTICULATED_GRADIENT_INFO_H

#include "ArticulatedBody.h"
#include "LeviCivita.h"
#include <functional>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

#include <Utils/ArticulatedBodyPragma.h>
template <typename T>
struct PBDArticulatedGradientInfoMap
{
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  PBDArticulatedGradientInfoMap();
  PBDArticulatedGradientInfoMap(const PBDArticulatedGradientInfoMap& other);
  void resetLambda(const ArticulatedBody& body,VecCM lambdaMap);
  void reset(const ArticulatedBody& body,VecCM xMap);
  Mat3X4T DTDLambda(sizeType k) const;
  //-------------------------------------------------------------toolTG
  T TG(const ArticulatedBody& body,Mat3XTCM G) const;
  //-------------------------------------------------------------toolDTG: G is modified
  void DTG(sizeType k,const ArticulatedBody& body,Mat3X4T GK,std::function<void(sizeType,T)> DTG) const;
  void DTG(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const;
  void DTGZ(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const;
  void DTGBF(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const;
  void DTGBFZ(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const;
  //-------------------------------------------------------------toolA: this is Ti, MRR,MRt,MtR,Mtt are modified
  //this is an interface whose input is 4 3X3 tensor
  void toolA(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3T MRR,Mat3T MRt,Mat3T MtR,Mat3T Mtt,std::function<void(sizeType,sizeType,T)> A) const;
  void toolA(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,std::function<void(sizeType,sizeType,T)> A) const;
  void toolANonRecursivePhase1Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1KiMwKj,const Mat3T& wK_1KiMtKj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolANonRecursivePhase1Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& tK_1KiMwKj,const Mat3T& tK_1KiMtKj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolANonRecursivePhase2Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMwK_1Kj,const Mat3T& tK_1iMwK_1Kj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolANonRecursivePhase2Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMtK_1Kj,const Mat3T& tK_1iMtK_1Kj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolARecursive(sizeType k,sizeType p,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt) const;
  void toolARecursiveInplace(sizeType k,sizeType p,const PBDArticulatedGradientInfoMap& Tj,Mat3T& MRR,Mat3T& MRt,Mat3T& MtR,Mat3T& Mtt) const;
  void toolAZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTM A) const;
  //this is an interface whose input is a 12X12 tensor
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  struct ToolAContractTensor
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,T val)
    {
      //MRR
      sizeType t=k*4;
      sizeType k3=k*3;
      MRR(LeviCivitaPositive<alpha,0>::value,k3+LeviCivitaPositive<gamma,0>::value)+=val*Ti._TM(LeviCivitaPositive<alpha,1>::value,t+theta)*Tj._TM(LeviCivitaPositive<gamma,1>::value,t+mu);
      MRR(LeviCivitaPositive<alpha,0>::value,k3+LeviCivitaNegative<gamma,0>::value)-=val*Ti._TM(LeviCivitaPositive<alpha,1>::value,t+theta)*Tj._TM(LeviCivitaNegative<gamma,1>::value,t+mu);
      MRR(LeviCivitaNegative<alpha,0>::value,k3+LeviCivitaPositive<gamma,0>::value)-=val*Ti._TM(LeviCivitaNegative<alpha,1>::value,t+theta)*Tj._TM(LeviCivitaPositive<gamma,1>::value,t+mu);
      MRR(LeviCivitaNegative<alpha,0>::value,k3+LeviCivitaNegative<gamma,0>::value)+=val*Ti._TM(LeviCivitaNegative<alpha,1>::value,t+theta)*Tj._TM(LeviCivitaNegative<gamma,1>::value,t+mu);
    }
  };
  template <sizeType alpha,sizeType gamma,sizeType mu>
  struct ToolAContractTensor<alpha,3,gamma,mu>
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,T val)
    {
      //MtR
      sizeType t=k*4;
      sizeType k3=k*3;
      MtR(alpha,k3+LeviCivitaPositive<gamma,0>::value)+=val*Tj._TM(LeviCivitaPositive<gamma,1>::value,t+mu);
      MtR(alpha,k3+LeviCivitaNegative<gamma,0>::value)-=val*Tj._TM(LeviCivitaNegative<gamma,1>::value,t+mu);
    }
  };
  template <sizeType alpha,sizeType theta,sizeType gamma>
  struct ToolAContractTensor<alpha,theta,gamma,3>
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,T val)
    {
      //MRt
      sizeType t=k*4;
      sizeType k3=k*3;
      MRt(LeviCivitaPositive<alpha,0>::value,k3+gamma)+=val*Ti._TM(LeviCivitaPositive<alpha,1>::value,t+theta);
      MRt(LeviCivitaNegative<alpha,0>::value,k3+gamma)-=val*Ti._TM(LeviCivitaNegative<alpha,1>::value,t+theta);
    }
  };
  template <sizeType alpha,sizeType gamma>
  struct ToolAContractTensor<alpha,3,gamma,3>
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,T val)
    {
      //Mtt
      sizeType k3=k*3;
      Mtt(alpha,k3+gamma)+=val;
    }
  };
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  void toolAContractTensor(sizeType k,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,T val) const
  {
    ToolAContractTensor<alpha,theta,gamma,mu>::contract(k,*this,Tj,MRR,MRt,MtR,Mtt,val);
  }
  void toolAContactAll(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTCM M) const;
  void toolA(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const;
  void toolAZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const;
  //this is for internal force's separate term
  void toolALR(sizeType kL,sizeType kR,const ArticulatedBody& body,const Mat3T& tensor,const Vec3T& ptL,const Vec3T& ptR,std::function<void(sizeType,sizeType,T)> A) const;
  void toolALR(sizeType kL,sizeType kR,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM tensor,MatTM A) const;
  //this is a brute-force interface for testing
  void toolABF(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const;
  void toolABFZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const;
  void toolABF2(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTM A) const;
  void toolABF2Z(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTM A) const;
  void toolABF2Z(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const;
  //-------------------------------------------------------------toolB: G is modified
  void toolB(const ArticulatedBody& body,Mat3XTM G,std::function<void(sizeType,sizeType,T)> B) const;
  void toolBNonRecursivePhase1Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1KiMwKj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolBNonRecursivePhase1Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& tK_1KiMwKj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolBNonRecursivePhase2Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMwK_1Kj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolBNonRecursivePhase2Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMtK_1Kj,std::function<void(sizeType,sizeType,T)> A) const;
  void toolBRecursive(sizeType k,sizeType p,Mat3XTM G) const;
  void toolBRecursiveInplace(sizeType k,sizeType p,Mat3X4T& G) const;
  void toolBZ(const ArticulatedBody& body,Mat3XTM G,MatTM B) const;
  //this is an interface whose input is a 12X12 tensor
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  struct ToolBContractTensor
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,Mat3XTM G,T valC)
    {
      sizeType k2=k<<1,k4=k*4;
      G(gamma,mu+k4)+=(Ti._DTLambdaM(LeviCivitaPositive<alpha,0>::value,k2+0)*Ti._TM(LeviCivitaPositive<alpha,1>::value,k4+theta)-
                       Ti._DTLambdaM(LeviCivitaNegative<alpha,0>::value,k2+0)*Ti._TM(LeviCivitaNegative<alpha,1>::value,k4+theta))*valC;
    }
  };
  template <sizeType alpha,sizeType gamma,sizeType mu>
  struct ToolBContractTensor<alpha,3,gamma,mu>
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,Mat3XTM G,T valC)
    {
      sizeType k2=k<<1,k4=k*4;
      G(gamma,mu+k4)+=Ti._DTLambdaM(alpha,k2+1)*valC;
    }
  };
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  void toolBContractTensor(sizeType k,Mat3XTM G,T valC) const
  {
    ToolBContractTensor<alpha,theta,gamma,mu>::contract(k,*this,G,valC);
  }
  void toolBContactAll(const ArticulatedBody& body,MatTCM M,Mat3XTM G) const;
  //this is a brute-force interface for testing
  void toolBBF(const ArticulatedBody& body,MatTCM M,MatTM B) const;
  void toolBBFZ(const ArticulatedBody& body,MatTCM M,MatTM B) const;
  //-------------------------------------------------------------toolA,toolB combined
  //this is an interface whose input is 4 3X3 tensor
  void toolAB(sizeType k,const ArticulatedBody& body,Mat3T MRR,Mat3T MRt,Mat3T MtR,Mat3T Mtt,Mat3X4T G,std::function<void(sizeType row,sizeType col,T val)> AB) const;
  void toolAB(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,std::function<void(sizeType row,sizeType col,T val)> AB) const;
  void toolABZ(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,MatTM AB) const;
  //this is an interface whose input is a 12X12 tensor
  void toolAB(const ArticulatedBody& body,MatTCM M,Mat3XTM G,MatTM AB) const;
  void toolABZ(const ArticulatedBody& body,MatTCM M,Mat3XTM G,MatTM AB) const;
  void toolABBF(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,MatTM AB) const;
  void toolABBFZ(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,MatTM AB) const;
  void toolABBFZ(const ArticulatedBody& body,MatTCM M,Mat3XTM G,MatTM AB) const;
  //-------------------------------------------------------------toolA,toolC combined
  //this is an interface whose input is 4 3X3 tensor
  void toolAC(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,std::function<void(sizeType row,sizeType col,T val)> AC) const;
  void toolACNonRecursivePhase1Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1KiMwKj,const Mat3T& wK_1KiMtKj,const Mat3T& wK_1KiLambdaMwKj,const Mat3T& wK_1KiLambdaMtKj,std::function<void(sizeType row,sizeType col,T val)> AC) const;
  void toolACNonRecursivePhase2Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iLambdaMwK_1Kj,const Mat3T& tK_1iLambdaMwK_1Kj,std::function<void(sizeType row,sizeType col,T val)> AC) const;
  void toolACNonRecursivePhase2Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iLambdaMtK_1Kj,const Mat3T& tK_1iLambdaMtK_1Kj,std::function<void(sizeType row,sizeType col,T val)> AC) const;
  void toolACRecursive(sizeType k,sizeType p,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC) const;
  void toolACZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTM AC) const;
  //this is an interface whose input is a 12X12 tensor
  template <sizeType r>
  T CWLambdaR(sizeType k,sizeType c) const
  {
    return _DTLambdaM(LeviCivitaPositive<r,0>::value,k<<1)*_TM(LeviCivitaPositive<r,1>::value,c)-_DTLambdaM(LeviCivitaNegative<r,0>::value,k<<1)*_TM(LeviCivitaNegative<r,1>::value,c);
  }
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  struct ToolCContractTensor
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,T valC)
    {
      //MRR
      sizeType t=k*4;
      sizeType k3=k*3;
      MRRA(LeviCivitaPositive<alpha,0>::value,k3+LeviCivitaPositive<gamma,0>::value)+=valC*Ti.CWLambdaR<LeviCivitaPositive<alpha,1>::value>(k,t+theta)*Tj._TM(LeviCivitaPositive<gamma,1>::value,t+mu);
      MRRA(LeviCivitaPositive<alpha,0>::value,k3+LeviCivitaNegative<gamma,0>::value)-=valC*Ti.CWLambdaR<LeviCivitaPositive<alpha,1>::value>(k,t+theta)*Tj._TM(LeviCivitaNegative<gamma,1>::value,t+mu);
      MRRA(LeviCivitaNegative<alpha,0>::value,k3+LeviCivitaPositive<gamma,0>::value)-=valC*Ti.CWLambdaR<LeviCivitaNegative<alpha,1>::value>(k,t+theta)*Tj._TM(LeviCivitaPositive<gamma,1>::value,t+mu);
      MRRA(LeviCivitaNegative<alpha,0>::value,k3+LeviCivitaNegative<gamma,0>::value)+=valC*Ti.CWLambdaR<LeviCivitaNegative<alpha,1>::value>(k,t+theta)*Tj._TM(LeviCivitaNegative<gamma,1>::value,t+mu);
    }
  };
  template <sizeType alpha,sizeType gamma,sizeType mu>
  struct ToolCContractTensor<alpha,3,gamma,mu>
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,T valC) {}
  };
  template <sizeType alpha,sizeType theta,sizeType gamma>
  struct ToolCContractTensor<alpha,theta,gamma,3>
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,T valC)
    {
      //MRt
      sizeType t=k*4;
      sizeType k3=k*3;
      MRtA(LeviCivitaPositive<alpha,0>::value,k3+gamma)+=valC*Ti.CWLambdaR<LeviCivitaPositive<alpha,1>::value>(k,t+theta);
      MRtA(LeviCivitaNegative<alpha,0>::value,k3+gamma)-=valC*Ti.CWLambdaR<LeviCivitaNegative<alpha,1>::value>(k,t+theta);
    }
  };
  template <sizeType alpha,sizeType gamma>
  struct ToolCContractTensor<alpha,3,gamma,3>
  {
    FORCE_INLINE static void contract(sizeType k,const PBDArticulatedGradientInfoMap& Ti,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,T valC) {}
  };
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  void toolCContractTensor(sizeType k,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,T valC) const
  {
    toolAContractTensor<alpha,theta,gamma,mu>(k,Tj,MRRC,MRtC,MtRC,MttC,valC);
    ToolCContractTensor<alpha,theta,gamma,mu>::contract(k,*this,Tj,MRRA,MRtA,MtRA,MttA,valC);
  }
  void toolCContactAll(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC) const;
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  void toolACContractTensor(sizeType k,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,T valA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,T valC) const
  {
    toolAContractTensor<alpha,theta,gamma,mu>(k,Tj,MRRA,MRtA,MtRA,MttA,valA);
    toolCContractTensor<alpha,theta,gamma,mu>(k,Tj,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,valC);
  }
  void toolACContactAll(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,MatTCM MA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC) const;
  void toolAC(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM MA,MatTCM MC,MatTM AC) const;
  void toolACZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM MA,MatTCM MC,MatTM AC) const;
  //this is for internal force's separate term
  void toolCLR(sizeType kL,sizeType kR,const ArticulatedBody& body,const Mat3T& tensor,const Vec3T& ptL,const Vec3T& ptR,std::function<void(sizeType row,sizeType col,T val)> C) const;
  void toolCLR(sizeType kL,sizeType kR,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM tensor,MatTM A) const;
  //this is a brute-force interface for testing
  void toolCBF(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM C) const;
  void toolCBFZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM C) const;
  //-------------------------------------------------------------toolD using brute-force method for testing
  void toolDBF(const ArticulatedBody& body,Mat3XTM G,MatTM D) const;
  void toolDRecursive(sizeType k,sizeType p,Mat3XTM G,Mat3XTM GB,T coef) const;
  void toolDBFZ(const ArticulatedBody& body,Mat3XTM G,MatTM D) const;
  //-------------------------------------------------------------toolA,B,C,CT,D combined
  //this is an interface whose input is 4 3X3 tensor
  void toolABCCTD(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,Mat3XTM G,Mat3XTM GB,std::function<void(sizeType row,sizeType col,T val)> ABCCTD) const;
  void toolABCCTDZ(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,Mat3XTM G,Mat3XTM GB,MatTM ABCCTD) const;
  //this is an interface whose input is a 12X12 tensor
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  void toolBCCTDContractTensor(sizeType k,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,T valC,Mat3XTM GB) const
  {
    toolCContractTensor<alpha,theta,gamma,mu>(k,*this,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,valC);
    toolBContractTensor<alpha,theta,gamma,mu>(k,GB,valC);
  }
  void toolBCCTDContactAll(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC,Mat3XTM GB) const;
  template <sizeType alpha,sizeType theta,sizeType gamma,sizeType mu>
  void toolABCCTDContractTensor(sizeType k,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,T valA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,T valC,Mat3XTM GB) const
  {
    toolBCCTDContractTensor<alpha,theta,gamma,mu>(k,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,valC,GB);
    toolAContractTensor<alpha,theta,gamma,mu>(k,*this,MRRA,MRtA,MtRA,MttA,valA/2);
  }
  void toolABCCTDContactAll(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,MatTCM MA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC,Mat3XTM GB) const;
  void toolABCCTD(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const;
  void toolABCCTDZ(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const;
  //this is a brute-force interface for testing
  void toolABCCTDBF(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const;
  void toolABCCTDBFZ(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const;
  //-------------------------------------------------------------Jacobian RC
  void JRSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR) const;
  void JRMSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Mat3T&)> JR) const;
  void JCSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JC) const;
  void JRCSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR,std::function<void(sizeType,const Vec3T&)> JC) const;
  void JRSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Mat3T& dR) const;
  void JCSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Vec3T& dC) const;
  Mat3T JRSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const;
  Vec3T JCSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const;
  void JRSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dR) const;
  void JCSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dC) const;
  //-------------------------------------------------------------Jacobian RCLambda
  void JRILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR) const;
  void JRMILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Mat3T&)> JR) const;
  void JCILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JC) const;
  void JRCILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR,std::function<void(sizeType,const Vec3T&)> JC) const;
  void JRILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Mat3T& dR) const;
  void JCILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Vec3T& dC) const;
  Mat3T JRILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const;
  Vec3T JCILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const;
  void JRILambdaSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dR) const;
  void JCILambdaSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dC) const;
  //-------------------------------------------------------------Hessian RK_1K
  void HRK_1KSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,sizeType,const Mat3T&)> HR) const;
  //-------------------------------------------------------------Hessian RK_1KLambda
  void HRK_1KLambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,sizeType,const Mat3T&)> HR) const;
  //data
  VecM _xM;
  Mat3XTM _TM,_TK_1KM,_DTM,_RDTM,_DDTM,_JTransM;
  Mat3XTM _DTLambdaM,_DTK_1KLambdaM,_DTILambdaM,_RDTILambdaM,_DTIILambdaM,_RDTIILambdaM;
};
template <typename T>
struct PBDArticulatedGradientInfo : public PBDArticulatedGradientInfoMap<T>
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  using PBDArticulatedGradientInfoMap<T>::_xM;
  using PBDArticulatedGradientInfoMap<T>::_TM;
  using PBDArticulatedGradientInfoMap<T>::_TK_1KM;
  using PBDArticulatedGradientInfoMap<T>::_DTM;
  using PBDArticulatedGradientInfoMap<T>::_RDTM;
  using PBDArticulatedGradientInfoMap<T>::_DDTM;
  using PBDArticulatedGradientInfoMap<T>::_JTransM;
  using PBDArticulatedGradientInfoMap<T>::_DTLambdaM;
  using PBDArticulatedGradientInfoMap<T>::_DTK_1KLambdaM;
  using PBDArticulatedGradientInfoMap<T>::_DTILambdaM;
  using PBDArticulatedGradientInfoMap<T>::_RDTILambdaM;
  using PBDArticulatedGradientInfoMap<T>::_DTIILambdaM;
  using PBDArticulatedGradientInfoMap<T>::_RDTIILambdaM;
  PBDArticulatedGradientInfo();
  PBDArticulatedGradientInfo(const PBDArticulatedGradientInfo& other);
  PBDArticulatedGradientInfo(const ArticulatedBody& body,const Vec& x);
  PBDArticulatedGradientInfo& operator=(const PBDArticulatedGradientInfo& other);
  void resetLambda(const ArticulatedBody& body,const Vec& lambda);
  void reset(const ArticulatedBody& body,const Vec& x);
  void reorthogonalize(const ArticulatedBody& body);
  //debug
  static void debug(const ArticulatedBody& body);
private:
  void resetPtr();
  //data
  Vec _x;
  Mat3XT _T,_TK_1K,_DT,_RDT,_DDT,_JTrans;
  Mat3XT _DTLambda,_DTK_1KLambda,_DTILambda,_RDTILambda,_DTIILambda,_RDTIILambda;
};

PRJ_END

#endif
