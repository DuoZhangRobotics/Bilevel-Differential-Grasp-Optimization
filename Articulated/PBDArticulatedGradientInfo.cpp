#include "PBDArticulatedGradientInfo.h"
#include "TensorContractPragma.h"
#include "JointFunc.h"
#include <Utils/CrossSpatialUtil.h>
#include <Utils/DebugGradient.h>
#include <Utils/Scalar.h>

USE_PRJ_NAMESPACE

//block
#define MRRI MRR.template block<3,3>(0,k*3)
#define MRtI MRt.template block<3,3>(0,k*3)
#define MtRI MtR.template block<3,3>(0,k*3)
#define MttI Mtt.template block<3,3>(0,k*3)
//blockA
#define MRRAI MRRA.template block<3,3>(0,k*3)
#define MRtAI MRtA.template block<3,3>(0,k*3)
#define MtRAI MtRA.template block<3,3>(0,k*3)
#define MttAI MttA.template block<3,3>(0,k*3)
//blockC
#define MRRCI MRRC.template block<3,3>(0,k*3)
#define MRtCI MRtC.template block<3,3>(0,k*3)
#define MtRCI MtRC.template block<3,3>(0,k*3)
#define MttCI MttC.template block<3,3>(0,k*3)

//CPU function API
template <typename T>
PBDArticulatedGradientInfoMap<T>::PBDArticulatedGradientInfoMap()
  :_xM(NULL,0),
   _TM(NULL,3,0,Eigen::OuterStride<>(3)),
   _TK_1KM(NULL,3,0,Eigen::OuterStride<>(3)),
   _DTM(NULL,3,0,Eigen::OuterStride<>(3)),
   _RDTM(NULL,3,0,Eigen::OuterStride<>(3)),
   _DDTM(NULL,3,0,Eigen::OuterStride<>(3)),
   _JTransM(NULL,3,0,Eigen::OuterStride<>(3)),
   _DTLambdaM(NULL,3,0,Eigen::OuterStride<>(3)),
   _DTK_1KLambdaM(NULL,3,0,Eigen::OuterStride<>(3)),
   _DTILambdaM(NULL,3,0,Eigen::OuterStride<>(3)),
   _RDTILambdaM(NULL,3,0,Eigen::OuterStride<>(3)),
   _DTIILambdaM(NULL,3,0,Eigen::OuterStride<>(3)),
   _RDTIILambdaM(NULL,3,0,Eigen::OuterStride<>(3)) {}
template <typename T>
PBDArticulatedGradientInfoMap<T>::PBDArticulatedGradientInfoMap(const PBDArticulatedGradientInfoMap& other)
  :_xM((VecM&)other._xM),
   _TM((Mat3XTM&)other._TM),
   _TK_1KM((Mat3XTM&)other._TK_1KM),
   _DTM((Mat3XTM&)other._DTM),
   _RDTM((Mat3XTM&)other._RDTM),
   _DDTM((Mat3XTM&)other._DDTM),
   _JTransM((Mat3XTM&)other._JTransM),
   _DTLambdaM((Mat3XTM&)other._DTLambdaM),
   _DTK_1KLambdaM((Mat3XTM&)other._DTK_1KLambdaM),
   _DTILambdaM((Mat3XTM&)other._DTILambdaM),
   _RDTILambdaM((Mat3XTM&)other._DTILambdaM),
   _DTIILambdaM((Mat3XTM&)other._DTIILambdaM),
   _RDTIILambdaM((Mat3XTM&)other._RDTIILambdaM) {}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::resetLambda(const ArticulatedBody& body,VecCM lambdaMap)
{
  Vec2i begEnd;
  Mat3X4T JTrans;
  sizeType nrJ=body.nrJ();
  for(sizeType i=0; i<nrJ; i++) {
    const Joint& J=body.joint(i);
    JTrans=TRANSI(_JTransM,i);
    //compute DTLLambda
    begEnd=J.RBegEnd();
    DWDLI(_DTK_1KLambdaM,i)=_DTM.block(0,begEnd[0],3,begEnd[1]-begEnd[0])*lambdaMap.segment(begEnd[0],begEnd[1]-begEnd[0]);
    begEnd=J.CBegEnd();
    DTDLI(_DTK_1KLambdaM,i)=_DTM.block(0,begEnd[0],3,begEnd[1]-begEnd[0])*lambdaMap.segment(begEnd[0],begEnd[1]-begEnd[0]);
    //compute DTILambda
    JointFunc<T>::DDTLambda(J,_DTILambdaM,mapM2CM(_DDTM),lambdaMap);
    //compute DTIILambda
    JointFunc<T>::DDDTLambda(J,mapV2CV(_xM),lambdaMap,_DTIILambdaM);
    //compute RDTIILambda
    Eigen::Block<Mat3XTM> DILambdaMBlk=_DTILambdaM.block(0,J._offDOF,3,J.nrDOF());
    Eigen::Block<Mat3XTM> RDILambdaMBlk=_RDTILambdaM.block(0,J._offDOF,3,J.nrDOF());
    Eigen::Block<Mat3XTM> DIILambdaMBlk=_DTIILambdaM.block(0,J._offDDT,3,J.nrDDT());
    Eigen::Block<Mat3XTM> RDIILambdaMBlk=_RDTIILambdaM.block(0,J._offDDT,3,J.nrDDT());
    DIILambdaMBlk=ROTI(_JTransM,i)*DIILambdaMBlk;
    if(J._parent >= 0) {
      DWDLI(_DTLambdaM,i)=DWDLI(_DTLambdaM,J._parent)+ROTI(_TM,J._parent)*DWDLI(_DTK_1KLambdaM,i);
      DTDLI(_DTLambdaM,i)=DWDLI(_DTLambdaM,J._parent).cross(CTRI(_TM,i)-CTRI(_TM,J._parent))+DTDLI(_DTLambdaM,J._parent)+ROTI(_TM,J._parent)*DTDLI(_DTK_1KLambdaM,i);
      RDILambdaMBlk=ROTI(_TM,J._parent)*DILambdaMBlk;
      RDIILambdaMBlk=ROTI(_TM,J._parent)*DIILambdaMBlk;
    } else {
      DWDLI(_DTLambdaM,i)=DWDLI(_DTK_1KLambdaM,i);
      DTDLI(_DTLambdaM,i)=DTDLI(_DTK_1KLambdaM,i);
      RDILambdaMBlk=DILambdaMBlk;
      RDIILambdaMBlk=DIILambdaMBlk;
    }
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::reset(const ArticulatedBody& body,VecCM xMap)
{
  Mat3X4T JTrans;
  sizeType nrJ=body.nrJ();
  for(sizeType i=0; i<nrJ; i++) {
    const Joint& J=body.joint(i);
    JTrans=TRANSI(_JTransM,i);
    //compute local transformation
    GETTM_T(JLT,_TK_1KM,i)
    JLT=JointFunc<T>::TDTDDT(J,xMap,_DTM,_DDTM);
    APPLY_TRANS(JLT,JTrans,JLT);
    Eigen::Block<Mat3XTM> DTMBlk=_DTM.block(0,J._offDOF,3,J.nrDOF());
    Eigen::Block<Mat3XTM> RDTMBlk=_RDTM.block(0,J._offDOF,3,J.nrDOF());
    Eigen::Block<Mat3XTM> DDTMBlk=_DDTM.block(0,J._offDDT,3,J.nrDDT());
    DTMBlk=ROT(JTrans)*DTMBlk;
    DDTMBlk=ROT(JTrans)*DDTMBlk;
    //compute parent transformation, RDT
    GETTM_T(JT,_TM,i)
    if(J._parent >= 0) {
      GETTM_T(PJT,_TM,J._parent)
      APPLY_TRANS(JT,PJT,JLT)
      RDTMBlk=ROT(PJT)*DTMBlk;
    } else {
      JT=JLT;
      RDTMBlk=DTMBlk;
    }
  }
}
template <typename T>
typename PBDArticulatedGradientInfoMap<T>::Mat3X4T PBDArticulatedGradientInfoMap<T>::DTDLambda(sizeType k) const
{
  return concatCol(cross<T>(DWDLI(_DTLambdaM,k))*ROTI(_TM,k),DTDLI(_DTLambdaM,k));
}
//-------------------------------------------------------------toolTG
template <typename T>
T PBDArticulatedGradientInfoMap<T>::TG(const ArticulatedBody& body,Mat3XTCM G) const
{
  T ret=0;
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    ret+=(ROTI(_TM,k)*ROTI(G,k).transpose()).trace();
    ret+=CTRI(_TM,k).dot(CTRI(G,k));
  }
  return ret;
}
//-------------------------------------------------------------toolDTG: G is modified
template <typename T>
void PBDArticulatedGradientInfoMap<T>::DTG(sizeType k,const ArticulatedBody& body,Mat3X4T GK,std::function<void(sizeType,T)> DTG) const
{
  Vec2i begEnd;
  Mat3T coefOmega;
  Vec3T coefOmegaV,coefT;
  for(; k>=0; k=body.joint(k)._parent) {
    const Joint& J=body.joint(k);
    GETTCM_T(JLT,_TK_1KM,k)
    //gradient R,C
    coefOmega=ROT(JLT)*ROT(GK).transpose();
    coefT=CTR(GK);
    if(J._parent>=0) {
      coefOmega=coefOmega*ROTI(_TM,J._parent);
      coefT=ROTI(_TM,J._parent).transpose()*coefT;
    }
    coefOmegaV=invCrossMatTrace<T>(coefOmega);
    //assemble
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      DTG(begEnd[0],_DTM.col(begEnd[0]).dot(coefOmegaV));
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      DTG(begEnd[0],_DTM.col(begEnd[0]).dot(coefT));
    //recursion
    if(J._parent>=0)
      toolBRecursiveInplace(k,J._parent,GK);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::DTG(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const
{
  Vec2i begEnd;
  Mat3T coefOmega;
  Vec3T coefOmegaV,coefT;
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    const Joint& J=body.joint(k);
    GETTCM_T(JLT,_TK_1KM,k)
    GETTM_T(GK,G,k)
    //gradient R,C
    coefOmega=ROT(JLT)*ROT(GK).transpose();
    coefT=CTR(GK);
    if(J._parent>=0) {
      coefOmega=coefOmega*ROTI(_TM,J._parent);
      coefT=ROTI(_TM,J._parent).transpose()*coefT;
    }
    coefOmegaV=invCrossMatTrace<T>(coefOmega);
    //assemble
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      DTG[begEnd[0]]+=_DTM.col(begEnd[0]).dot(coefOmegaV);
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      DTG[begEnd[0]]+=_DTM.col(begEnd[0]).dot(coefT);
    //recursion
    if(J._parent>=0)
      toolBRecursive(k,J._parent,G);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::DTGZ(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const
{
  DTG.setZero();
  this->DTG(body,G,DTG);
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::DTGBF(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const
{
  sizeType nrJ=body.nrJ();
  for(sizeType k=0; k<nrJ; k++)
    this->DTG(k,body,G.template block<3,4>(0,k*4),[&](sizeType col,T val) {
    DTG[col]+=val;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::DTGBFZ(const ArticulatedBody& body,Mat3XTM G,VecM DTG) const
{
  DTG.setZero();
  this->DTGBF(body,G,DTG);
}
//-------------------------------------------------------------toolA: this is Ti, MRR,MRt,MtR,Mtt are modified
//this is an interface whose input is 4 3X3 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolA(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3T MRR,Mat3T MRt,Mat3T MtR,Mat3T Mtt,std::function<void(sizeType,sizeType,T)> A) const
{
  Mat3T wK_1KiMwKj,wK_1KiMtKj;
  Mat3T tK_1KiMwKj,tK_1KiMtKj;
  Mat3T wK_1iMwK_1Kj,tK_1iMwK_1Kj;
  Mat3T wK_1iMtK_1Kj,tK_1iMtK_1Kj;
  Mat3T RP,RPj,DC;
  for(; k>=0; k=body.joint(k)._parent) {
    const Joint& J=body.joint(k);
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    RPj=J._parent>=0?ROTI(Tj._TM,J._parent):Mat3T(Mat3T::Identity());
    DC=cross<T>(J._parent>=0?CTRI(_TM,k)-CTRI(_TM,J._parent):Vec3T(Vec3T::Zero()));
    if(J.isRotational()) {
      wK_1KiMwKj=RP.transpose()*MRR;
      wK_1KiMtKj=RP.transpose()*MRt;
      wK_1iMwK_1Kj=(MRR+DC*MtR)*RPj;
      tK_1iMwK_1Kj=MtR*RPj;
      toolANonRecursivePhase1Rotational(k,body,Tj,wK_1KiMwKj,wK_1KiMtKj,A);
      toolANonRecursivePhase2Rotational(k,body,Tj,wK_1iMwK_1Kj,tK_1iMwK_1Kj,A);
    } else {
      tK_1KiMwKj=RP.transpose()*MtR;
      tK_1KiMtKj=RP.transpose()*Mtt;
      wK_1iMtK_1Kj=(MRt+DC*Mtt)*RPj;
      tK_1iMtK_1Kj=Mtt*RPj;
      toolANonRecursivePhase1Translational(k,body,Tj,tK_1KiMwKj,tK_1KiMtKj,A);
      toolANonRecursivePhase2Translational(k,body,Tj,wK_1iMtK_1Kj,tK_1iMtK_1Kj,A);
    }
    //recursion
    if(J._parent>=0)
      toolARecursiveInplace(k,J._parent,Tj,MRR,MRt,MtR,Mtt);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolA(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,std::function<void(sizeType,sizeType,T)> A) const
{
  Mat3T wK_1KiMwKj,wK_1KiMtKj;
  Mat3T tK_1KiMwKj,tK_1KiMtKj;
  Mat3T wK_1iMwK_1Kj,tK_1iMwK_1Kj;
  Mat3T wK_1iMtK_1Kj,tK_1iMtK_1Kj;
  Mat3T RP,RPj,DC;
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    const Joint& J=body.joint(k);
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    RPj=J._parent>=0?ROTI(Tj._TM,J._parent):Mat3T(Mat3T::Identity());
    DC=cross<T>(J._parent>=0?CTRI(_TM,k)-CTRI(_TM,J._parent):Vec3T(Vec3T::Zero()));
    if(J.isRotational()) {
      wK_1KiMwKj=RP.transpose()*MRRI;
      wK_1KiMtKj=RP.transpose()*MRtI;
      wK_1iMwK_1Kj=(MRRI+DC*MtRI)*RPj;
      tK_1iMwK_1Kj=MtRI*RPj;
      toolANonRecursivePhase1Rotational(k,body,Tj,wK_1KiMwKj,wK_1KiMtKj,A);
      toolANonRecursivePhase2Rotational(k,body,Tj,wK_1iMwK_1Kj,tK_1iMwK_1Kj,A);
    } else {
      tK_1KiMwKj=RP.transpose()*MtRI;
      tK_1KiMtKj=RP.transpose()*MttI;
      wK_1iMtK_1Kj=(MRtI+DC*MttI)*RPj;
      tK_1iMtK_1Kj=MttI*RPj;
      toolANonRecursivePhase1Translational(k,body,Tj,tK_1KiMwKj,tK_1KiMtKj,A);
      toolANonRecursivePhase2Translational(k,body,Tj,wK_1iMtK_1Kj,tK_1iMtK_1Kj,A);
    }
    //recursion
    if(J._parent>=0)
      toolARecursive(k,J._parent,Tj,MRR,MRt,MtR,Mtt);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolANonRecursivePhase1Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1KiMwKj,const Mat3T& wK_1KiMtKj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  Tj.JRCSparse(body,k,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(begEnd[0],col,_DTM.col(begEnd[0]).dot(wK_1KiMwKj*v));
  },[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(begEnd[0],col,_DTM.col(begEnd[0]).dot(wK_1KiMtKj*v));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolANonRecursivePhase1Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& tK_1KiMwKj,const Mat3T& tK_1KiMtKj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  Tj.JRCSparse(body,k,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(begEnd[0],col,_DTM.col(begEnd[0]).dot(tK_1KiMwKj*v));
  },[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(begEnd[0],col,_DTM.col(begEnd[0]).dot(tK_1KiMtKj*v));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolANonRecursivePhase2Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMwK_1Kj,const Mat3T& tK_1iMwK_1Kj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  JRCSparse(body,J._parent,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(col,begEnd[0],v.dot(wK_1iMwK_1Kj*Tj._DTM.col(begEnd[0])));
  },[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(col,begEnd[0],v.dot(tK_1iMwK_1Kj*Tj._DTM.col(begEnd[0])));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolANonRecursivePhase2Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMtK_1Kj,const Mat3T& tK_1iMtK_1Kj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  JRCSparse(body,J._parent,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(col,begEnd[0],v.dot(wK_1iMtK_1Kj*Tj._DTM.col(begEnd[0])));
  },[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(col,begEnd[0],v.dot(tK_1iMtK_1Kj*Tj._DTM.col(begEnd[0])));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolARecursive(sizeType k,sizeType p,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt) const
{
  Mat3T DC=cross<T>(CTRI(_TM,k)-CTRI(_TM,p));
  Mat3T DCj=cross<T>(CTRI(Tj._TM,k)-CTRI(Tj._TM,p));
  MRR.template block<3,3>(0,p*3)+=MRRI+MRtI*DCj.transpose()+DC*(MtRI+MttI*DCj.transpose());
  MRt.template block<3,3>(0,p*3)+=MRtI+DC*MttI;
  MtR.template block<3,3>(0,p*3)+=MtRI+MttI*DCj.transpose();
  Mtt.template block<3,3>(0,p*3)+=MttI;
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolARecursiveInplace(sizeType k,sizeType p,const PBDArticulatedGradientInfoMap& Tj,Mat3T& MRR,Mat3T& MRt,Mat3T& MtR,Mat3T& Mtt) const
{
  Mat3T DC=cross<T>(CTRI(_TM,k)-CTRI(_TM,p));
  Mat3T DCj=cross<T>(CTRI(Tj._TM,k)-CTRI(Tj._TM,p));
  MRR+=MRt*DCj.transpose()+DC*(MtR+Mtt*DCj.transpose());
  MRt+=DC*Mtt;
  MtR+=Mtt*DCj.transpose();
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTM A) const
{
  A.setZero();
  toolA(body,Tj,MRR,MRt,MtR,Mtt,[&](sizeType row,sizeType col,T val) {
    A(row,col)+=val;
  });
}
//this is an interface whose input is a 12X12 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAContactAll(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTCM M) const
{
#define CONTRACT(R0,R1,C0,C1) toolAContractTensor<R0,R1,C0,C1>(k,Tj,MRR,MRt,MtR,Mtt,tensor(R0+R1*3,C0+C1*3));
  MRR.setZero();
  MRt.setZero();
  MtR.setZero();
  Mtt.setZero();
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    Eigen::Block<MatTCM,12,12> tensor=M.template block<12,12>(0,k*12);
    CONTRACTRC
  }
#undef CONTRACT
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolA(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const
{
  Mat3XT MRR,MRt,MtR,Mtt;
  MRR.resize(3,body.nrJ()*3);
  MRt.resize(3,body.nrJ()*3);
  MtR.resize(3,body.nrJ()*3);
  Mtt.resize(3,body.nrJ()*3);
  toolAContactAll(body,Tj,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),M);
  toolA(body,Tj,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),[&](sizeType row,sizeType col,T val) {
    A(row,col)+=val;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const
{
  A.setZero();
  toolA(body,Tj,M,A);
}
//this is for internal force's separate term
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolALR(sizeType kL,sizeType kR,const ArticulatedBody& body,const Mat3T& tensor,const Vec3T& ptL,const Vec3T& ptR,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec3T RPtL=ROTI(_TM,kL)*ptL,RPtR=ROTI(_TM,kR)*ptR,DRPtL;
  JRCSparse(body,kL,[&](sizeType row,const Vec3T& JRL) {
    DRPtL=JRL.cross(RPtL);
    JRCSparse(body,kR,[&](sizeType col,const Vec3T& JRR) {
      A(row,col,DRPtL.dot(tensor*JRR.cross(RPtR)));
    },[&](sizeType col,const Vec3T& JCR) {
      A(row,col,DRPtL.dot(tensor*JCR));
    });
  },[&](sizeType row,const Vec3T& JCL) {
    JRCSparse(body,kR,[&](sizeType col,const Vec3T& JRR) {
      A(row,col,JCL.dot(tensor*JRR.cross(RPtR)));
    },[&](sizeType col,const Vec3T& JCR) {
      A(row,col,JCL.dot(tensor*JCR));
    });
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolALR(sizeType kL,sizeType kR,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM tensor,MatTM A) const
{
#define CONTRACT(R0,R1,C0,C1) val+=RT(R0,R1)*RTj(C0,C1)*tensor(R0+R1*3,C0+C1*3);
  sizeType nrDOF=body.nrDOF();
  T val;
  Mat3X4T RT,RTj;
  Mat3XT w,t,wj,tj;
  w.setZero(3,nrDOF);
  t.setZero(3,nrDOF);
  wj.setZero(3,nrDOF);
  tj.setZero(3,nrDOF);
  JRSparse(body,kL,mapM(w));
  JCSparse(body,kL,mapM(t));
  Tj.JRSparse(body,kR,mapM(wj));
  Tj.JCSparse(body,kR,mapM(tj));
  for(sizeType r=0; r<nrDOF; r++)
    for(sizeType c=0; c<nrDOF; c++) {
      RT=concatCol(cross<T>(w.col(r))*ROTI(_TM,kL),t.col(r));
      RTj=concatCol(cross<T>(wj.col(c))*ROTI(Tj._TM,kR),tj.col(c));
      val=0;
      CONTRACTRC
      A(r,c)+=val;
    }
#undef CONTRACT
}
//this is a brute-force interface for testing
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABF(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const
{
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    MatTCM MBlk(&(M.coeffRef(0,k*12)),12,12,Eigen::OuterStride<>(M.outerStride()));
    toolALR(k,k,body,Tj,MBlk,A);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABFZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const
{
  A.setZero();
  toolABF(body,Tj,M,A);
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABF2(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTM A) const
{
  sizeType nrJ=body.nrJ();
  for(sizeType k=0; k<nrJ; k++)
    toolA(k,body,Tj,MRRI,MRtI,MtRI,MttI,[&](sizeType row,sizeType col,T val) {
    A(row,col)+=val;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABF2Z(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,MatTM A) const
{
  A.setZero();
  toolABF2(body,Tj,MRR,MRt,MtR,Mtt,A);
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABF2Z(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM A) const
{
  Mat3XT MRR,MRt,MtR,Mtt;
  MRR.resize(3,body.nrJ()*3);
  MRt.resize(3,body.nrJ()*3);
  MtR.resize(3,body.nrJ()*3);
  Mtt.resize(3,body.nrJ()*3);
  toolAContactAll(body,Tj,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),M);
  toolABF2Z(body,Tj,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),A);
}
//-------------------------------------------------------------toolB: G is modified
//this is an interface whose input is 4 3X3 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolB(const ArticulatedBody& body,Mat3XTM G,std::function<void(sizeType,sizeType,T)> B) const
{
  Mat3T RP,wK_1KiMwKj,tK_1KiMwKj;
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    const Joint& J=body.joint(k);
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    if(J.isRotational()) {
      wK_1KiMwKj=RP.transpose()*invDoubleCrossMatTrace<T>(ROTI(_TM,k)*ROTI(G,k).transpose());
      toolBNonRecursivePhase1Rotational(k,body,*this,wK_1KiMwKj,B);
      toolBNonRecursivePhase2Rotational(k,body,*this,wK_1KiMwKj.transpose(),B);
      JointFunc<T>::DDT(J,B,mapM2CM(_DDTM),RP.transpose()*invCrossMatTrace<T>(ROTI(_TM,k)*ROTI(G,k).transpose()));
    } else {
      tK_1KiMwKj=RP.transpose()*cross<T>(CTRI(G,k));
      toolBNonRecursivePhase1Translational(k,body,*this,tK_1KiMwKj,B);
      toolBNonRecursivePhase2Translational(k,body,*this,tK_1KiMwKj.transpose(),B);
    }
    //recursion
    if(J._parent>=0)
      toolBRecursive(k,J._parent,G);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBNonRecursivePhase1Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1KiMwKj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  Tj.JRSparse(body,k,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(begEnd[0],col,_DTM.col(begEnd[0]).dot(wK_1KiMwKj*v));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBNonRecursivePhase1Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& tK_1KiMwKj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  Tj.JRSparse(body,k,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(begEnd[0],col,_DTM.col(begEnd[0]).dot(tK_1KiMwKj*v));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBNonRecursivePhase2Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMwK_1Kj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  JRSparse(body,J._parent,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(col,begEnd[0],v.dot(wK_1iMwK_1Kj*Tj._DTM.col(begEnd[0])));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBNonRecursivePhase2Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iMtK_1Kj,std::function<void(sizeType,sizeType,T)> A) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  JRSparse(body,J._parent,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      A(col,begEnd[0],v.dot(wK_1iMtK_1Kj*Tj._DTM.col(begEnd[0])));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBRecursive(sizeType k,sizeType p,Mat3XTM G) const
{
  GETTCM_T(JLT,_TK_1KM,k)
  GETTM_T(GK,G,k)
  GETTM_T(GKP,G,p)
  ROT(GKP)+=ROT(GK)*ROT(JLT).transpose()+CTR(GK)*CTR(JLT).transpose();
  CTR(GKP)+=CTR(GK);
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBRecursiveInplace(sizeType k,sizeType p,Mat3X4T& G) const
{
  GETTCM_T(JLT,_TK_1KM,k)
  ROT(G)=ROT(G)*ROT(JLT).transpose()+CTR(G)*CTR(JLT).transpose();
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBZ(const ArticulatedBody& body,Mat3XTM G,MatTM B) const
{
  B.setZero();
  toolB(body,G,[&](sizeType row,sizeType col,T val) {
    B(row,col)+=val;
  });
}
//this is an interface whose input is a 12X12 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBContactAll(const ArticulatedBody& body,MatTCM M,Mat3XTM G) const
{
#define CONTRACT(R0,R1,C0,C1) toolBContractTensor<R0,R1,C0,C1>(k,G,tensor(R0+R1*3,C0+C1*3));
  G.setZero();
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    Eigen::Block<MatTCM,12,12> tensor=M.template block<12,12>(0,k*12);
    CONTRACTRC
  }
#undef CONTRACT
}
//this is a brute-force interface for testing
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBBF(const ArticulatedBody& body,MatTCM M,MatTM B) const
{
#define CONTRACT(R0,R1,C0,C1) G(C0,C1+k*4)+=tensor(R0+R1*3,C0+C1*3)*TG1(R0,R1);
  Mat3XT G=Mat3XT::Zero(3,body.nrJ()*4);
  for(sizeType k=body.nrJ()-1; k>=0; k--) {
    Mat3X4T TG1=concatCol(cross<T>(DWDLI(_DTLambdaM,k))*ROTI(_TM,k),DTDLI(_DTLambdaM,k));
    Eigen::Block<MatTCM,12,12> tensor=M.template block<12,12>(0,k*12);
    CONTRACTRC
  }
  toolB(body,mapM(G),[&](sizeType row,sizeType col,T val) {
    B(row,col)+=val;
  });
#undef CONTRACT
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBBFZ(const ArticulatedBody& body,MatTCM M,MatTM B) const
{
  B.setZero();
  toolBBF(body,M,B);
}
//-------------------------------------------------------------toolA,toolB combined
//this is an interface whose input is 4 3X3 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAB(sizeType k,const ArticulatedBody& body,Mat3T MRR,Mat3T MRt,Mat3T MtR,Mat3T Mtt,Mat3X4T G,std::function<void(sizeType row,sizeType col,T val)> AB) const
{
  Mat3T wK_1KiMwKj,wK_1KiMtKj;
  Mat3T tK_1KiMwKj,tK_1KiMtKj;
  Mat3T wK_1iMwK_1Kj,tK_1iMwK_1Kj;
  Mat3T wK_1iMtK_1Kj,tK_1iMtK_1Kj;
  Mat3T RP,BContrib,DC;
  for(; k>=0; k=body.joint(k)._parent) {
    const Joint& J=body.joint(k);
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    DC=cross<T>(J._parent>=0?CTRI(_TM,k)-CTRI(_TM,J._parent):Vec3T(Vec3T::Zero()));
    if(J.isRotational()) {
      BContrib=RP.transpose()*invDoubleCrossMatTrace<T>(ROTI(_TM,k)*ROT(G).transpose());
      wK_1KiMwKj=RP.transpose()*MRR+BContrib;
      wK_1KiMtKj=RP.transpose()*MRt;
      wK_1iMwK_1Kj=(MRR+DC*MtR)*RP+BContrib.transpose();
      tK_1iMwK_1Kj=MtR*RP;
      toolANonRecursivePhase1Rotational(k,body,*this,wK_1KiMwKj,wK_1KiMtKj,AB);
      toolANonRecursivePhase2Rotational(k,body,*this,wK_1iMwK_1Kj,tK_1iMwK_1Kj,AB);
      JointFunc<T>::DDT(J,AB,mapM2CM(_DDTM),RP.transpose()*invCrossMatTrace<T>(ROTI(_TM,k)*ROT(G).transpose()));
    } else {
      BContrib=RP.transpose()*cross<T>(CTR(G));
      tK_1KiMwKj=RP.transpose()*MtR+BContrib;
      tK_1KiMtKj=RP.transpose()*Mtt;
      wK_1iMtK_1Kj=(MRt+DC*Mtt)*RP+BContrib.transpose();
      tK_1iMtK_1Kj=Mtt*RP;
      toolANonRecursivePhase1Translational(k,body,*this,tK_1KiMwKj,tK_1KiMtKj,AB);
      toolANonRecursivePhase2Translational(k,body,*this,wK_1iMtK_1Kj,tK_1iMtK_1Kj,AB);
    }
    //recursion
    if(J._parent>=0) {
      toolARecursiveInplace(k,J._parent,*this,MRR,MRt,MtR,Mtt);
      toolBRecursiveInplace(k,J._parent,G);
    }
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAB(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,std::function<void(sizeType row,sizeType col,T val)> AB) const
{
  Mat3T wK_1KiMwKj,wK_1KiMtKj;
  Mat3T tK_1KiMwKj,tK_1KiMtKj;
  Mat3T wK_1iMwK_1Kj,tK_1iMwK_1Kj;
  Mat3T wK_1iMtK_1Kj,tK_1iMtK_1Kj;
  Mat3T RP,BContrib,DC;
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    const Joint& J=body.joint(k);
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    DC=cross<T>(J._parent>=0?CTRI(_TM,k)-CTRI(_TM,J._parent):Vec3T(Vec3T::Zero()));
    if(J.isRotational()) {
      BContrib=RP.transpose()*invDoubleCrossMatTrace<T>(ROTI(_TM,k)*ROTI(G,k).transpose());
      wK_1KiMwKj=RP.transpose()*MRRI+BContrib;
      wK_1KiMtKj=RP.transpose()*MRtI;
      wK_1iMwK_1Kj=(MRRI+DC*MtRI)*RP+BContrib.transpose();
      tK_1iMwK_1Kj=MtRI*RP;
      toolANonRecursivePhase1Rotational(k,body,*this,wK_1KiMwKj,wK_1KiMtKj,AB);
      toolANonRecursivePhase2Rotational(k,body,*this,wK_1iMwK_1Kj,tK_1iMwK_1Kj,AB);
      JointFunc<T>::DDT(J,AB,mapM2CM(_DDTM),RP.transpose()*invCrossMatTrace<T>(ROTI(_TM,k)*ROTI(G,k).transpose()));
    } else {
      BContrib=RP.transpose()*cross<T>(CTRI(G,k));
      tK_1KiMwKj=RP.transpose()*MtRI+BContrib;
      tK_1KiMtKj=RP.transpose()*MttI;
      wK_1iMtK_1Kj=(MRtI+DC*MttI)*RP+BContrib.transpose();
      tK_1iMtK_1Kj=MttI*RP;
      toolANonRecursivePhase1Translational(k,body,*this,tK_1KiMwKj,tK_1KiMtKj,AB);
      toolANonRecursivePhase2Translational(k,body,*this,wK_1iMtK_1Kj,tK_1iMtK_1Kj,AB);
    }
    //recursion
    if(J._parent>=0) {
      toolARecursive(k,J._parent,*this,MRR,MRt,MtR,Mtt);
      toolBRecursive(k,J._parent,G);
    }
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABZ(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,MatTM AB) const
{
  AB.setZero();
  toolAB(body,MRR,MRt,MtR,Mtt,G,[&](sizeType row,sizeType col,T val) {
    AB(row,col)+=val;
  });
}
//this is an interface whose input is a 12X12 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAB(const ArticulatedBody& body,MatTCM M,Mat3XTM G,MatTM AB) const
{
  Mat3XT MRR,MRt,MtR,Mtt;
  MRR.resize(3,body.nrJ()*3);
  MRt.resize(3,body.nrJ()*3);
  MtR.resize(3,body.nrJ()*3);
  Mtt.resize(3,body.nrJ()*3);
  toolAContactAll(body,*this,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),M);
  toolAB(body,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),G,[&](sizeType row,sizeType col,T val) {
    AB(row,col)+=val;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABZ(const ArticulatedBody& body,MatTCM M,Mat3XTM G,MatTM AB) const
{
  AB.setZero();
  toolAB(body,M,G,AB);
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABBF(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,MatTM AB) const
{
  sizeType nrJ=body.nrJ();
  for(sizeType k=0; k<nrJ; k++)
    toolAB(k,body,MRRI,MRtI,MtRI,MttI,G.template block<3,4>(0,k*4),[&](sizeType row,sizeType col,T val) {
    AB(row,col)+=val;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABBFZ(const ArticulatedBody& body,Mat3XTM MRR,Mat3XTM MRt,Mat3XTM MtR,Mat3XTM Mtt,Mat3XTM G,MatTM AB) const
{
  AB.setZero();
  toolABBF(body,MRR,MRt,MtR,Mtt,G,AB);
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABBFZ(const ArticulatedBody& body,MatTCM M,Mat3XTM G,MatTM AB) const
{
  Mat3XT MRR,MRt,MtR,Mtt;
  MRR.resize(3,body.nrJ()*3);
  MRt.resize(3,body.nrJ()*3);
  MtR.resize(3,body.nrJ()*3);
  Mtt.resize(3,body.nrJ()*3);
  toolAContactAll(body,*this,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),M);
  toolABBFZ(body,mapM(MRR),mapM(MRt),mapM(MtR),mapM(Mtt),G,AB);
}
//-------------------------------------------------------------toolA,toolC combined
//this is an interface whose input is 4 3X3 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAC(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,std::function<void(sizeType row,sizeType col,T val)> AC) const
{
  Mat3T wK_1KiMwKj,wK_1KiMtKj,wK_1KiLambdaMwKj,wK_1KiLambdaMtKj;
  Mat3T tK_1KiMwKj,tK_1KiMtKj;
  Mat3T wK_1iMwK_1Kj,tK_1iMwK_1Kj;
  Mat3T wK_1iMtK_1Kj,tK_1iMtK_1Kj;
  Mat3T wK_1iLambdaMwK_1Kj,tK_1iLambdaMwK_1Kj;
  Mat3T wK_1iLambdaMtK_1Kj,tK_1iLambdaMtK_1Kj;
  Mat3T RP,RPj,DC,DCj,DWDLP,DDTDLP;
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    const Joint& J=body.joint(k);
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    RPj=J._parent>=0?ROTI(Tj._TM,J._parent):Mat3T(Mat3T::Identity());
    DC=cross<T>(J._parent>=0?CTRI(_TM,k)-CTRI(_TM,J._parent):Vec3T(Vec3T::Zero()));
    DCj=cross<T>(J._parent>=0?CTRI(Tj._TM,k)-CTRI(Tj._TM,J._parent):Vec3T(Vec3T::Zero()));
    DWDLP=cross<T>(J._parent>=0?DWDLI(_DTLambdaM,J._parent):Vec3T(Vec3T::Zero()));
    DDTDLP=cross<T>(J._parent>=0?DTDLI(_DTLambdaM,k)-DTDLI(_DTLambdaM,J._parent):Vec3T(Vec3T::Zero()));
    if(J.isRotational()) {
      wK_1KiMwKj=RP.transpose()*(MRRAI+DWDLP.transpose()*MRRCI);
      wK_1KiMtKj=RP.transpose()*(MRtAI+DWDLP.transpose()*MRtCI);
      wK_1KiLambdaMwKj=RP.transpose()*MRRCI;
      wK_1KiLambdaMtKj=RP.transpose()*MRtCI;
      wK_1iMwK_1Kj=(MRRAI+DC*MtRAI+DDTDLP*MtRCI)*RPj;
      tK_1iMwK_1Kj=MtRAI*RPj;
      wK_1iLambdaMwK_1Kj=(MRRCI+DC*MtRCI)*RPj;
      tK_1iLambdaMwK_1Kj=MtRCI*RPj;
      toolACNonRecursivePhase1Rotational(k,body,Tj,wK_1KiMwKj,wK_1KiMtKj,wK_1KiLambdaMwKj,wK_1KiLambdaMtKj,AC);
      toolANonRecursivePhase2Rotational(k,body,Tj,wK_1iMwK_1Kj,tK_1iMwK_1Kj,AC);
      toolACNonRecursivePhase2Rotational(k,body,Tj,wK_1iLambdaMwK_1Kj,tK_1iLambdaMwK_1Kj,AC);
    } else {
      tK_1KiMwKj=RP.transpose()*(MtRAI+DWDLP.transpose()*MtRCI);
      tK_1KiMtKj=RP.transpose()*(MttAI+DWDLP.transpose()*MttCI);
      wK_1iMtK_1Kj=(MRtAI+DC*MttAI+DDTDLP*MttCI)*RPj;
      tK_1iMtK_1Kj=MttAI*RPj;
      wK_1iLambdaMtK_1Kj=(MRtCI+DC*MttCI)*RPj;
      tK_1iLambdaMtK_1Kj=MttCI*RPj;
      toolANonRecursivePhase1Translational(k,body,Tj,tK_1KiMwKj,tK_1KiMtKj,AC);
      toolANonRecursivePhase2Translational(k,body,Tj,wK_1iMtK_1Kj,tK_1iMtK_1Kj,AC);
      toolACNonRecursivePhase2Translational(k,body,Tj,wK_1iLambdaMtK_1Kj,tK_1iLambdaMtK_1Kj,AC);
    }
    //recursion
    if(J._parent>=0)
      toolACRecursive(k,J._parent,Tj,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolACNonRecursivePhase1Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1KiMwKj,const Mat3T& wK_1KiMtKj,const Mat3T& wK_1KiLambdaMwKj,const Mat3T& wK_1KiLambdaMtKj,std::function<void(sizeType row,sizeType col,T val)> AC) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  Tj.JRCSparse(body,k,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      AC(begEnd[0],col,_DTM.col(begEnd[0]).dot(wK_1KiMwKj*v)+_DTILambdaM.col(begEnd[0]).dot(wK_1KiLambdaMwKj*v));
  },[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      AC(begEnd[0],col,_DTM.col(begEnd[0]).dot(wK_1KiMtKj*v)+_DTILambdaM.col(begEnd[0]).dot(wK_1KiLambdaMtKj*v));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolACNonRecursivePhase2Rotational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iLambdaMwK_1Kj,const Mat3T& tK_1iLambdaMwK_1Kj,std::function<void(sizeType row,sizeType col,T val)> AC) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  JRCILambdaSparse(body,J._parent,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      AC(col,begEnd[0],v.dot(wK_1iLambdaMwK_1Kj*Tj._DTM.col(begEnd[0])));
  },[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      AC(col,begEnd[0],v.dot(tK_1iLambdaMwK_1Kj*Tj._DTM.col(begEnd[0])));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolACNonRecursivePhase2Translational(sizeType k,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,const Mat3T& wK_1iLambdaMtK_1Kj,const Mat3T& tK_1iLambdaMtK_1Kj,std::function<void(sizeType row,sizeType col,T val)> AC) const
{
  Vec2i begEnd;
  const Joint& J=body.joint(k);
  JRCILambdaSparse(body,J._parent,[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      AC(col,begEnd[0],v.dot(wK_1iLambdaMtK_1Kj*Tj._DTM.col(begEnd[0])));
  },[&](sizeType col,const Vec3T& v) {
    for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++)
      AC(col,begEnd[0],v.dot(tK_1iLambdaMtK_1Kj*Tj._DTM.col(begEnd[0])));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolACRecursive(sizeType k,sizeType p,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC) const
{
  Mat3T DC=cross<T>(CTRI(_TM,k)-CTRI(_TM,p));
  Mat3T DCj=cross<T>(CTRI(Tj._TM,k)-CTRI(Tj._TM,p));
  Mat3T DDTDLP=cross<T>(DTDLI(_DTLambdaM,k)-DTDLI(_DTLambdaM,p));
  MRRC.template block<3,3>(0,p*3)+=MRRCI+MRtCI*DCj.transpose()+DC*(MtRCI+MttCI*DCj.transpose());
  MRtC.template block<3,3>(0,p*3)+=MRtCI+DC*MttCI;
  MtRC.template block<3,3>(0,p*3)+=MtRCI+MttCI*DCj.transpose();
  MttC.template block<3,3>(0,p*3)+=MttCI;
  MRRA.template block<3,3>(0,p*3)+=DDTDLP*(MtRCI+MttCI*DCj.transpose());
  MRtA.template block<3,3>(0,p*3)+=DDTDLP*MttCI;
  toolARecursive(k,p,Tj,MRRA,MRtA,MtRA,MttA);
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolACZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTM AC) const
{
  AC.setZero();
  toolAC(body,Tj,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,[&](sizeType row,sizeType col,T val) {
    AC(row,col)+=val;
  });
}
//this is an interface whose input is a 12X12 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolCContactAll(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC) const
{
#define CONTRACT(R0,R1,C0,C1) toolCContractTensor<R0,R1,C0,C1>(k,Tj,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,tensor(R0+R1*3,C0+C1*3));
  MRRA.setZero();
  MRtA.setZero();
  MtRA.setZero();
  MttA.setZero();
  MRRC.setZero();
  MRtC.setZero();
  MtRC.setZero();
  MttC.setZero();
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    Eigen::Block<MatTCM,12,12> tensor=MC.template block<12,12>(0,k*12);
    CONTRACTRC
  }
#undef CONTRACT
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolACContactAll(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,MatTCM MA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC) const
{
#define CONTRACT(R0,R1,C0,C1) toolACContractTensor<R0,R1,C0,C1>(k,Tj,MRRA,MRtA,MtRA,MttA,tensorA(R0+R1*3,C0+C1*3),MRRC,MRtC,MtRC,MttC,tensorC(R0+R1*3,C0+C1*3));
  MRRA.setZero();
  MRtA.setZero();
  MtRA.setZero();
  MttA.setZero();
  MRRC.setZero();
  MRtC.setZero();
  MtRC.setZero();
  MttC.setZero();
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    Eigen::Block<MatTCM,12,12> tensorA=MA.template block<12,12>(0,k*12);
    Eigen::Block<MatTCM,12,12> tensorC=MC.template block<12,12>(0,k*12);
    CONTRACTRC
  }
#undef CONTRACT
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolAC(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM MA,MatTCM MC,MatTM AC) const
{
  Mat3XT MRRA,MRtA,MtRA,MttA;
  Mat3XT MRRC,MRtC,MtRC,MttC;
  MRRA.resize(3,body.nrJ()*3);
  MRtA.resize(3,body.nrJ()*3);
  MtRA.resize(3,body.nrJ()*3);
  MttA.resize(3,body.nrJ()*3);
  MRRC.resize(3,body.nrJ()*3);
  MRtC.resize(3,body.nrJ()*3);
  MtRC.resize(3,body.nrJ()*3);
  MttC.resize(3,body.nrJ()*3);
  toolACContactAll(body,Tj,mapM(MRRA),mapM(MRtA),mapM(MtRA),mapM(MttA),MA,mapM(MRRC),mapM(MRtC),mapM(MtRC),mapM(MttC),MC);
  toolAC(body,Tj,mapM(MRRA),mapM(MRtA),mapM(MtRA),mapM(MttA),mapM(MRRC),mapM(MRtC),mapM(MtRC),mapM(MttC),[&](sizeType row,sizeType col,T val) {
    AC(row,col)+=val;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolACZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM MA,MatTCM MC,MatTM AC) const
{
  AC.setZero();
  toolAC(body,Tj,MA,MC,AC);
}
//this is for internal force's separate term
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolCLR(sizeType kL,sizeType kR,const ArticulatedBody& body,const Mat3T& tensor,const Vec3T& ptL,const Vec3T& ptR,std::function<void(sizeType row,sizeType col,T val)> C) const
{
  Vec3T RPtL=ROTI(_TM,kL)*ptL;
  Vec3T CLRPtL=DWDLI(_DTLambdaM,kL).cross(RPtL);
  Vec3T RPtR=ROTI(_TM,kR)*ptR,DRPtL;
  JRCILambdaSparse(body,kL,[&](sizeType row,const Vec3T& JRL) {
    DRPtL=JRL.cross(RPtL)+_RDTM.col(row).cross(CLRPtL);
    JRCSparse(body,kR,[&](sizeType col,const Vec3T& JRR) {
      C(row,col,DRPtL.dot(tensor*JRR.cross(RPtR)));
    },[&](sizeType col,const Vec3T& JCR) {
      C(row,col,DRPtL.dot(tensor*JCR));
    });
  },[&](sizeType row,const Vec3T& JCL) {
    JRCSparse(body,kR,[&](sizeType col,const Vec3T& JRR) {
      C(row,col,JCL.dot(tensor*JRR.cross(RPtR)));
    },[&](sizeType col,const Vec3T& JCR) {
      C(row,col,JCL.dot(tensor*JCR));
    });
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolCLR(sizeType kL,sizeType kR,const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM tensor,MatTM A) const
{
#define CONTRACT(R0,R1,C0,C1) val+=RT(R0,R1)*RTj(C0,C1)*tensor(R0+R1*3,C0+C1*3);
  sizeType nrDOF=body.nrDOF();
  T val;
  Mat3X4T RT,RTj;
  Mat3XT w,wIL,tIL,wj,tj;
  w.setZero(3,nrDOF);
  wIL.setZero(3,nrDOF);
  tIL.setZero(3,nrDOF);
  wj.setZero(3,nrDOF);
  tj.setZero(3,nrDOF);
  JRSparse(body,kL,mapM(w));
  JRILambdaSparse(body,kL,mapM(wIL));
  JCILambdaSparse(body,kL,mapM(tIL));
  Tj.JRSparse(body,kR,mapM(wj));
  Tj.JCSparse(body,kR,mapM(tj));
  for(sizeType r=0; r<nrDOF; r++)
    for(sizeType c=0; c<nrDOF; c++) {
      RT=concatCol((cross<T>(wIL.col(r))+cross<T>(w.col(r))*cross<T>(DWDLI(_DTLambdaM,kL)))*ROTI(_TM,kL),tIL.col(r));
      RTj=concatCol(cross<T>(wj.col(c))*ROTI(Tj._TM,kR),tj.col(c));
      val=0;
      CONTRACTRC
      A(r,c)+=val;
    }
#undef CONTRACT
}
//this is a brute-force interface for testing
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolCBF(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM C) const
{
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    MatTCM MBlk(&(M.coeffRef(0,k*12)),12,12,Eigen::OuterStride<>(M.outerStride()));
    toolCLR(k,k,body,Tj,MBlk,C);
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolCBFZ(const ArticulatedBody& body,const PBDArticulatedGradientInfoMap& Tj,MatTCM M,MatTM C) const
{
  C.setZero();
  toolCBF(body,Tj,M,C);
}
//-------------------------------------------------------------toolD using brute-force method for testing
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolDBF(const ArticulatedBody& body,Mat3XTM G,MatTM D) const
{
  Vec2i begEnd;
  sizeType nrJ=body.nrJ();
  Mat3T RP,C,wK_1KiMwKj,tK_1KiMwKj;
  Mat3XT GB=Mat3XT::Zero(3,nrJ*4);
  for(sizeType k=nrJ-1; k>=0; k--) {
    const Joint& J=body.joint(k);
    //term1: third order derivatives of T_{K-1}^K
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    HRK_1KLambdaSparse(body,k,[&](sizeType row,sizeType col,const Mat3T& HR) {
      D(row,col)+=(RP*HR*ROTI(G,k).transpose()).trace();
    });
    //term2: absorb into toolB
    if(J.isRotational()) {
      wK_1KiMwKj=RP.transpose()*invDoubleCrossMatTrace<T>(ROTI(_TM,k)*ROTI(GB,k).transpose());
      toolBNonRecursivePhase1Rotational(k,body,*this,wK_1KiMwKj,[&](sizeType row,sizeType col,T val) {
        D(row,col)+=val;
      });
      toolBNonRecursivePhase2Rotational(k,body,*this,wK_1KiMwKj.transpose(),[&](sizeType row,sizeType col,T val) {
        D(row,col)+=val;
      });
      JointFunc<T>::DDT(J,D,mapM2CM(_DDTM),RP.transpose()*invCrossMatTrace<T>(ROTI(_TM,k)*ROTI(GB,k).transpose()));
    } else {
      tK_1KiMwKj=RP.transpose()*cross<T>(CTRI(GB,k));
      toolBNonRecursivePhase1Translational(k,body,*this,tK_1KiMwKj,[&](sizeType row,sizeType col,T val) {
        D(row,col)+=val;
      });
      toolBNonRecursivePhase2Translational(k,body,*this,tK_1KiMwKj.transpose(),[&](sizeType row,sizeType col,T val) {
        D(row,col)+=val;
      });
    }
    if(J._parent>=0) {
      //term3: second order derivatives of T_{K-1}^K
      HRK_1KSparse(body,k,[&](sizeType row,sizeType col,const Mat3T& HR) {
        D(row,col)+=(cross<T>(DWDLI(_DTLambdaM,J._parent))*RP*HR*ROTI(G,k).transpose()).trace();
      });
      //term4: lambdaI|J, lambdaJ|I
      JRMILambdaSparse(body,J._parent,[&](sizeType col,const Mat3T& R) {
        for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++) {
          T val=(R*cross<T>(_DTM.col(begEnd[0]))*ROTI(_TK_1KM,k)*ROTI(G,k).transpose()).trace();
          D(begEnd[0],col)+=val;
          D(col,begEnd[0])+=val;
        }
        for(begEnd=J.CBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++) {
          T val=(R*_DTM.col(begEnd[0])).dot(CTRI(G,k));
          D(begEnd[0],col)+=val;
          D(col,begEnd[0])+=val;
        }
      });
      //term5: I|lambdaJ, J|lambdaI
      JRMSparse(body,J._parent,[&](sizeType col,const Mat3T& R) {
        for(begEnd=J.RBegEnd(); begEnd[0]<begEnd[1]; begEnd[0]++) {
          C=cross<T>(_DTM.col(begEnd[0]))*cross<T>(DWDLI(_DTK_1KLambdaM,k))+cross<T>(_DTILambdaM.col(begEnd[0]));
          T val=(R*C*ROTI(_TK_1KM,k)*ROTI(G,k).transpose()).trace();
          D(begEnd[0],col)+=val;
          D(col,begEnd[0])+=val;
        }
      });
      //recursion
      toolDRecursive(k,J._parent,G,mapM(GB),1);
      toolBRecursive(k,J._parent,mapM(GB));
      toolBRecursive(k,J._parent,G);
    }
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolDRecursive(sizeType k,sizeType p,Mat3XTM G,Mat3XTM GB,T coef) const
{
  Mat3T DRLLambda=cross<T>(DWDLI(_DTK_1KLambdaM,k))*ROTI(_TK_1KM,k);
  ROTI(GB,p)+=(ROTI(G,k)*DRLLambda.transpose()+CTRI(G,k)*DTDLI(_DTK_1KLambdaM,k).transpose())*coef;
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolDBFZ(const ArticulatedBody& body,Mat3XTM G,MatTM D) const
{
  D.setZero();
  toolDBF(body,G,D);
}
//-------------------------------------------------------------toolA,B,C,CT,D combined
//this is an interface whose input is 4 3X3 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABCCTD(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,Mat3XTM G,Mat3XTM GB,std::function<void(sizeType row,sizeType col,T val)> ABCCTD) const
{
  GB/=2;
  Mat3T wK_1KiMwKjB,tK_1KiMwKjB;
  Mat3T wK_1KiMwKj,wK_1KiMtKj,wK_1KiLambdaMwKj,wK_1KiLambdaMtKj;
  Mat3T tK_1KiMwKj,tK_1KiMtKj;
  Mat3T wK_1iMwK_1Kj,tK_1iMwK_1Kj;
  Mat3T wK_1iMtK_1Kj,tK_1iMtK_1Kj;
  Mat3T wK_1iLambdaMwK_1Kj,tK_1iLambdaMwK_1Kj;
  Mat3T wK_1iLambdaMtK_1Kj,tK_1iLambdaMtK_1Kj;
  Mat3T RP,CDWDLP,CDWDLL,DC,DDTDLP;
  Vec3T DWDLP,DWDLL;
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    const Joint& J=body.joint(k);
    RP=J._parent>=0?ROTI(_TM,J._parent):Mat3T(Mat3T::Identity());
    DC=cross<T>(J._parent>=0?CTRI(_TM,k)-CTRI(_TM,J._parent):Vec3T(Vec3T::Zero()));
    DWDLP=J._parent>=0?DWDLI(_DTLambdaM,J._parent):Vec3T(Vec3T::Zero());
    DWDLL=DWDLI(_DTLambdaM,k)-DWDLP;
    CDWDLP=cross<T>(DWDLP);
    CDWDLL=cross<T>(DWDLL);
    DDTDLP=cross<T>(J._parent>=0?DTDLI(_DTLambdaM,k)-DTDLI(_DTLambdaM,J._parent):Vec3T(Vec3T::Zero()));
    if(J.isRotational()) {
      wK_1KiMwKj=RP.transpose()*(MRRAI+CDWDLP.transpose()*MRRCI);
      wK_1KiMwKjB=RP.transpose()*invDoubleCrossMatTrace<T>(ROTI(_TM,k)*ROTI(GB,k).transpose());
      wK_1KiMtKj=RP.transpose()*(MRtAI+CDWDLP.transpose()*MRtCI);
      wK_1KiLambdaMwKj=RP.transpose()*(MRRCI+invDoubleCrossMatTrace<T>(ROTI(_TM,k)*ROTI(G,k).transpose()));
      wK_1KiLambdaMtKj=RP.transpose()*MRtCI;
      wK_1iMwK_1Kj=(MRRAI+DC*MtRAI+DDTDLP*MtRCI+invDoubleCrossMatTrace<T>(DWDLP,ROTI(_TM,k)*ROTI(G,k).transpose()).transpose()-
                    invDoubleCrossMatTrace<T>(ROTI(G,k)*ROTI(_TM,k).transpose()*CDWDLL))*RP;
      tK_1iMwK_1Kj=MtRAI*RP;
      wK_1iLambdaMwK_1Kj=(MRRCI+DC*MtRCI+invDoubleCrossMatTrace<T>(ROTI(G,k)*ROTI(_TM,k).transpose()))*RP;
      tK_1iLambdaMwK_1Kj=MtRCI*RP;
      toolACNonRecursivePhase1Rotational(k,body,*this,wK_1KiMwKj+wK_1KiMwKjB,wK_1KiMtKj,wK_1KiLambdaMwKj,wK_1KiLambdaMtKj,ABCCTD);
      toolANonRecursivePhase2Rotational(k,body,*this,wK_1iMwK_1Kj+wK_1KiMwKjB.transpose(),tK_1iMwK_1Kj,ABCCTD);
      toolACNonRecursivePhase2Rotational(k,body,*this,wK_1iLambdaMwK_1Kj,tK_1iLambdaMwK_1Kj,ABCCTD);
      //diagonal part
      Mat3T& RGT=wK_1KiMwKj=ROTI(_TM,k)*ROTI(G,k).transpose()/2;
      Mat3T& CWRGT_RGTCWP=wK_1KiMwKjB=CDWDLL*RGT+RGT*CDWDLP;
      Mat3T& IDC_RGT=wK_1KiMtKj=invDoubleCrossMatTrace<T>(RGT);
      Mat3T& IDC_CWRGT_RGTCWP=wK_1KiLambdaMwKj=invDoubleCrossMatTrace<T>(CWRGT_RGTCWP);
      for(Vec2i RI=J.RBegEnd(); RI[0]<RI[1]; RI[0]++)
        for(Vec2i RJ=J.RBegEnd(); RJ[0]<RJ[1]; RJ[0]++)
          ABCCTD(RI[0],RJ[0],_RDTM.col(RJ[0]).dot(IDC_RGT*_RDTILambdaM.col(RI[0]))-_RDTILambdaM.col(RJ[0]).dot(IDC_RGT*_RDTM.col(RI[0]))+_RDTM.col(RJ[0]).dot(IDC_CWRGT_RGTCWP*_RDTM.col(RI[0])));
      JointFunc<T>::DDT(J,ABCCTD,mapM2CM(_DTIILambdaM),RP.transpose()*invCrossMatTrace<T>(RGT));
      JointFunc<T>::DDT(J,ABCCTD,mapM2CM(_DDTM),RP.transpose()*invCrossMatTrace<T>(CWRGT_RGTCWP+ROTI(_TM,k)*ROTI(GB,k).transpose()));
    } else {
      tK_1KiMwKj=RP.transpose()*(MtRAI+CDWDLP.transpose()*MtRCI);
      tK_1KiMwKjB=RP.transpose()*cross<T>(CTRI(GB,k));
      tK_1KiMtKj=RP.transpose()*(MttAI+CDWDLP.transpose()*MttCI);
      wK_1iMtK_1Kj=(MRtAI+DC*MttAI+DDTDLP*MttCI-cross<T>(CTRI(G,k))*CDWDLP)*RP;
      tK_1iMtK_1Kj=MttAI*RP;
      wK_1iLambdaMtK_1Kj=(MRtCI+DC*MttCI-cross<T>(CTRI(G,k)))*RP;
      tK_1iLambdaMtK_1Kj=MttCI*RP;
      toolANonRecursivePhase1Translational(k,body,*this,tK_1KiMwKj+tK_1KiMwKjB,tK_1KiMtKj,ABCCTD);
      toolANonRecursivePhase2Translational(k,body,*this,wK_1iMtK_1Kj+tK_1KiMwKjB.transpose(),tK_1iMtK_1Kj,ABCCTD);
      toolACNonRecursivePhase2Translational(k,body,*this,wK_1iLambdaMtK_1Kj,tK_1iLambdaMtK_1Kj,ABCCTD);
    }
    if(J._parent>=0) {
      //recursion
      toolDRecursive(k,J._parent,G,GB,0.5f);
      toolBRecursive(k,J._parent,G);
      toolBRecursive(k,J._parent,GB);
      toolACRecursive(k,J._parent,*this,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC);
    }
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABCCTDZ(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,Mat3XTM G,Mat3XTM GB,MatTM ABCCTD) const
{
  ABCCTD.setZero();
  toolABCCTD(body,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,G,GB,[&](sizeType row,sizeType col,T val) {
    ABCCTD(row,col)+=val;
    ABCCTD(col,row)+=val;
  });
}
//this is an interface whose input is a 12X12 tensor
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolBCCTDContactAll(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC,Mat3XTM GB) const
{
#define CONTRACT(R0,R1,C0,C1) toolBCCTDContractTensor<R0,R1,C0,C1>(k,MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,tensor(R0+R1*3,C0+C1*3),GB);
  MRRA.setZero();
  MRtA.setZero();
  MtRA.setZero();
  MttA.setZero();
  MRRC.setZero();
  MRtC.setZero();
  MtRC.setZero();
  MttC.setZero();
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    Eigen::Block<MatTCM,12,12> tensor=MC.template block<12,12>(0,k*12);
    CONTRACTRC
  }
#undef CONTRACT
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABCCTDContactAll(const ArticulatedBody& body,Mat3XTM MRRA,Mat3XTM MRtA,Mat3XTM MtRA,Mat3XTM MttA,MatTCM MA,Mat3XTM MRRC,Mat3XTM MRtC,Mat3XTM MtRC,Mat3XTM MttC,MatTCM MC,Mat3XTM GB) const
{
#define CONTRACT(R0,R1,C0,C1) toolABCCTDContractTensor<R0,R1,C0,C1>(k,MRRA,MRtA,MtRA,MttA,tensorA(R0+R1*3,C0+C1*3),MRRC,MRtC,MtRC,MttC,tensorC(R0+R1*3,C0+C1*3),GB);
  MRRA.setZero();
  MRtA.setZero();
  MtRA.setZero();
  MttA.setZero();
  MRRC.setZero();
  MRtC.setZero();
  MtRC.setZero();
  MttC.setZero();
  sizeType nrJ=body.nrJ();
  for(sizeType k=nrJ-1; k>=0; k--) {
    Eigen::Block<MatTCM,12,12> tensorA=MA.template block<12,12>(0,k*12);
    Eigen::Block<MatTCM,12,12> tensorC=MC.template block<12,12>(0,k*12);
    CONTRACTRC
  }
#undef CONTRACT
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABCCTD(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const
{
  Mat3XT MRRA,MRtA,MtRA,MttA,MRRC,MRtC,MtRC,MttC,GB;
  MRRA.setZero(3,body.nrJ()*3);
  MRtA.setZero(3,body.nrJ()*3);
  MtRA.setZero(3,body.nrJ()*3);
  MttA.setZero(3,body.nrJ()*3);
  MRRC.setZero(3,body.nrJ()*3);
  MRtC.setZero(3,body.nrJ()*3);
  MtRC.setZero(3,body.nrJ()*3);
  MttC.setZero(3,body.nrJ()*3);
  GB.setZero(3,body.nrJ()*4);
  toolABCCTDContactAll(body,mapM(MRRA),mapM(MRtA),mapM(MtRA),mapM(MttA),MA,mapM(MRRC),mapM(MRtC),mapM(MtRC),mapM(MttC),MC,mapM(GB));
  toolABCCTD(body,mapM(MRRA),mapM(MRtA),mapM(MtRA),mapM(MttA),mapM(MRRC),mapM(MRtC),mapM(MtRC),mapM(MttC),G,mapM(GB),[&](sizeType row,sizeType col,T val) {
    ABCCTD(row,col)+=val;
    ABCCTD(col,row)+=val;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABCCTDZ(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const
{
  ABCCTD.setZero();
  toolABCCTD(body,MA,MC,G,ABCCTD);
}
//this is a brute-force interface for testing
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABCCTDBF(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const
{
  toolCBF(body,*this,MC,ABCCTD);
  ABCCTD=(ABCCTD+ABCCTD.transpose()).eval();
  toolA(body,*this,MA,ABCCTD);
  toolBBF(body,MC,ABCCTD);
  toolDBF(body,G,ABCCTD);
#undef CONTRACT
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::toolABCCTDBFZ(const ArticulatedBody& body,MatTCM MA,MatTCM MC,Mat3XTM G,MatTM ABCCTD) const
{
  ABCCTD.setZero();
  toolABCCTDBF(body,MA,MC,G,ABCCTD);
}
//-------------------------------------------------------------Jacobian RC
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR) const
{
  for(sizeType j=JID; j>=0; j=body.joint(j)._parent) {
    const Joint& J=body.joint(j);
    for(Vec2i R=J.RBegEnd(); R[0]!=R[1]; R[0]++)
      JR(R[0],_RDTM.col(R[0]));
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRMSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Mat3T&)> JR) const
{
  JRSparse(body,JID,[&](sizeType c,const Vec3T& R) {
    JR(c,cross<T>(R)*ROTI(_TM,JID));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JCSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JC) const
{
  for(sizeType j=JID; j>=0; j=body.joint(j)._parent) {
    const Joint& J=body.joint(j);
    for(Vec2i C=J.CBegEnd(); C[0]!=C[1]; C[0]++)
      JC(C[0],_RDTM.col(C[0]));
    Vec3T DC=CTRI(_TM,JID)-CTRI(_TM,j);
    for(Vec2i R=J.RBegEnd(); R[0]!=R[1]; R[0]++)
      JC(R[0],_RDTM.col(R[0]).cross(DC));
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRCSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR,std::function<void(sizeType,const Vec3T&)> JC) const
{
  for(sizeType j=JID; j>=0; j=body.joint(j)._parent) {
    const Joint& J=body.joint(j);
    for(Vec2i C=J.CBegEnd(); C[0]!=C[1]; C[0]++)
      JC(C[0],_RDTM.col(C[0]));
    Vec3T DC=CTRI(_TM,JID)-CTRI(_TM,j);
    for(Vec2i R=J.RBegEnd(); R[0]!=R[1]; R[0]++) {
      JR(R[0],_RDTM.col(R[0]));
      JC(R[0],_RDTM.col(R[0]).cross(DC));
    }
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Mat3T& dR) const
{
  dR.setZero();
  JRMSparse(body,JID,[&](sizeType c,const Mat3T& R) {
    dR+=dx[c]*R;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JCSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Vec3T& dC) const
{
  dC.setZero();
  JCSparse(body,JID,[&](sizeType c,const Vec3T& v) {
    dC+=dx[c]*v;
  });
}
template <typename T>
typename PBDArticulatedGradientInfoMap<T>::Mat3T PBDArticulatedGradientInfoMap<T>::JRSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const
{
  Mat3T dR;
  JRSparse(body,JID,dx,dR);
  return dR;
}
template <typename T>
typename PBDArticulatedGradientInfoMap<T>::Vec3T PBDArticulatedGradientInfoMap<T>::JCSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const
{
  Vec3T dC;
  JCSparse(body,JID,dx,dC);
  return dC;
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dR) const
{
  dR.setZero();
  JRSparse(body,JID,[&](sizeType c,const Vec3T& v) {
    dR.col(c)+=v;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JCSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dC) const
{
  dC.setZero();
  JCSparse(body,JID,[&](sizeType c,const Vec3T& v) {
    dC.col(c)+=v;
  });
}
//-------------------------------------------------------------Jacobian RCLambda
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR) const
{
  Vec3T CP;
  for(sizeType j=JID; j>=0; j=body.joint(j)._parent) {
    const Joint& J=body.joint(j);
    CP=J._parent>=0?DWDLI(_DTLambdaM,J._parent):Vec3T(Vec3T::Zero());
    for(Vec2i R=J.RBegEnd(); R[0]!=R[1]; R[0]++)
      JR(R[0],CP.cross(_RDTM.col(R[0]))+_RDTILambdaM.col(R[0]));
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRMILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Mat3T&)> JR) const
{
  JRILambdaSparse(body,JID,[&](sizeType c,const Vec3T& R) {
    JR(c,(cross<T>(_RDTM.col(c))*cross<T>(DWDLI(_DTLambdaM,JID))+cross<T>(R))*ROTI(_TM,JID));
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JCILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JC) const
{
  Vec3T CP,DDTDL,DC,WIL;
  for(sizeType j=JID; j>=0; j=body.joint(j)._parent) {
    const Joint& J=body.joint(j);
    CP=J._parent>=0?DWDLI(_DTLambdaM,J._parent):Vec3T(Vec3T::Zero());
    DDTDL=DTDLI(_DTLambdaM,JID)-DTDLI(_DTLambdaM,j);
    DC=CTRI(_TM,JID)-CTRI(_TM,j);
    for(Vec2i R=J.RBegEnd(); R[0]!=R[1]; R[0]++) {
      WIL=CP.cross(_RDTM.col(R[0]))+_RDTILambdaM.col(R[0]);
      JC(R[0],_RDTM.col(R[0]).cross(DDTDL)+WIL.cross(DC));
    }
    for(Vec2i C=J.CBegEnd(); C[0]!=C[1]; C[0]++)
      JC(C[0],CP.cross(_RDTM.col(C[0])));
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRCILambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,const Vec3T&)> JR,std::function<void(sizeType,const Vec3T&)> JC) const
{
  Vec3T CP,DDTDL,DC,WIL;
  for(sizeType j=JID; j>=0; j=body.joint(j)._parent) {
    const Joint& J=body.joint(j);
    CP=J._parent>=0?DWDLI(_DTLambdaM,J._parent):Vec3T(Vec3T::Zero());
    DDTDL=DTDLI(_DTLambdaM,JID)-DTDLI(_DTLambdaM,j);
    DC=CTRI(_TM,JID)-CTRI(_TM,j);
    for(Vec2i R=J.RBegEnd(); R[0]!=R[1]; R[0]++) {
      WIL=CP.cross(_RDTM.col(R[0]))+_RDTILambdaM.col(R[0]);
      JR(R[0],WIL);
      JC(R[0],_RDTM.col(R[0]).cross(DDTDL)+WIL.cross(DC));
    }
    for(Vec2i C=J.CBegEnd(); C[0]!=C[1]; C[0]++)
      JC(C[0],CP.cross(_RDTM.col(C[0])));
  }
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Mat3T& dR) const
{
  dR.setZero();
  JRMILambdaSparse(body,JID,[&](sizeType c,const Mat3T& R) {
    dR+=dx[c]*R;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JCILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx,Vec3T& dC) const
{
  dC.setZero();
  JCILambdaSparse(body,JID,[&](sizeType c,const Vec3T& v) {
    dC+=dx[c]*v;
  });
}
template <typename T>
typename PBDArticulatedGradientInfoMap<T>::Mat3T PBDArticulatedGradientInfoMap<T>::JRILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const
{
  Mat3T dR;
  JRILambdaSparse(body,JID,dx,dR);
  return dR;
}
template <typename T>
typename PBDArticulatedGradientInfoMap<T>::Vec3T PBDArticulatedGradientInfoMap<T>::JCILambdaSparse(const ArticulatedBody& body,sizeType JID,VecM dx) const
{
  Vec3T dC;
  JCILambdaSparse(body,JID,dx,dC);
  return dC;
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JRILambdaSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dR) const
{
  dR.setZero();
  JRILambdaSparse(body,JID,[&](sizeType c,const Vec3T& v) {
    dR.col(c)+=v;
  });
}
template <typename T>
void PBDArticulatedGradientInfoMap<T>::JCILambdaSparse(const ArticulatedBody& body,sizeType JID,Mat3XTM dC) const
{
  dC.setZero();
  JCILambdaSparse(body,JID,[&](sizeType c,const Vec3T& v) {
    dC.col(c)+=v;
  });
}
//-------------------------------------------------------------Hessian RK_1K
template <typename T>
void PBDArticulatedGradientInfoMap<T>::HRK_1KSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,sizeType,const Mat3T&)> HR) const
{
  Mat3T DRIJ;
  const Joint& J=body.joint(JID);
  for(Vec2i RI=J.RBegEnd(); RI[0]<RI[1]; RI[0]++)
    for(Vec2i RJ=J.RBegEnd(); RJ[0]<RJ[1]; RJ[0]++) {
      DRIJ=cross<T>(_DTM.col(RI[0]))*cross<T>(_DTM.col(RJ[0]))+cross<T>(JointFunc<T>::DDT(J,mapM2CM(_DDTM),RJ[0],RI[0]));
      HR(RI[0],RJ[0],DRIJ*ROTI(_TK_1KM,JID));
    }
}
//-------------------------------------------------------------Hessian RK_1KLambda
template <typename T>
void PBDArticulatedGradientInfoMap<T>::HRK_1KLambdaSparse(const ArticulatedBody& body,sizeType JID,std::function<void(sizeType,sizeType,const Mat3T&)> HR) const
{
  Mat3T DRIJ;
  const Joint& J=body.joint(JID);
  for(Vec2i RI=J.RBegEnd(); RI[0]<RI[1]; RI[0]++) {
    for(Vec2i RJ=J.RBegEnd(); RJ[0]<RJ[1]; RJ[0]++) {
      DRIJ=(cross<T>(_DTM.col(RI[0]))*cross<T>(_DTM.col(RJ[0]))+cross<T>(JointFunc<T>::DDT(J,mapM2CM(_DDTM),RJ[0],RI[0])))*cross<T>(DWDLI(_DTK_1KLambdaM,JID));
      DRIJ+=cross<T>(_DTILambdaM.col(RI[0]))*cross<T>(_DTM.col(RJ[0]))+cross<T>(_DTM.col(RI[0]))*cross<T>(_DTILambdaM.col(RJ[0]))+cross<T>(JointFunc<T>::DDT(J,mapM2CM(_DTIILambdaM),RJ[0],RI[0]));
      HR(RI[0],RJ[0],DRIJ*ROTI(_TK_1KM,JID));
    }
  }
}
//ArticulatedGradientInfo
template <typename T>
PBDArticulatedGradientInfo<T>::PBDArticulatedGradientInfo() {}
template <typename T>
PBDArticulatedGradientInfo<T>::PBDArticulatedGradientInfo(const PBDArticulatedGradientInfo& other)
{
  operator=(other);
}
template <typename T>
PBDArticulatedGradientInfo<T>::PBDArticulatedGradientInfo(const ArticulatedBody& body,const Vec& x)
{
  reset(body,x);
}
template <typename T>
PBDArticulatedGradientInfo<T>& PBDArticulatedGradientInfo<T>::operator=(const PBDArticulatedGradientInfo& other)
{
  _x=other._x;
  _T=other._T;
  _TK_1K=other._TK_1K;
  _DT=other._DT;
  _RDT=other._RDT;
  _DDT=other._DDT;
  _DTLambda=other._DTLambda;
  _DTK_1KLambda=other._DTK_1KLambda;
  _DTILambda=other._DTILambda;
  _RDTILambda=other._RDTILambda;
  _DTIILambda=other._DTIILambda;
  _RDTIILambda=other._RDTIILambda;
  resetPtr();
  return *this;
}
template <typename T>
void PBDArticulatedGradientInfo<T>::resetLambda(const ArticulatedBody& body,const Vec& lambda)
{
  PBDArticulatedGradientInfoMap<T>::resetLambda(body,mapV(lambda));
}
template <typename T>
void PBDArticulatedGradientInfo<T>::reset(const ArticulatedBody& body,const Vec& x)
{
  sizeType nrJ=body.nrJ();
  sizeType nrDOF=body.nrDOF();
  sizeType nrDDT=body.nrDDT();
  _x=x;
  _T.resize(3,nrJ*4);
  _TK_1K.resize(3,nrJ*4);
  _DT.resize(3,nrDOF);
  _RDT.resize(3,nrDOF);
  _DDT.resize(3,nrDDT);
  reorthogonalize(body);
  _DTLambda.resize(3,nrJ*2);
  _DTK_1KLambda.resize(3,nrJ*2);
  _DTILambda.resize(3,nrDOF);
  _RDTILambda.resize(3,nrDOF);
  _DTIILambda.resize(3,nrDDT);
  _RDTIILambda.resize(3,nrDDT);
  resetPtr();
  PBDArticulatedGradientInfoMap<T>::reset(body,mapV(x));
}
template <typename T>
void PBDArticulatedGradientInfo<T>::reorthogonalize(const ArticulatedBody& body)
{
  sizeType nrJ=body.nrJ();
  //this re-orthogonalize JTrans, which is essential for MPFR precision
  //this is because JTrans was originally calculated using low-precision
  if(_JTrans.rows()!=3 || _JTrans.cols()!=nrJ*4) {
    _JTrans.resize(3,nrJ*4);
    for(sizeType j=0; j<nrJ; j++) {
      Mat3X4T JTrans=body.joint(j)._trans.template cast<T>();
      JTrans.col(0)/=std::sqrt(JTrans.col(0).squaredNorm());
      JTrans.col(1)-=JTrans.col(1).dot(JTrans.col(0))*JTrans.col(0);
      JTrans.col(1)/=std::sqrt(JTrans.col(1).squaredNorm());
      JTrans.col(2)=JTrans.col(0).cross(JTrans.col(1));
      TRANSI(_JTrans,j)=JTrans;
    }
  }
}
//debug
template <typename T>
void PBDArticulatedGradientInfo<T>::debug(const ArticulatedBody& body)
{
  sizeType nrJ=body.nrJ();
  sizeType nrDOF=body.nrDOF();

  Vec x,x2,dx,lambda,g,g2,gBF;
  x.resize(nrDOF);
  x2.resize(nrDOF);
  dx.resize(nrDOF);
  lambda.resize(nrDOF);
  g.resize(nrDOF);
  g2.resize(nrDOF);
  gBF.resize(nrDOF);

  MatT ABF,ABF2,A,B,B2,AB,AB2,AB2BF,AC,ACBF,DBF,ABCCTD,ABCCTDBF,M,M2;
  ABF.resize(nrDOF,nrDOF);
  ABF2.resize(nrDOF,nrDOF);
  A.resize(nrDOF,nrDOF);
  B.resize(nrDOF,nrDOF);
  B2.resize(nrDOF,nrDOF);
  AB.resize(nrDOF,nrDOF);
  AB2.resize(nrDOF,nrDOF);
  AB2BF.resize(nrDOF,nrDOF);
  AC.resize(nrDOF,nrDOF);
  ACBF.resize(nrDOF,nrDOF);
  DBF.resize(nrDOF,nrDOF);
  ABCCTD.resize(nrDOF,nrDOF);
  ABCCTDBF.resize(nrDOF,nrDOF);
  M.resize(12,nrJ*12);
  M2.resize(12,nrJ*12);

  Mat3XT JRCF0,JRCF,HRss[3][3],HRss2[3][3];
  JRCF0.resize(3,nrJ*4);
  JRCF.resize(3,nrJ*4);

  PBDArticulatedGradientInfo G,G2,G3,G4;
  DEFINE_NUMERIC_DELTA_T(T)
  for(sizeType JID=0; JID<nrJ; JID++) {
    INFO("-------------------------------------------------------------DebugPBDArticulatedGradientInfo")
    //calculate info
    x.setRandom();
    x2.setRandom();
    dx.setRandom();
    lambda.setRandom();
    JRCF0.setRandom();
    G.reset(body,x);
    G2.reset(body,x+dx*DELTA);
    G3.reset(body,x2);
    G4.reset(body,x+lambda*DELTA);
    //debug JR,JC
    DEBUG_GRADIENT("dR",std::sqrt(G.JRSparse(body,JID,mapV(dx)).squaredNorm()),std::sqrt((G.JRSparse(body,JID,mapV(dx))-(ROTI(G2._TM,JID)-ROTI(G._TM,JID))/DELTA).squaredNorm()))
    DEBUG_GRADIENT("dC",std::sqrt(G.JCSparse(body,JID,mapV(dx)).squaredNorm()),std::sqrt((G.JCSparse(body,JID,mapV(dx))-(CTRI(G2._TM,JID)-CTRI(G._TM,JID))/DELTA).squaredNorm()))
    //debug JRILambda,JCILambda
    G.resetLambda(body,lambda);
    DEBUG_GRADIENT("dRILambda",std::sqrt(G.JRILambdaSparse(body,JID,mapV(dx)).squaredNorm()),std::sqrt((G.JRILambdaSparse(body,JID,mapV(dx))-(G4.JRSparse(body,JID,mapV(dx))-G.JRSparse(body,JID,mapV(dx)))/DELTA).squaredNorm()))
    DEBUG_GRADIENT("dCILambda",std::sqrt(G.JCILambdaSparse(body,JID,mapV(dx)).squaredNorm()),std::sqrt((G.JCILambdaSparse(body,JID,mapV(dx))-(G4.JCSparse(body,JID,mapV(dx))-G.JCSparse(body,JID,mapV(dx)))/DELTA).squaredNorm()))
    //debug HR
    G.HRK_1KSparse(body,JID,[&](sizeType row,sizeType col,const Mat3T& HR) {
      PBDArticulatedGradientInfo G5(body,x+Vec::Unit(x.size(),col)*DELTA);
      Mat3T DR=cross<T>(G._DTM.col(row))*ROTI(G._TK_1KM,JID);
      Mat3T DR2=cross<T>(G5._DTM.col(row))*ROTI(G5._TK_1KM,JID);
      DEBUG_GRADIENT("HR",std::sqrt(HR.squaredNorm()),std::sqrt((HR-(DR2-DR)/DELTA).squaredNorm()))
      row-=body.joint(JID)._offDOF;
      col-=body.joint(JID)._offDOF;
      HRss[row][col]=HR;
    });
    G4.HRK_1KSparse(body,JID,[&](sizeType row,sizeType col,const Mat3T& HR) {
      row-=body.joint(JID)._offDOF;
      col-=body.joint(JID)._offDOF;
      HRss2[row][col]=HR;
    });
    G.HRK_1KLambdaSparse(body,JID,[&](sizeType row,sizeType col,const Mat3T HRLambda) {
      row-=body.joint(JID)._offDOF;
      col-=body.joint(JID)._offDOF;
      DEBUG_GRADIENT("HRLambda",std::sqrt(HRLambda.squaredNorm()),std::sqrt((HRLambda-(HRss2[row][col]-HRss[row][col])/DELTA).squaredNorm()))
    });
    //debug DTG
    T E=G.TG(body,mapCM(JRCF=JRCF0));
    T E2=G2.TG(body,mapCM(JRCF=JRCF0));
    G.DTGZ(body,mapM(JRCF=JRCF0),mapV(g));
    G2.DTGZ(body,mapM(JRCF=JRCF0),mapV(g2));
    DEBUG_GRADIENT("DTG",g.dot(dx),(g.dot(dx)-(E2-E)/DELTA))
    G.DTGBFZ(body,mapM(JRCF=JRCF0),mapV(gBF));
    DEBUG_GRADIENT("DTGBF",std::sqrt(g.squaredNorm()),std::sqrt((g-gBF).squaredNorm()))
    //debug toolA
    M.setRandom();
    G.toolAZ(body,G3,mapCM(M),mapM(A));
    G.toolABFZ(body,G3,mapCM(M),mapM(ABF));
    G.toolABF2Z(body,G3,mapCM(M),mapM(ABF2));
    DEBUG_GRADIENT("toolABF",std::sqrt(ABF.squaredNorm()),std::sqrt((ABF-A).squaredNorm()))
    DEBUG_GRADIENT("toolABF2",std::sqrt(ABF2.squaredNorm()),std::sqrt((ABF2-A).squaredNorm()))
    //debug toolB
    G.toolBZ(body,mapM(JRCF=JRCF0),mapM(B));
    DEBUG_GRADIENT("toolB",std::sqrt((B*dx).squaredNorm()),std::sqrt((B*dx-(g2-g)/DELTA).squaredNorm()))
    //debug toolAB
    G.toolAZ(body,G,mapCM(M),mapM(AB));
    G.toolBZ(body,mapM(JRCF=JRCF0),mapM(AB2));
    AB+=AB2;
    G.toolABZ(body,mapCM(M),mapM(JRCF=JRCF0),mapM(AB2));
    DEBUG_GRADIENT("toolAB",std::sqrt(AB.squaredNorm()),std::sqrt((AB-AB2).squaredNorm()))
    G.toolABBFZ(body,mapCM(M),mapM(JRCF=JRCF0),mapM(AB2BF));
    DEBUG_GRADIENT("toolABBF",std::sqrt(AB2.squaredNorm()),std::sqrt((AB2-AB2BF).squaredNorm()))
    //debug toolAC
    M2.setRandom();
    G.toolACZ(body,G3,mapCM(M2),mapCM(M),mapM(AC));
    G.toolCBFZ(body,G3,mapCM(M),mapM(ACBF));
    G.toolABF(body,G3,mapCM(M2),mapM(ACBF));
    DEBUG_GRADIENT("toolACBF",std::sqrt(ACBF.squaredNorm()),std::sqrt((ACBF-AC).squaredNorm()))
    //debug toolDBF
    G.toolDBFZ(body,mapM(JRCF=JRCF0),mapM(DBF));
    G4.toolBZ(body,mapM(JRCF=JRCF0),mapM(B2));
    DEBUG_GRADIENT("toolDBF",std::sqrt(DBF.squaredNorm()),std::sqrt((DBF-(B2-B)/DELTA).squaredNorm()))
    //debug toolABCCTD
    for(sizeType i=0; i<M.cols(); i+=12)
      M.template block<12,12>(0,i)=(M.template block<12,12>(0,i)+M.template block<12,12>(0,i).transpose()).eval();
    G.toolABCCTDZ(body,mapCM(M),mapCM(M2),mapM(JRCF=JRCF0),mapM(ABCCTD));
    G.toolABCCTDBFZ(body,mapCM(M),mapCM(M2),mapM(JRCF=JRCF0),mapM(ABCCTDBF));
    DEBUG_GRADIENT("toolABCCTD",std::sqrt(ABCCTD.squaredNorm()),std::sqrt((ABCCTD-ABCCTDBF).squaredNorm()))
  }
}
template <typename T>
void PBDArticulatedGradientInfo<T>::resetPtr()
{
  new (&_xM)           VecM(mapV(_x));
  new (&_TM)           Mat3XTM(mapM(_T));
  new (&_TK_1KM)       Mat3XTM(mapM(_TK_1K));
  new (&_DTM)          Mat3XTM(mapM(_DT));
  new (&_RDTM)         Mat3XTM(mapM(_RDT));
  new (&_DDTM)         Mat3XTM(mapM(_DDT));
  new (&_JTransM)      Mat3XTM(mapM(_JTrans));
  new (&_DTLambdaM)    Mat3XTM(mapM(_DTLambda));
  new (&_DTK_1KLambdaM)Mat3XTM(mapM(_DTK_1KLambda));
  new (&_DTILambdaM)   Mat3XTM(mapM(_DTILambda));
  new (&_RDTILambdaM)  Mat3XTM(mapM(_RDTILambda));
  new (&_DTIILambdaM)  Mat3XTM(mapM(_DTIILambda));
  new (&_RDTIILambdaM) Mat3XTM(mapM(_RDTIILambda));
}
//instance
PRJ_BEGIN
template struct PBDArticulatedGradientInfoMap<double>;
template struct PBDArticulatedGradientInfo<double>;
#ifdef ALL_TYPES
template struct PBDArticulatedGradientInfoMap<__float128>;
template struct PBDArticulatedGradientInfo<__float128>;
template struct PBDArticulatedGradientInfoMap<mpfr::mpreal>;
template struct PBDArticulatedGradientInfo<mpfr::mpreal>;
#endif
PRJ_END
