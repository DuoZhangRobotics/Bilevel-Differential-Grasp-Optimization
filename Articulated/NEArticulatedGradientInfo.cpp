#include "NEArticulatedGradientInfo.h"
#include "PBDArticulatedGradientInfo.h"
#include "TensorContractPragma.h"
#include "JointFunc.h"
#include <Utils/CrossSpatialUtil.h>
#include <Utils/DebugGradient.h>
#include <Utils/Scalar.h>

USE_PRJ_NAMESPACE

//CPU function API
template <typename T>
NEArticulatedGradientInfoMap<T>::NEArticulatedGradientInfoMap()
  :_sparsityM(NULL,0),
   _nonZeroM(NULL,0),
   _qM(NULL,0),
   _dqM(NULL,0),
   _ddqM(NULL,0),
   _tauM(NULL,0),

   _invTPM(NULL,6,0,Eigen::OuterStride<>(6)),
   _invTM(NULL,6,0,Eigen::OuterStride<>(6)),
   _SM(NULL,6,0,Eigen::OuterStride<>(6)),
   _UM(NULL,6,0,Eigen::OuterStride<>(6)),
   _vM(NULL,6,0,Eigen::OuterStride<>(6)),
   _vPM(NULL,6,0,Eigen::OuterStride<>(6)),
   _aM(NULL,6,0,Eigen::OuterStride<>(6)),
   _fM(NULL,6,0,Eigen::OuterStride<>(6)),
   _ICM(NULL,6,0,Eigen::OuterStride<>(6)),
   _JTransM(NULL,6,0,Eigen::OuterStride<>(6)),
   _IM(NULL,6,0,Eigen::OuterStride<>(6)),
   _DvJDqM(NULL,6,0,Eigen::OuterStride<>(6)),
   _DdvJDqM(NULL,6,0,Eigen::OuterStride<>(6)),
   _DdvJDdqM(NULL,6,0,Eigen::OuterStride<>(6)),

   _DtauDqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _DtauDdqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _DvDqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _DvDdqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _DaDqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _DaDdqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _DfDqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _DfDdqM(NULL,0,0,Eigen::OuterStride<>(0)),
   _HM(NULL,0,0,Eigen::OuterStride<>(0)),
   _invHM(NULL,0,0,Eigen::OuterStride<>(0)) {}
template <typename T>
NEArticulatedGradientInfoMap<T>::NEArticulatedGradientInfoMap(const NEArticulatedGradientInfoMap& other)
  :_sparsityM(NULL,0),
   _nonZeroM(NULL,0),
   _qM((VecM&)other._qM),
   _dqM((VecM&)other._dqM),
   _ddqM((VecM&)other._ddqM),
   _tauM((VecM&)other._tauM),

   _invTPM((Mat6XTM&)other._invTPM),
   _invTM((Mat6XTM&)other._invTM),
   _SM((Mat6XTM&)other._SM),
   _UM((Mat6XTM&)other._UM),
   _vM((Mat6XTM&)other._vM),
   _vPM((Mat6XTM&)other._vPM),
   _aM((Mat6XTM&)other._aM),
   _fM((Mat6XTM&)other._fM),
   _ICM((Mat6XTM&)other._IM),
   _JTransM((Mat6XTM&)other._JTransM),
   _IM((Mat6XTM&)other._IM),
   _DvJDqM((Mat6XTM&)other._DvJDqM),
   _DdvJDqM((Mat6XTM&)other._DdvJDqM),
   _DdvJDdqM((Mat6XTM&)other._DdvJDdqM),

   _DtauDqM((MatTM&)other._DtauDqM),
   _DtauDdqM((MatTM&)other._DtauDdqM),
   _DvDqM((MatTM&)other._DvDqM),
   _DvDdqM((MatTM&)other._DvDdqM),
   _DaDqM((MatTM&)other._DaDqM),
   _DaDdqM((MatTM&)other._DaDdqM),
   _DfDqM((MatTM&)other._DfDqM),
   _DfDdqM((MatTM&)other._DfDdqM),
   _HM((MatTM&)other._HM),
   _invHM((MatTM&)other._invHM) {}
template <typename T>
void NEArticulatedGradientInfoMap<T>::reset(const ArticulatedBody& body,VecCM qMap,VecCM dqMap,VecCM ddqMap)
{
  reset(body,qMap,dqMap);
  _ddqM=ddqMap;
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::reset(const ArticulatedBody& body,VecCM qMap,VecCM dqMap)
{
  reset(body,qMap);
  _dqM=dqMap;
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::reset(const ArticulatedBody& body,VecCM qMap)
{
  _qM=qMap;
  Mat6T TJ;
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    JointFunc<T>::JCALC(J,mapV2CV(_qM),mapCV<Vec>((Vec*)NULL),mapCV<Vec>((Vec*)NULL),mapM(TJ),_SM,mapV<Vec6T>((Vec6T*)NULL),mapV<Vec6T>((Vec6T*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapM<Mat6XT>((Mat6XT*)NULL));
    TRANSI6(_invTPM,i)=spatialInv<T>(TRANSI6(_JTransM,i)*TJ);
    if(J._parent>=0)
      TRANSI6(_invTM,i)=TRANSI6(_invTPM,i)*TRANSI6(_invTM,J._parent);
    else TRANSI6(_invTM,i)=TRANSI6(_invTPM,i);
  }
  _updateH=_updateInvH=false;
}
//-------------------------------------------------------------RNEA
template <typename T>
void NEArticulatedGradientInfoMap<T>::RNEA(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,bool useDDQ)
{
  Vec6T vJ,dvJ;
  VecCM ddqM=useDDQ?mapV2CV<Vec>(_ddqM):mapCV<Vec>((Vec*)NULL);
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    JointFunc<T>::JCALC(J,mapV2CV(_qM),mapV2CV(_dqM),ddqM,mapM<Mat6T>((Mat6T*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapV(vJ),mapV(dvJ),mapM<Mat6XT>((Mat6XT*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapM<Mat6XT>((Mat6XT*)NULL));
    if(J._parent>=0) {
      _vM.col(i)=TRANSI6(_invTPM,i)*_vM.col(J._parent);
      _aM.col(i)=TRANSI6(_invTPM,i)*_aM.col(J._parent);
    } else {
      _vM.col(i).setZero();
      if(a0.data())
        _aM.col(i)=TRANSI6(_invTPM,i)*a0;
      else _aM.col(i).setZero();
    }
    _vM.col(i)+=vJ;
    _aM.col(i)+=dvJ+spatialCross<T>(_vM.col(i),vJ);
    _fM.col(i)=TRANSI6(_IM,i)*_aM.col(i);
    _fM.col(i)+=spatialCrossStar<T>(_vM.col(i),TRANSI6(_IM,i)*_vM.col(i));
    if(fx.data())
      _fM.col(i)-=spatialXStar<T>(TRANSI6(_invTM,i))*fx.col(i);
  }
  for(sizeType i=body.nrJ()-1; i>=0; i--) {
    const Joint& J=body.joint(i);
    _tauM.segment(J._offDOF,J.nrDOF())=_SM.block(0,J._offDOF,6,J.nrDOF()).transpose()*_fM.col(i);
    if(J._parent>=0)
      _fM.col(J._parent)+=TRANSI6(_invTPM,i).transpose()*_fM.col(i);
  }
}
//-------------------------------------------------------------RNEADerivatives
template <typename T>
void NEArticulatedGradientInfoMap<T>::DvDdq(const ArticulatedBody& body,sizeType i)
{
  const Joint& J=body.joint(i);
  sizeType nrDOF=J.nrDOF();
  _DvDdqM.block(i*6,J._offDOF,6,nrDOF)+=_SM.block(0,J._offDOF,6,nrDOF);
  for(sizeType j=body.joint(i)._parent,p=j; j>=0; j=body.joint(j)._parent) {
    const Joint& JP=body.joint(j);
    sizeType nrDOFP=JP.nrDOF();
    _DvDdqM.block(i*6,JP._offDOF,6,nrDOFP)+=TRANSI6(_invTPM,i)*_DvDdqM.block(p*6,JP._offDOF,6,nrDOFP);
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DvDq(const ArticulatedBody& body,sizeType i)
{
  const Joint& J=body.joint(i);
  sizeType nrDOF=J.nrDOF();
  for(sizeType d=0; d<nrDOF; d++)
    _DvDqM.block(i*6,J._offDOF+d,6,1)+=_DvJDqM.col(J._offDOF+d)-spatialCross<T>(_SM.col(J._offDOF+d),_vM.col(i));
  for(sizeType j=body.joint(i)._parent,p=j; j>=0; j=body.joint(j)._parent) {
    const Joint& JP=body.joint(j);
    sizeType nrDOFP=JP.nrDOF();
    _DvDqM.block(i*6,JP._offDOF,6,nrDOFP)+=TRANSI6(_invTPM,i)*_DvDqM.block(p*6,JP._offDOF,6,nrDOFP);
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DaDdq(const ArticulatedBody& body,sizeType i,Vec6TCM vJ)
{
  const Joint& J=body.joint(i);
  sizeType nrDOF=J.nrDOF();
  for(sizeType d=0; d<nrDOF; d++)
    _DaDdqM.block(i*6,J._offDOF+d,6,1)+=spatialCross<T>(_vM.col(i),_SM.col(J._offDOF+d))+_DdvJDdqM.col(J._offDOF+d);
  for(sizeType j=body.joint(i)._parent,p=j; j>=0; j=body.joint(j)._parent) {
    const Joint& JP=body.joint(j);
    sizeType nrDOFP=JP.nrDOF();
    _DaDdqM.block(i*6,JP._offDOF,6,nrDOFP)+=TRANSI6(_invTPM,i)*_DaDdqM.block(p*6,JP._offDOF,6,nrDOFP);
  }
  for(sizeType j=i; j>=0; j=body.joint(j)._parent) {
    const Joint& JP=body.joint(j);
    sizeType nrDOFP=JP.nrDOF();
    for(sizeType d=0; d<nrDOFP; d++)
      _DaDdqM.block(i*6,JP._offDOF+d,6,1)+=spatialCross<T>(_DvDdqM.block(i*6,JP._offDOF+d,6,1),vJ);
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DaDq(const ArticulatedBody& body,sizeType i,Vec6TCM vJ)
{
  const Joint& J=body.joint(i);
  sizeType nrDOF=J.nrDOF();
  _DaDqM.block(i*6,J._offDOF,6,nrDOF)+=_DdvJDqM.block(0,J._offDOF,6,nrDOF);
  if(vJ.data()) //otherwise we are computing DHDq
    for(sizeType d=0; d<nrDOF; d++) {
      _DaDqM.block(i*6,J._offDOF+d,6,1)+=spatialCross<T>(vJ,spatialCross<T>(_vM.col(i),_SM.col(J._offDOF+d)));
      _DaDqM.block(i*6,J._offDOF+d,6,1)+=spatialCross<T>(_vM.col(i),_DvJDqM.col(J._offDOF+d));
    }
  for(sizeType d=0; d<nrDOF; d++)
    _DaDqM.block(i*6,J._offDOF+d,6,1)-=spatialCross<T>(_SM.col(J._offDOF+d),_aM.col(i));
  for(sizeType j=body.joint(i)._parent,p=j; j>=0; j=body.joint(j)._parent) {
    const Joint& JP=body.joint(j);
    sizeType nrDOFP=JP.nrDOF();
    _DaDqM.block(i*6,JP._offDOF,6,nrDOFP)+=TRANSI6(_invTPM,i)*_DaDqM.block(p*6,JP._offDOF,6,nrDOFP);
  }
  if(vJ.data()) //otherwise we are computing DHDq
    for(sizeType j=i; j>=0; j=body.joint(j)._parent) {
      const Joint& JP=body.joint(j);
      sizeType nrDOFP=JP.nrDOF();
      for(sizeType d=0; d<nrDOFP; d++)
        _DaDqM.block(i*6,JP._offDOF+d,6,1)+=spatialCross<T>(_DvDqM.block(i*6,JP._offDOF+d,6,1),vJ);
    }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::Dfx(const ArticulatedBody& body,sizeType i,Vec6TCM fx)
{
  Mat6T XStar=spatialXStar<T>(TRANSI6(_invTM,i));
  const Joint& J=body.joint(i);
  sizeType nrDOF=J.nrDOF();
  _vPM.block(0,J._offDOF,6,nrDOF)=XStar.transpose()*_SM.block(0,J._offDOF,6,nrDOF);
  for(sizeType j=i; j>=0; j=body.joint(j)._parent) {
    const Joint& JP=body.joint(j);
    sizeType nrDOFP=JP.nrDOF();
    for(sizeType d=0; d<nrDOFP; d++)
      _DfDqM.block(i*6,JP._offDOF+d,6,1)+=XStar*spatialCrossStar<T>(_vPM.col(JP._offDOF+d),fx);
  }
  _fM.col(i)-=XStar*fx;
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::Df(const ArticulatedBody& body,sizeType i,Vec6TCM vJ,Mat6XTCM IMCustom)
{
  if(!vJ.data()) {
    for(sizeType j=i; j>=0; j=body.joint(j)._parent) {
      const Joint& JP=body.joint(j);
      sizeType nrDOFP=JP.nrDOF();
      _DfDqM.block(i*6,JP._offDOF,6,nrDOFP)+=TRANSI6(IMCustom,i)*_DaDqM.block(i*6,JP._offDOF,6,nrDOFP);
    }
  } else {
    Vec6T IvR=TRANSI6(IMCustom,i)*_vM.col(i);
    Mat6T IvL=spatialCrossStar<T>(_vM.col(i))*TRANSI6(IMCustom,i);
    for(sizeType j=i; j>=0; j=body.joint(j)._parent) {
      const Joint& JP=body.joint(j);
      sizeType nrDOFP=JP.nrDOF();
      for(sizeType d=0; d<nrDOFP; d++) {
        _DfDqM.block(i*6,JP._offDOF+d,6,1)+=TRANSI6(IMCustom,i)*_DaDqM.block(i*6,JP._offDOF+d,6,1);
        _DfDqM.block(i*6,JP._offDOF+d,6,1)+=spatialCrossStar<T>(_DvDqM.block(i*6,JP._offDOF+d,6,1),IvR);
        _DfDqM.block(i*6,JP._offDOF+d,6,1)+=IvL*_DvDqM.block(i*6,JP._offDOF+d,6,1);

        _DfDdqM.block(i*6,JP._offDOF+d,6,1)+=TRANSI6(IMCustom,i)*_DaDdqM.block(i*6,JP._offDOF+d,6,1);
        _DfDdqM.block(i*6,JP._offDOF+d,6,1)+=spatialCrossStar<T>(_DvDdqM.block(i*6,JP._offDOF+d,6,1),IvR);
        _DfDdqM.block(i*6,JP._offDOF+d,6,1)+=IvL*_DvDdqM.block(i*6,JP._offDOF+d,6,1);
      }
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::Dtau(const ArticulatedBody& body,sizeType i,sizeType j)
{
  const Joint& J=body.joint(i);
  const Joint& JJ=body.joint(j);
  sizeType nrDOF=J.nrDOF();
  sizeType nrDOFJ=JJ.nrDOF();
  _DtauDqM.block(J._offDOF,JJ._offDOF,nrDOF,nrDOFJ)+=_SM.block(0,J._offDOF,6,nrDOF).transpose()*_DfDqM.block(i*6,JJ._offDOF,6,nrDOFJ);
  _DtauDdqM.block(J._offDOF,JJ._offDOF,nrDOF,nrDOFJ)+=_SM.block(0,J._offDOF,6,nrDOF).transpose()*_DfDdqM.block(i*6,JJ._offDOF,6,nrDOFJ);
  if(j==i) {
    JointFunc<T>::DSTDqf(JJ,mapV2CV(_qM),_DtauDqM,_fM.col(i));
    if(JJ._parent>=0)
      Dtau(body,i,JJ._parent);
    for(sizeType c:JJ._children)
      Dtau(body,i,c);
  } else if(j<i) {
    if(JJ._parent>=0)
      Dtau(body,i,JJ._parent);
  } else {
    for(sizeType c:JJ._children)
      Dtau(body,i,c);
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DHDq(const ArticulatedBody& body,sizeType i,sizeType j,MatTM dHdq) const
{
  const Joint& J=body.joint(i);
  const Joint& JJ=body.joint(j);
  sizeType nrDOF=J.nrDOF();
  sizeType nrDOFJ=JJ.nrDOF();
  dHdq.block(J._offDOF,JJ._offDOF,nrDOF,nrDOFJ)+=_SM.block(0,J._offDOF,6,nrDOF).transpose()*_DfDqM.block(i*6,JJ._offDOF,6,nrDOFJ);
  if(j==i) {
    JointFunc<T>::DSTDqf(JJ,mapV2CV(_qM),dHdq,_fM.col(i));
    if(JJ._parent>=0)
      DHDq(body,i,JJ._parent,dHdq);
    for(sizeType c:JJ._children)
      DHDq(body,i,c,dHdq);
  } else if(j<i) {
    if(JJ._parent>=0)
      DHDq(body,i,JJ._parent,dHdq);
  } else {
    for(sizeType c:JJ._children)
      DHDq(body,i,c,dHdq);
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DfP(const ArticulatedBody& body,sizeType i,sizeType j,bool dq)
{
  sizeType p=body.joint(i)._parent;
  const Joint& J=body.joint(j);
  sizeType nrDOFJ=J.nrDOF();
  _DfDqM.block(p*6,J._offDOF,6,nrDOFJ)+=TRANSI6(_invTPM,i).transpose()*_DfDqM.block(i*6,J._offDOF,6,nrDOFJ);
  if(dq)
    _DfDdqM.block(p*6,J._offDOF,6,nrDOFJ)+=TRANSI6(_invTPM,i).transpose()*_DfDdqM.block(i*6,J._offDOF,6,nrDOFJ);
  if(j==i) {
    for(sizeType d=0; d<nrDOFJ; d++)
      _DfDqM.block(p*6,J._offDOF+d,6,1)+=TRANSI6(_invTPM,i).transpose()*spatialCrossStar<T>(_SM.col(J._offDOF+d),_fM.col(i));
    if(J._parent>=0)
      DfP(body,i,J._parent,dq);
    for(sizeType c:J._children)
      DfP(body,i,c,dq);
  } else if(j<i) {
    if(J._parent>=0)
      DfP(body,i,J._parent,dq);
  } else {
    for(sizeType c:J._children)
      DfP(body,i,c,dq);
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::RNEADerivatives(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,bool useDDQ)
{
  Vec6T vJ,dvJ;
  _DvDqM.setZero();
  _DvDdqM.setZero();
  _DaDqM.setZero();
  _DaDdqM.setZero();
  _DfDqM.setZero();
  _DfDdqM.setZero();
  _DtauDqM.setZero();
  _DtauDdqM.setZero();

  _DvJDqM.setZero();
  _DdvJDqM.setZero();
  _DdvJDdqM.setZero();
  VecCM ddqM=useDDQ?mapV2CV<Vec>(_ddqM):mapCV<Vec>((Vec*)NULL);
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    JointFunc<T>::JCALC(J,mapV2CV(_qM),mapV2CV(_dqM),ddqM,mapM<Mat6T>((Mat6T*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapV(vJ),mapV(dvJ),_DvJDqM,_DdvJDqM,_DdvJDdqM);
    if(J._parent>=0) {
      _vM.col(i)=TRANSI6(_invTPM,i)*_vM.col(J._parent);
      _aM.col(i)=TRANSI6(_invTPM,i)*_aM.col(J._parent);
    } else {
      _vM.col(i).setZero();
      if(a0.data())
        _aM.col(i)=TRANSI6(_invTPM,i)*a0;
      else _aM.col(i).setZero();
    }
    _vM.col(i)+=vJ;
    _aM.col(i)+=dvJ+spatialCross<T>(_vM.col(i),vJ);
    _fM.col(i)=TRANSI6(_IM,i)*_aM.col(i);
    _fM.col(i)+=spatialCrossStar<T>(_vM.col(i),TRANSI6(_IM,i)*_vM.col(i));
    if(fx.data())
      Dfx(body,i,Vec6TCM(&(fx.coeffRef(0,i))));
    DvDq(body,i);
    DvDdq(body,i);
    DaDq(body,i,mapCV(vJ));
    DaDdq(body,i,mapCV(vJ));
    Df(body,i,mapCV(vJ),mapM2CM(_IM));
  }
  for(sizeType i=body.nrJ()-1; i>=0; i--) {
    const Joint& J=body.joint(i);
    _tauM.segment(J._offDOF,J.nrDOF())=_SM.block(0,J._offDOF,6,J.nrDOF()).transpose()*_fM.col(i);
    Dtau(body,i,i);
    if(J._parent>=0) {
      _fM.col(J._parent)+=TRANSI6(_invTPM,i).transpose()*_fM.col(i);
      DfP(body,i,i,true);
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::RNEADerivativesFD(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecM q1,VecM dq1,VecM ddq1,T delta,bool useDDQ)
{
  q1=_qM;
  dq1=_dqM;
  ddq1=_ddqM;
  for(sizeType i=0; i<q1.size(); i++) {
    reset(body,mapV2CV(q1),mapV2CV(dq1),mapV2CV(ddq1));
    //base
    RNEA(body,a0,fx,useDDQ);
    _DtauDqM.col(i)=_tauM;
    //perturbed
    _qM[i]+=delta;
    RNEA(body,a0,fx,useDDQ);
    _DtauDqM.col(i)=(_tauM-_DtauDqM.col(i))/delta;
  }
  for(sizeType i=0; i<q1.size(); i++) {
    reset(body,mapV2CV(q1),mapV2CV(dq1),mapV2CV(ddq1));
    //base
    RNEA(body,a0,fx,useDDQ);
    _DtauDdqM.col(i)=_tauM;
    //perturbed
    _dqM[i]+=delta;
    RNEA(body,a0,fx,useDDQ);
    _DtauDdqM.col(i)=(_tauM-_DtauDdqM.col(i))/delta;
  }
}
template <typename T>
typename NEArticulatedGradientInfoMap<T>::MatTCM NEArticulatedGradientInfoMap<T>::getDtauDddq() const
{
  return mapM2CM(_HM);
}
//-------------------------------------------------------------MDP
template <typename T>
void NEArticulatedGradientInfoMap<T>::Sc(const ArticulatedBody& body,ColiCM joints,MatTM sc) const
{
  sc.setZero();
  for(sizeType i=0; i<(sizeType)joints.size(); i++) {
    sizeType ip=joints[i];
    while(ip>=0) {
      const Joint& JIP=body.joint(ip);
      sizeType nrDOF=JIP.nrDOF();
      Mat6XTCM SMBlk(&(_SM.coeffRef(0,JIP._offDOF)),6,nrDOF,Eigen::OuterStride<>(_SM.outerStride()));
      sc.block(i*6,JIP._offDOF,6,nrDOF)=spatialInv<T>(TRANSI6(_invTM,ip))*SMBlk;
      ip=JIP._parent;
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DScDqT(const ArticulatedBody& body,ColiCM joints,MatTM DscDqTr,VecCM r)
{
  _vM.setZero();
  _nonZeroM.setZero();
  for(sizeType i=0; i<(sizeType)joints.size(); i++) {
    sizeType ip=joints[i];
    while(ip>=0) {
      const Joint& JIP=body.joint(ip);
      _vM.col(ip)+=r.template segment<6>(i*6);
      _nonZeroM[ip]=1;
      ip=JIP._parent;
    }
  }
  Mat6T TJ,TJP;
  Mat6X3T TJSc,TJScp;
  DscDqTr.setZero();
  for(sizeType i=0; i<body.nrJ(); i++)
    if(_nonZeroM[i]==1) {
      const Joint& J=body.joint(i);
      sizeType nrDOF=J.nrDOF();
      TJ=spatialInv<T>(TRANSI6(_invTM,i));
      Mat6XTM TJScBlk(TJSc.data(),6,nrDOF,Eigen::OuterStride<>(TJSc.outerStride()));
      TJScBlk=TJ*_SM.block(0,J._offDOF,6,nrDOF);
      JointFunc<T>::DSTDqf(J,mapV2CV(_qM),DscDqTr,TJ.transpose()*_vM.col(i));
      //row
      sizeType ip=i;
      while(ip>=0) {
        const Joint& JIP=body.joint(ip);
        sizeType nrDOFP=JIP.nrDOF();
        TJP=spatialInv<T>(TRANSI6(_invTM,ip));
        Mat6XTM TJScBlkp(TJScp.data(),6,nrDOFP,Eigen::OuterStride<>(TJScp.outerStride()));
        MatTM DscDqTrBlkp(&(DscDqTr.coeffRef(J._offDOF,JIP._offDOF)),nrDOF,nrDOFP,Eigen::OuterStride<>(DscDqTr.outerStride()));
        TJScBlkp=TJP*_SM.block(0,JIP._offDOF,6,nrDOFP);
        for(sizeType d=0; d<nrDOFP; d++)
          DscDqTrBlkp.col(d)-=TJScBlk.transpose()*spatialCrossStar<T>(TJScBlkp.col(d),_vM.col(i));
        ip=JIP._parent;
      }
    }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DScDq(const ArticulatedBody& body,ColiCM joints,MatTM DscDqr,VecCM r) const
{
  Mat6T TJ;
  Vec6T Scr;
  Mat6X3T DSDqf;
  DscDqr.setZero();
  for(sizeType i=0; i<(sizeType)joints.size(); i++) {
    sizeType ip=joints[i];
    Scr.setZero();
    while(ip>=0) {
      const Joint& JIP=body.joint(ip);
      sizeType nrDOF=JIP.nrDOF();
      Mat6XTCM SMBlk(&(_SM.coeffRef(0,JIP._offDOF)),6,nrDOF,Eigen::OuterStride<>(_SM.outerStride()));
      MatTM DscDqrBlk(&(DscDqr.coeffRef(i*6,JIP._offDOF)),6,nrDOF,Eigen::OuterStride<>(DscDqr.outerStride()));
      TJ=spatialInv<T>(TRANSI6(_invTM,ip));

      DSDqf.setZero();
      MatTM DSDqfM(DSDqf.data(),6,nrDOF,Eigen::OuterStride<>(DSDqf.outerStride()));
      JointFunc<T>::DSDqf(JIP,mapV2CV(_qM),DSDqfM,VecCM(&(r.coeffRef(JIP._offDOF)),nrDOF));
      DscDqrBlk=TJ*DSDqfM;

      Scr+=TJ*(SMBlk*r.segment(JIP._offDOF,nrDOF));
      for(sizeType d=0; d<nrDOF; d++)
        DscDqrBlk.col(d)+=spatialCross<T>(TJ*SMBlk.col(d),Scr);
      ip=JIP._parent;
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::ScInvH(const ArticulatedBody& body,ColiCM joints,MatTM sc,MatTM scInvH) const
{
  scInvH.setZero();
  Sc(body,joints,sc);
  for(sizeType i=0; i<(sizeType)joints.size(); i++) {
    sizeType ip=joints[i];
    while(ip>=0) {
      const Joint& JIP=body.joint(ip);
      sizeType nrDOF=JIP.nrDOF();
      Eigen::Block<const MatTM> invHMBlk=_invHM.block(JIP._offDOF,0,nrDOF,_invHM.cols());
      Eigen::Block<MatTM> scBlk=sc.block(i*6,JIP._offDOF,6,nrDOF);
      scInvH.block(i*6,0,6,_invHM.cols())+=scBlk*invHMBlk;
      ip=JIP._parent;
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::ScInvHScT(const ArticulatedBody& body,ColiCM joints,MatTM sc,MatTM scInvH,MatTM scInvHscT) const
{
  scInvHscT.setZero();
  ScInvH(body,joints,sc,scInvH);
  for(sizeType i=0; i<(sizeType)joints.size(); i++) {
    sizeType ip=joints[i];
    while(ip>=0) {
      const Joint& JIP=body.joint(ip);
      sizeType nrDOF=JIP.nrDOF();
      Eigen::Block<MatTM> scBlk=sc.block(i*6,JIP._offDOF,6,nrDOF);
      Eigen::Block<MatTM> scInvHBlk=scInvH.block(0,JIP._offDOF,scInvH.rows(),nrDOF);
      scInvHscT.block(i*6,0,6,scInvHscT.cols())+=scBlk*scInvHBlk.transpose();
      ip=JIP._parent;
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DTG(sizeType k,const ArticulatedBody& body,Mat3X4T GK,std::function<void(sizeType,T)> dtg) const
{
  sizeType ip=k;
  Mat3X4T T0=getTrans(k);
  Vec3T wCoef=invCrossMatTrace<T>(ROT(T0)*ROT(GK).transpose());
  while(ip>=0) {
    const Joint& JIP=body.joint(ip);
    sizeType nrDOF=JIP.nrDOF();
    Mat6XTCM SMBlk(&(_SM.coeffRef(0,JIP._offDOF)),6,nrDOF,Eigen::OuterStride<>(_SM.outerStride()));
    for(sizeType d=0; d<nrDOF; d++) {
      Vec6T TJSMBlk=spatialInv<T>(TRANSI6(_invTM,ip))*SMBlk.col(d);
      Vec3T w=TJSMBlk.template segment<3>(0);
      Vec3T v=TJSMBlk.template segment<3>(3)+w.cross(CTR(T0));
      dtg(JIP._offDOF+d,w.dot(wCoef)+v.dot(CTR(GK)));
    }
    ip=JIP._parent;
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DTG(sizeType k,const ArticulatedBody& body,Mat3X4T GK,VecM dtg) const
{
  DTG(k,body,GK,[&](sizeType off,T val) {
    dtg[off]+=val;
  });
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::DHDq(const ArticulatedBody& body,MatTM dHdq,VecCM r,Mat6XTCM IMCustom)
{
  if(!IMCustom.data())
    new (&IMCustom) Mat6XTCM(mapM2CM(_IM));
  Vec6T dvJ;
  dHdq.setZero();
  _DvDqM.setZero();
  _DaDqM.setZero();
  _DfDqM.setZero();
  _DdvJDqM.setZero();
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    JointFunc<T>::JCALC(J,mapV2CV(_qM),mapCV<Vec>((Vec*)NULL),r,mapM<Mat6T>((Mat6T*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapV<Vec6T>((Vec6T*)NULL),mapV(dvJ),mapM<Mat6XT>((Mat6XT*)NULL),_DdvJDqM,mapM<Mat6XT>((Mat6XT*)NULL));
    if(J._parent>=0)
      _aM.col(i)=TRANSI6(_invTPM,i)*_aM.col(J._parent);
    else _aM.col(i).setZero();
    _aM.col(i)+=dvJ;
    _fM.col(i)=TRANSI6(IMCustom,i)*_aM.col(i);
    DaDq(body,i,mapCV<Vec6T>((Vec6T*)NULL));
    Df(body,i,mapCV<Vec6T>((Vec6T*)NULL),IMCustom);
  }
  for(sizeType i=body.nrJ()-1; i>=0; i--) {
    const Joint& J=body.joint(i);
    DHDq(body,i,i,dHdq);
    if(J._parent>=0) {
      _fM.col(J._parent)+=TRANSI6(_invTPM,i).transpose()*_fM.col(i);
      DfP(body,i,i,false);
    }
  }
}
template <typename T>
typename NEArticulatedGradientInfoMap<T>::Mat3X4T NEArticulatedGradientInfoMap<T>::getTrans(sizeType JID) const
{
  return fromSpatial<T>(spatialInv<T>(TRANSI6(_invTM,JID)));
}
//-------------------------------------------------------------CRBA
template <typename T>
void NEArticulatedGradientInfoMap<T>::CRBA(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0,bool useLTDL)
{
  RNEA(body,a0,fx,false);
  calcH(body);
  _invHM=_HM;
  if(useLTDL) {
    LTDL();
    LTDLSolve(_ddqM=tau0-_tauM);
  } else {
    LTL();
    LTLSolve(_ddqM=tau0-_tauM);
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::calcHInvH(const ArticulatedBody& body,bool useLTDL)
{
  calcH(body);
  if(_updateInvH)
    return;
  _invHM=_HM;
  _DtauDqM.setIdentity();
  if(useLTDL) {
    LTDL();
    LTDLSolve(_DtauDqM);
  } else {
    LTL();
    LTLSolve(_DtauDqM);
  }
  _invHM=_DtauDqM;
  _updateInvH=true;
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::calcHInner(const ArticulatedBody& body,MatTM H,Mat6XTCM IMCustom)
{
  if(!IMCustom.data())
    new (&IMCustom) Mat6XTCM(mapM2CM(_IM));
  _ICM=IMCustom;
  Mat6X3T F;
  H.setZero();
  for(sizeType i=body.nrJ()-1; i>=0; i--) {
    const Joint& J=body.joint(i);
    sizeType nrDOF=J.nrDOF();
    if(J._parent>=0)
      TRANSI6(_ICM,J._parent)+=TRANSI6(_invTPM,i).transpose()*(TRANSI6(_ICM,i)*TRANSI6(_invTPM,i));
    Eigen::Block<Mat6XTM> Si=_SM.template block(0,J._offDOF,6,nrDOF);
    Eigen::Block<Mat6X3T> Fb=F.template block(0,0,6,nrDOF);
    Fb=TRANSI6(_ICM,i)*Si;
    H.template block(J._offDOF,J._offDOF,nrDOF,nrDOF)=Si.transpose()*Fb;
    for(sizeType j=i; body.joint(j)._parent>=0;) {
      Fb=(TRANSI6(_invTPM,j).transpose()*Fb).eval();
      j=body.joint(j)._parent;
      const Joint& Jj=body.joint(j);
      sizeType nrDOFj=Jj.nrDOF();
      Eigen::Block<Mat6XTM> Sj=_SM.template block(0,Jj._offDOF,6,nrDOFj);
      H.template block(J._offDOF,Jj._offDOF,nrDOF,nrDOFj)=Fb.transpose()*Sj;
      if(nrDOF>0 && nrDOFj>0)
        H.template block(Jj._offDOF,J._offDOF,nrDOFj,nrDOF)=H.template block(J._offDOF,Jj._offDOF,nrDOF,nrDOFj).transpose();
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::calcH(const ArticulatedBody& body)
{
  if(_updateH)
    return;
  calcHInner(body,_HM);
  _updateH=true;
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::LTDL()
{
  //factorize
  T a;
  for(sizeType k=_invHM.rows()-1,i,j; k>=0; k--) {
    i=_sparsityM[k];
    while(i>=0) {
      a=_invHM(k,i)/_invHM(k,k);
      j=i;
      while(j>=0) {
        _invHM(i,j)-=a*_invHM(k,j);
        j=_sparsityM[j];
      }
      _invHM(k,i)=a;
      i=_sparsityM[i];
    }
  }
}
template <typename T>
template <typename M>
void NEArticulatedGradientInfoMap<T>::LTDLSolve(M x)
{
  //solve LT
  for(sizeType i=_invHM.rows()-1,j; i>=0; i--) {
    j=_sparsityM[i];
    while(j>=0) {
      x.row(j)-=_invHM(i,j)*x.row(i);
      j=_sparsityM[j];
    }
  }
  //solve D
  for(sizeType i=0; i<_invHM.rows(); i++)
    x.row(i)/=_invHM(i,i);
  //solve L
  for(sizeType i=0,j; i<_invHM.rows(); i++) {
    j=_sparsityM[i];
    while(j>=0) {
      x.row(i)-=_invHM(i,j)*x.row(j);
      j=_sparsityM[j];
    }
  }
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::LTL()
{
  //factorize
  for(sizeType k=_invHM.rows()-1,i,j; k>=0; k--) {
    _invHM(k,k)=std::sqrt(_invHM(k,k));
    i=_sparsityM[k];
    while(i>=0) {
      _invHM(k,i)/=_invHM(k,k);
      i=_sparsityM[i];
    }
    i=_sparsityM[k];
    while(i>=0) {
      j=i;
      while(j>=0) {
        _invHM(i,j)-=_invHM(k,i)*_invHM(k,j);
        j=_sparsityM[j];
      }
      i=_sparsityM[i];
    }
  }
}
template <typename T>
template <typename M>
void NEArticulatedGradientInfoMap<T>::LTLSolve(M x)
{
  //solve LT
  for(sizeType i=_invHM.rows()-1,j; i>=0; i--) {
    x.row(i)/=_invHM(i,i);
    j=_sparsityM[i];
    while(j>=0) {
      x.row(j)-=_invHM(i,j)*x.row(i);
      j=_sparsityM[j];
    }
  }
  //solve L
  for(sizeType i=0,j; i<_invHM.rows(); i++) {
    j=_sparsityM[i];
    while(j>=0) {
      x.row(i)-=_invHM(i,j)*x.row(j);
      j=_sparsityM[j];
    }
    x.row(i)/=_invHM(i,i);
  }
}
//-------------------------------------------------------------ABA
template <typename T>
void NEArticulatedGradientInfoMap<T>::ABA(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0)
{
  //reuse _ICM as IA
  _ICM=_IM;
  if(tau0.data())
    _tauM=tau0;
  else _tauM.setZero();
  Mat6T Ia;
  Vec6T vJ,dvJ,pa,aPrime;
  VecCM ddqM=mapCV<Vec>((Vec*)NULL);
  //pass-1
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    JointFunc<T>::JCALC(J,mapV2CV(_qM),mapV2CV(_dqM),ddqM,mapM<Mat6T>((Mat6T*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapV(vJ),mapV(dvJ),mapM<Mat6XT>((Mat6XT*)NULL),mapM<Mat6XT>((Mat6XT*)NULL),mapM<Mat6XT>((Mat6XT*)NULL));
    if(J._parent>=0)
      _vM.col(i)=TRANSI6(_invTPM,i)*_vM.col(J._parent);
    else _vM.col(i).setZero();
    _vM.col(i)+=vJ;
    //reuse _aM as c
    _aM.col(i)=dvJ+spatialCross<T>(_vM.col(i),vJ);
    //reuse _fM as pA
    _fM.col(i)=spatialCrossStar<T>(_vM.col(i),TRANSI6(_IM,i)*_vM.col(i));
    if(fx.data())
      _fM.col(i)-=spatialXStar<T>(TRANSI6(_invTM,i))*fx.col(i);
  }
  //pass-2
  for(sizeType i=body.nrJ()-1; i>=0; i--) {
    const Joint& J=body.joint(i);
    sizeType nrDOF=J.nrDOF();
    Eigen::Block<Mat6XTM> Si=_SM.template block(0,J._offDOF,6,nrDOF);
    Eigen::Block<Mat6XTM> Ui=_UM.template block(0,J._offDOF,6,nrDOF);
    Ui=TRANSI6(_ICM,i)*Si;
    Eigen::Block<MatTM> Di=_HM.template block(J._offDOF,J._offDOF,nrDOF,nrDOF);
    invert(Si.transpose()*Ui,Di);
    Eigen::Block<VecM,-1,1> taui=_tauM.template segment(J._offDOF,nrDOF);
    taui-=Si.transpose()*_fM.col(i);
    if(J._parent>=0) {
      Ia=TRANSI6(_ICM,i)-Ui*(Di*Ui.transpose());
      pa=_fM.col(i)+Ia*_aM.col(i)+Ui*(Di*taui);
      TRANSI6(_ICM,J._parent)+=TRANSI6(_invTPM,i).transpose()*Ia*TRANSI6(_invTPM,i);
      _fM.col(J._parent)+=TRANSI6(_invTPM,i).transpose()*pa;
    }
  }
  //pass-3
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    sizeType nrDOF=J.nrDOF();
    if(J._parent>=0) {
      aPrime=TRANSI6(_invTPM,i)*_aM.col(J._parent);
    } else {
      if(a0.data())
        aPrime=TRANSI6(_invTPM,i)*a0;
      else aPrime.setZero();
    }
    aPrime+=_aM.col(i);
    Eigen::Block<Mat6XTM> Si=_SM.template block(0,J._offDOF,6,nrDOF);
    Eigen::Block<Mat6XTM> Ui=_UM.template block(0,J._offDOF,6,nrDOF);
    Eigen::Block<MatTM> Di=_HM.template block(J._offDOF,J._offDOF,nrDOF,nrDOF);
    Eigen::Block<VecM,-1,1> taui=_tauM.template segment(J._offDOF,nrDOF);
    _ddqM.template segment(J._offDOF,nrDOF)=Di*(taui-Ui.transpose()*aPrime);
    _aM.col(i)=aPrime+Si*_ddqM.template segment(J._offDOF,nrDOF);
  }
}
template <typename T>
template <typename M,typename IM>
void NEArticulatedGradientInfoMap<T>::invert(M in,IM inv)
{
  if(in.rows()==0) {
  } else if(in.rows()==1) {
    inv(0,0)=1/in(0,0);
  } else if(in.rows()==2) {
    T det=1/(in(0,0)*in(1,1)-in(0,1)*in(1,0));
    inv(0,0)=in(1,1)*det;
    inv(1,1)=in(0,0)*det;
    inv(0,1)=in(0,1)*-det;
    inv(1,0)=in(1,0)*-det;
  } else {
#define get(I,J) in(I-1,J-1)
#define set(I,J,v) inv(I-1,J-1)=v;
    T detA=0,buffer;
    T& result = detA;
    result += get(1, 1) * get(2, 2) * get(3, 3);
    result += get(1, 2) * get(2, 3) * get(3, 1);
    result += get(1, 3) * get(2, 1) * get(3, 2);
    result -= get(1, 3) * get(2, 2) * get(3, 1);
    result -= get(1, 2) * get(2, 1) * get(3, 3);
    result -= get(1, 1) * get(2, 3) * get(3, 2);

    buffer = get(2, 2) * get(3, 3) - get(2, 3) * get(3, 2);
    set(1, 1, buffer / detA);
    buffer = get(1, 3) * get(3, 2) - get(1, 2) * get(3, 3);
    set(1, 2, buffer / detA);
    buffer = get(1, 2) * get(2, 3) - get(1, 3) * get(2, 2);
    set(1, 3, buffer / detA);

    buffer = get(2, 3) * get(3, 1) - get(2, 1) * get(3, 3);
    set(2, 1, buffer / detA);
    buffer = get(1, 1) * get(3, 3) - get(1, 3) * get(3, 1);
    set(2, 2, buffer / detA);
    buffer = get(1, 3) * get(2, 1) - get(1, 1) * get(2, 3);
    set(2, 3, buffer / detA);

    buffer = get(2, 1) * get(3, 2) - get(2, 2) * get(3, 1);
    set(3, 1, buffer / detA);
    buffer = get(1, 2) * get(3, 1) - get(1, 1) * get(3, 2);
    set(3, 2, buffer / detA);
    buffer = get(1, 1) * get(2, 2) - get(1, 2) * get(2, 1);
    set(3, 3, buffer / detA);
#undef get
#undef set
  }
}
//-------------------------------------------------------------ABADerivatives
template <typename T>
void NEArticulatedGradientInfoMap<T>::ABADerivatives(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0,bool useLTDL)
{
  RNEA(body,a0,fx,false);
  calcHInvH(body,useLTDL);
  if(tau0.data())
    _ddqM=_invHM*(tau0-_tauM);  //now we have ddq
  else _ddqM=_invHM*-_tauM;     //now we have ddq
  //compute derivatives of RNEA
  RNEADerivatives(body,a0,fx,true);
  _DtauDqM=-_invHM*_DtauDqM;
  _DtauDdqM=-_invHM*_DtauDdqM;
}
template <typename T>
void NEArticulatedGradientInfoMap<T>::ABADerivativesFD(const ArticulatedBody& body,Vec6TCM a0,Mat6XTCM fx,VecCM tau0,VecM q1,VecM dq1,VecM ddq1,T delta,bool useLTDL)
{
  RNEA(body,a0,fx,false);
  calcH(body);
  _invHM=_HM;
  _DtauDqM.setIdentity();
  if(useLTDL) {
    LTDL();
    LTDLSolve(_DtauDqM);
  } else {
    LTL();
    LTLSolve(_DtauDqM);
  }
  _invHM=_DtauDqM;
  if(tau0.data())
    _ddqM=_invHM*(tau0-_tauM);  //now we have ddq
  else _ddqM=_invHM*-_tauM;     //now we have ddq
  //compute derivatives of RNEA
  RNEADerivativesFD(body,a0,fx,q1,dq1,ddq1,delta,true);
  _DtauDqM=-_invHM*_DtauDqM;
  _DtauDdqM=-_invHM*_DtauDdqM;
}
template <typename T>
typename NEArticulatedGradientInfoMap<T>::MatTCM NEArticulatedGradientInfoMap<T>::getDddqDq() const
{
  return mapM2CM(_DtauDqM);
}
template <typename T>
typename NEArticulatedGradientInfoMap<T>::MatTCM NEArticulatedGradientInfoMap<T>::getDddqDdq() const
{
  return mapM2CM(_DtauDdqM);
}
template <typename T>
typename NEArticulatedGradientInfoMap<T>::MatTCM NEArticulatedGradientInfoMap<T>::getDddqDtau0() const
{
  return mapM2CM(_invHM);
}
template <typename T>
typename NEArticulatedGradientInfoMap<T>::MatTCM NEArticulatedGradientInfoMap<T>::getH(const ArticulatedBody& body)
{
  calcH(body);
  return mapM2CM(_HM);
}
template <typename T>
typename NEArticulatedGradientInfoMap<T>::MatTCM NEArticulatedGradientInfoMap<T>::getInvH(const ArticulatedBody& body)
{
  calcHInvH(body);
  return mapM2CM(_invHM);
}
//ArticulatedGradientInfo
template <typename T>
NEArticulatedGradientInfo<T>::NEArticulatedGradientInfo() {}
template <typename T>
NEArticulatedGradientInfo<T>::NEArticulatedGradientInfo(const NEArticulatedGradientInfo& other)
{
  operator=(other);
}
template <typename T>
NEArticulatedGradientInfo<T>::NEArticulatedGradientInfo(const ArticulatedBody& body,const Vec& q,const Vec& dq,const Vec& ddq)
{
  reset(body,q,dq,ddq);
}
template <typename T>
NEArticulatedGradientInfo<T>::NEArticulatedGradientInfo(const ArticulatedBody& body,const Vec& q,const Vec& dq)
{
  reset(body,q,dq);
}
template <typename T>
NEArticulatedGradientInfo<T>::NEArticulatedGradientInfo(const ArticulatedBody& body,const Vec& q)
{
  reset(body,q);
}
template <typename T>
NEArticulatedGradientInfo<T>& NEArticulatedGradientInfo<T>::operator=(const NEArticulatedGradientInfo& other)
{
  _sparsity=other._sparsity;
  _nonZero=other._nonZero;
  _q=other._q;
  _dq=other._dq;
  _ddq=other._ddq;
  _tau=other._tau;

  _invTP=other._invTP;
  _invT=other._invT;
  _S=other._S;
  _U=other._U;
  _v=other._v;
  _vP=other._vP;
  _a=other._a;
  _f=other._f;
  _IC=other._IC;
  _DvJDq=other._DvJDq;
  _DdvJDq=other._DdvJDq;
  _DdvJDdq=other._DdvJDdq;
  //JTrans,I
  _JTrans=other._JTrans;
  _I=other._I;

  _DtauDq=other._DtauDq;
  _DtauDdq=other._DtauDdq;
  _DvDq=other._DvDq;
  _DvDdq=other._DvDdq;
  _DaDq=other._DaDq;
  _DaDdq=other._DaDdq;
  _DfDq=other._DfDq;
  _DfDdq=other._DfDdq;
  _H=other._H;
  _invH=other._invH;
  resetPtr();
  return *this;
}
template <typename T>
void NEArticulatedGradientInfo<T>::reset(const ArticulatedBody& body)
{
  sizeType nrJ=body.nrJ();
  sizeType nrDOF=body.nrDOF();
  _nonZero.resize(nrJ);
  _q.resize(nrDOF);
  _dq.resize(nrDOF);
  _ddq.resize(nrDOF);
  _tau.resize(nrDOF);

  _invTP.resize(6,nrJ*6);
  _invT.resize(6,nrJ*6);
  _S.resize(6,nrDOF);
  _U.resize(6,nrDOF);
  _v.resize(6,nrJ);
  _vP.resize(6,nrDOF);
  _a.resize(6,nrJ);
  _f.resize(6,nrJ);
  _IC.resize(6,nrJ*6);
  _DvJDq.resize(6,nrDOF);
  _DdvJDq.resize(6,nrDOF);
  _DdvJDdq.resize(6,nrDOF);

  _DtauDq.resize(nrDOF,nrDOF);
  _DtauDdq.resize(nrDOF,nrDOF);
  _DvDq.resize(nrJ*6,nrDOF);
  _DvDdq.resize(nrJ*6,nrDOF);
  _DaDq.resize(nrJ*6,nrDOF);
  _DaDdq.resize(nrJ*6,nrDOF);
  _DfDq.resize(nrJ*6,nrDOF);
  _DfDdq.resize(nrJ*6,nrDOF);
  _H.resize(nrDOF,nrDOF);
  _invH.resize(nrDOF,nrDOF);
  reorthogonalize(body);
  resetPtr();
}
template <typename T>
void NEArticulatedGradientInfo<T>::reset(const ArticulatedBody& body,const Vec& q,const Vec& dq,const Vec& ddq)
{
  reset(body);
  NEArticulatedGradientInfoMap<T>::reset(body,mapV(q),mapV(dq),mapV(ddq));
}
template <typename T>
void NEArticulatedGradientInfo<T>::reset(const ArticulatedBody& body,const Vec& q,const Vec& dq)
{
  reset(body);
  NEArticulatedGradientInfoMap<T>::reset(body,mapV(q),mapV(dq));
}
template <typename T>
void NEArticulatedGradientInfo<T>::reset(const ArticulatedBody& body,const Vec& q)
{
  reset(body);
  NEArticulatedGradientInfoMap<T>::reset(body,mapV(q));
}
template <typename T>
void NEArticulatedGradientInfo<T>::reorthogonalize(const ArticulatedBody& body)
{
  sizeType nrJ=body.nrJ();
  sizeType nrDOF=body.nrDOF();
  //this re-orthogonalize JTrans, which is essential for MPFR precision
  //this is because JTrans was originally calculated using low-precision
  if(_JTrans.rows()!=6 || _JTrans.cols()!=nrJ*6) {
    _JTrans.resize(6,nrJ*6);
    for(sizeType j=0; j<nrJ; j++) {
      Mat3X4T JTrans=body.joint(j)._trans.template cast<T>();
      JTrans.col(0)/=std::sqrt(JTrans.col(0).squaredNorm());
      JTrans.col(1)-=JTrans.col(1).dot(JTrans.col(0))*JTrans.col(0);
      JTrans.col(1)/=std::sqrt(JTrans.col(1).squaredNorm());
      JTrans.col(2)=JTrans.col(0).cross(JTrans.col(1));
      TRANSI6(_JTrans,j)=toSpatial<T>(JTrans);
    }
  }
  //this computes and stores inertial matrix
  if(_I.rows()!=6 || _I.cols()!=nrJ*6) {
    Mat6T P;
    _I.resize(6,nrJ*6);
    P.setZero();
    P.template block<3,3>(0,3).setIdentity();
    P.template block<3,3>(3,0).setIdentity();
    for(sizeType j=0; j<nrJ; j++)
      TRANSI6(_I,j)=P*body.joint(j).getMass().template cast<T>()*P;
  }
  //this computes and stores sparsity
  _sparsity.setConstant(nrDOF,-1);
  for(sizeType i=body.nrJ()-1; i>=0; i--) {
    const Joint& J=body.joint(i);
    for(sizeType d=J._offDOF+J.nrDOF()-1; d>J._offDOF; d--)
      _sparsity[d]=d-1;
    if(J._parent>=0 && J._offDOF>0 && J._typeJoint!=Joint::FIX_JOINT) {
      sizeType ip=J._parent;
      while(ip>=0)
        if(body.joint(ip)._typeJoint!=Joint::FIX_JOINT)
          break;
        else ip=body.joint(ip)._parent;
      if(ip>=0) {
        const Joint& JP=body.joint(ip);
        _sparsity[J._offDOF]=JP._offDOF+JP.nrDOF()-1;
      }
    }
  }
}
template <typename T>
typename NEArticulatedGradientInfo<T>::Mat3XT NEArticulatedGradientInfo<T>::getTrans() const
{
  Mat3XT ret=Mat3XT::Zero(3,_invTM.cols()/6*4);
  for(sizeType i=0; i<_invTM.cols()/6; i++)
    TRANSI(ret,i)=NEArticulatedGradientInfoMap<T>::getTrans(i);
  return ret;
}
template <typename T>
void NEArticulatedGradientInfo<T>::debug(const ArticulatedBody& body)
{
  sizeType nrJ=body.nrJ();
  sizeType nrDOF=body.nrDOF();

  Vec6T a,a2,a0;
  Vec3T x,xG,xG2,v,v2;
  MatT Sc,Sc2,ScInvH,ScInvHScT,DHDqr,DScDqr,DScDqTr;
  Vec q,dq,ddq,tau0,deltaq,r,r2;
  Vec q1,dq1,ddq1,dfx,dtg;
  Mat6XT fx,fx2;
  Mat3X4T TCoef;
  fx.resize(6,nrJ);
  q.resize(nrDOF);
  dq.resize(nrDOF);
  ddq.resize(nrDOF);
  dtg.resize(nrDOF);
  tau0.resize(nrDOF);
  deltaq.resize(nrDOF);
  r.resize(nrDOF);
  q1.resize(nrDOF);
  dq1.resize(nrDOF);
  ddq1.resize(nrDOF);
  DHDqr.resize(nrDOF,nrDOF);

  std::vector<sizeType> joints;
  PBDArticulatedGradientInfo<T> GPBD;
  NEArticulatedGradientInfo G,G2;
  DEFINE_NUMERIC_DELTA_T(T)
  T DELTA_MULT=1E3f*DELTA;
  for(sizeType JID=0; JID<nrJ; JID++) {
    INFO("-------------------------------------------------------------DebugNEArticulatedGradientInfo")
    //calculate info
    q.setRandom();
    dq.setRandom();
    ddq.setRandom();
    tau0.setRandom();
    deltaq.setRandom();
    TCoef.setRandom();
    dtg.setZero();

    G.reset(body,q,dq,ddq);
    G.DTG(JID,body,TCoef,mapV(dtg));
    G.RNEA(body,mapCV<Vec6T>((Vec6T*)NULL),mapCM<Mat6XT>((Mat6XT*)NULL));
    GPBD.reset(body,q);
    T errT=(spatialInv<T>(TRANSI6(G._invTM,JID))-toSpatial<T>(TRANSI(GPBD._TM,JID))).squaredNorm();
    DEBUG_GRADIENT("T",std::sqrt(TRANSI6(G._invTM,JID).squaredNorm()),std::sqrt(errT))
    G2.reset(body,q+deltaq*DELTA);
    DEBUG_GRADIENT("DTG",dtg.dot(deltaq),dtg.dot(deltaq)-((G2.NEArticulatedGradientInfoMap<T>::getTrans(JID)-G.NEArticulatedGradientInfoMap<T>::getTrans(JID))/DELTA*TCoef.transpose()).trace())

    G2.reset(body,q+dq*DELTA,dq+ddq*DELTA,ddq);
    G2.RNEA(body,mapCV<Vec6T>((Vec6T*)NULL),mapCM<Mat6XT>((Mat6XT*)NULL));
    x.setRandom(3);
    v=spatialVel<T>(G._vM.col(JID),x);
    xG=spatialApplyTransInv<T>(TRANSI6(G._invTM,JID),x);
    xG2=spatialApplyTransInv<T>(TRANSI6(G2._invTM,JID),x);
    v2=(spatialApplyTrans<T>(TRANSI6(G._invTM,JID),xG2)-spatialApplyTrans<T>(TRANSI6(G._invTM,JID),xG))/DELTA;
    DEBUG_GRADIENT("v",std::sqrt(v.squaredNorm()),std::sqrt((v-v2).squaredNorm()))

    a=G._aM.col(JID);
    a2=(G2._vM.col(JID)-G._vM.col(JID))/DELTA;
    DEBUG_GRADIENT("a",std::sqrt(q.squaredNorm()),std::sqrt((a-a2).squaredNorm()))

    a0.setRandom();
    G.reset(body,q,dq,ddq);
    G.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>((Mat6XT*)NULL));
    for(sizeType j=0; j<nrJ; j++)
      fx.col(j)=-TRANSI6(G._invTM,j).transpose()*TRANSI6(G._IM,j)*TRANSI6(G._invTM,j)*a0;
    G2.reset(body,q,dq,ddq);
    G2.RNEA(body,mapCV<Vec6T>(NULL),mapCM<Mat6XT>(fx));
    DEBUG_GRADIENT("g",std::sqrt(G._tauM.squaredNorm()),std::sqrt((G._tauM-G2._tauM).squaredNorm()))
    a0.setRandom();
    fx.setRandom();

    G.reset(body,q,dq);
    G.CRBA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0),false);
    G.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx));
    DEBUG_GRADIENT("CRBA-LTL",std::sqrt(tau0.squaredNorm()),std::sqrt((G._tau-tau0).squaredNorm()))

    G.reset(body,q,dq);
    G.CRBA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0),true);
    G.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx));
    DEBUG_GRADIENT("CRBA-LTDL",std::sqrt(tau0.squaredNorm()),std::sqrt((G._tau-tau0).squaredNorm()))

    G.reset(body,q,dq);
    G.ABA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0));
    G.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx));
    DEBUG_GRADIENT("ABA",std::sqrt(tau0.squaredNorm()),std::sqrt((G._tau-tau0).squaredNorm()))

    G.reset(body,q,dq);
    G.ABADerivatives(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0),true);
    //G.ABADerivativesFD(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0),mapV(q1),mapV(dq1),mapV(ddq1),DELTA,true);
    G2.reset(body,q+deltaq*DELTA_MULT,dq);
    G2.ABA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0));
    DEBUG_GRADIENT("ABADerivatives-q",std::sqrt((G.getDddqDq()*deltaq).squaredNorm()),std::sqrt((G.getDddqDq()*deltaq-(G2._ddq-G._ddq)/DELTA_MULT).squaredNorm()))

    G.reset(body,q,dq);
    G.ABADerivatives(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0),true);
    //G.ABADerivativesFD(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0),mapV(q1),mapV(dq1),mapV(ddq1),DELTA,true);
    G2.reset(body,q,dq+deltaq*DELTA_MULT);
    G2.ABA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0));
    DEBUG_GRADIENT("ABADerivatives-dq",std::sqrt((G.getDddqDdq()*deltaq).squaredNorm()),std::sqrt((G.getDddqDdq()*deltaq-(G2._ddq-G._ddq)/DELTA_MULT).squaredNorm()))

    G.reset(body,q,dq,ddq);
    G.RNEADerivatives(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),true);
    G2.reset(body,q+deltaq*DELTA,dq,ddq);
    G2.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),true);
    DEBUG_GRADIENT("DvDq",std::sqrt((G._DvDqM.block(JID*6,0,6,nrDOF)*deltaq).squaredNorm()),std::sqrt((G._DvDqM.block(JID*6,0,6,nrDOF)*deltaq-(G2._vM.col(JID)-G._vM.col(JID))/DELTA).squaredNorm()));
    DEBUG_GRADIENT("DaDq",std::sqrt((G._DaDqM.block(JID*6,0,6,nrDOF)*deltaq).squaredNorm()),std::sqrt((G._DaDqM.block(JID*6,0,6,nrDOF)*deltaq-(G2._aM.col(JID)-G._aM.col(JID))/DELTA).squaredNorm()));
    DEBUG_GRADIENT("DfDq",std::sqrt((G._DfDqM.block(JID*6,0,6,nrDOF)*deltaq).squaredNorm()),std::sqrt((G._DfDqM.block(JID*6,0,6,nrDOF)*deltaq-(G2._fM.col(JID)-G._fM.col(JID))/DELTA).squaredNorm()));
    DEBUG_GRADIENT("DtauDq",std::sqrt((G._DtauDqM*deltaq).squaredNorm()),std::sqrt((G._DtauDqM*deltaq-(G2._tauM-G._tauM)/DELTA).squaredNorm()));
    G2.reset(body,q,dq+deltaq*DELTA,ddq);
    G2.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),true);
    DEBUG_GRADIENT("DvDdq",std::sqrt((G._DvDdqM.block(JID*6,0,6,nrDOF)*deltaq).squaredNorm()),std::sqrt((G._DvDdqM.block(JID*6,0,6,nrDOF)*deltaq-(G2._vM.col(JID)-G._vM.col(JID))/DELTA).squaredNorm()));
    DEBUG_GRADIENT("DaDdq",std::sqrt((G._DaDdqM.block(JID*6,0,6,nrDOF)*deltaq).squaredNorm()),std::sqrt((G._DaDdqM.block(JID*6,0,6,nrDOF)*deltaq-(G2._aM.col(JID)-G._aM.col(JID))/DELTA).squaredNorm()));
    DEBUG_GRADIENT("DfDdq",std::sqrt((G._DfDdqM.block(JID*6,0,6,nrDOF)*deltaq).squaredNorm()),std::sqrt((G._DfDdqM.block(JID*6,0,6,nrDOF)*deltaq-(G2._fM.col(JID)-G._fM.col(JID))/DELTA).squaredNorm()));
    DEBUG_GRADIENT("DtauDdq",std::sqrt((G._DtauDdqM*deltaq).squaredNorm()),std::sqrt((G._DtauDdqM*deltaq-(G2._tauM-G._tauM)/DELTA).squaredNorm()));

    joints.clear();
    for(sizeType j=0; j<JID; j++)
      if(RandEngine::randR01()>0.5)
        joints.push_back(j);
    joints.push_back(JID);
    ColiCM jointsM(&joints[0],(sizeType)joints.size());
    std::ostringstream oss;
    oss << jointsM.transpose();
    INFOV("Testing random joints: %s!",oss.str().c_str())

    r.setRandom();
    dfx.setRandom(jointsM.size()*6);
    fx2=fx;
    for(sizeType j=0; j<(sizeType)joints.size(); j++)
      fx2.col(jointsM[j])+=dfx.template segment<6>(j*6)*DELTA;
    Sc.resize(jointsM.size()*6,nrDOF);
    Sc2.resize(jointsM.size()*6,nrDOF);
    ScInvH.resize(jointsM.size()*6,nrDOF);
    ScInvHScT.resize(jointsM.size()*6,jointsM.size()*6);
    G.reset(body,q,dq,ddq);
    G2.reset(body,q,dq,ddq);
    G.Sc(body,jointsM,mapM(Sc));
    G.DHDq(body,mapM(DHDqr),mapCV(r));
    G.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),true);
    G2.RNEA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx2),true);
    DEBUG_GRADIENT("Sc-tau",std::sqrt((-Sc.transpose()*dfx).squaredNorm()),std::sqrt((-Sc.transpose()*dfx-(G2._tau-G._tau)/DELTA).squaredNorm()))
    G.calcHInvH(body);
    G.ScInvH(body,jointsM,mapM(Sc),mapM(ScInvH));
    G.ScInvHScT(body,jointsM,mapM(Sc),mapM(ScInvH),mapM(ScInvHScT));
    G.ABA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx),mapCV(tau0));
    G2.ABA(body,mapCV<Vec6T>(a0),mapCM<Mat6XT>(fx2),mapCV(tau0));
    DEBUG_GRADIENT("ScInvH",std::sqrt(ScInvH.squaredNorm()),std::sqrt((ScInvH-Sc*G._invH).squaredNorm()))
    DEBUG_GRADIENT("ScInvHScT",std::sqrt(ScInvHScT.squaredNorm()),std::sqrt((ScInvHScT-Sc*G._invH*Sc.transpose()).squaredNorm()))
    DEBUG_GRADIENT("ScInvH-ddq",std::sqrt((ScInvH.transpose()*dfx).squaredNorm()),std::sqrt((ScInvH.transpose()*dfx-(G2._ddq-G._ddq)/DELTA).squaredNorm()))

    DScDqr.setRandom(jointsM.size()*6,nrDOF);
    DScDqTr.setRandom(nrDOF,nrDOF);
    G.reset(body,q,dq,ddq);
    G.calcH(body);
    G.Sc(body,jointsM,mapM(Sc));
    G.DHDq(body,mapM(DHDqr),mapCV(r));
    G.DScDq(body,jointsM,mapM(DScDqr),mapCV(r));
    G.DScDqT(body,jointsM,mapM(DScDqTr),mapCV(dfx));
    G2.reset(body,q+deltaq*DELTA,dq,ddq);
    G2.calcH(body);
    G2.Sc(body,jointsM,mapM(Sc2));
    DEBUG_GRADIENT("DHDq",std::sqrt((DHDqr*deltaq).squaredNorm()),std::sqrt((DHDqr*deltaq-(G2._HM-G._HM)*r/DELTA).squaredNorm()))
    DEBUG_GRADIENT("DScDqr",std::sqrt((DScDqr*deltaq).squaredNorm()),std::sqrt((DScDqr*deltaq-(Sc2-Sc)*r/DELTA).squaredNorm()))
    DEBUG_GRADIENT("DScDqTr",std::sqrt((DScDqTr*deltaq).squaredNorm()),std::sqrt((DScDqTr*deltaq-(Sc2-Sc).transpose()*dfx/DELTA).squaredNorm()))
  }
}
template <typename T>
void NEArticulatedGradientInfo<T>::resetPtr()
{
  new (&_sparsityM)    ColiM(mapV(_sparsity));
  new (&_nonZeroM)     ColiM(mapV(_nonZero));
  new (&_qM)           VecM(mapV(_q));
  new (&_dqM)          VecM(mapV(_dq));
  new (&_ddqM)         VecM(mapV(_ddq));
  new (&_tauM)         VecM(mapV(_tau));

  new (&_invTPM)       Mat6XTM(mapM(_invTP));
  new (&_invTM)        Mat6XTM(mapM(_invT));
  new (&_SM)           Mat6XTM(mapM(_S));
  new (&_UM)           Mat6XTM(mapM(_U));
  new (&_vM)           Mat6XTM(mapM(_v));
  new (&_vPM)          Mat6XTM(mapM(_vP));
  new (&_aM)           Mat6XTM(mapM(_a));
  new (&_fM)           Mat6XTM(mapM(_f));
  new (&_ICM)          Mat6XTM(mapM(_IC));
  new (&_DvJDqM)       Mat6XTM(mapM(_DvJDq));
  new (&_DdvJDqM)      Mat6XTM(mapM(_DdvJDq));
  new (&_DdvJDdqM)     Mat6XTM(mapM(_DdvJDdq));
  //JTrans,I
  new (&_JTransM)      Mat6XTM(mapM(_JTrans));
  new (&_IM)           Mat6XTM(mapM(_I));

  new (&_DtauDqM)      MatTM(mapM(_DtauDq));
  new (&_DtauDdqM)     MatTM(mapM(_DtauDdq));
  new (&_DvDqM)        MatTM(mapM(_DvDq));
  new (&_DvDdqM)       MatTM(mapM(_DvDdq));
  new (&_DaDqM)        MatTM(mapM(_DaDq));
  new (&_DaDdqM)       MatTM(mapM(_DaDdq));
  new (&_DfDqM)        MatTM(mapM(_DfDq));
  new (&_DfDdqM)       MatTM(mapM(_DfDdq));
  new (&_HM)           MatTM(mapM(_H));
  new (&_invHM)        MatTM(mapM(_invH));
}
//instance
PRJ_BEGIN
template struct NEArticulatedGradientInfoMap<double>;
template struct NEArticulatedGradientInfo<double>;
#ifdef ALL_TYPES
template struct NEArticulatedGradientInfoMap<__float128>;
template struct NEArticulatedGradientInfo<__float128>;
template struct NEArticulatedGradientInfoMap<mpfr::mpreal>;
template struct NEArticulatedGradientInfo<mpfr::mpreal>;
#endif
PRJ_END
