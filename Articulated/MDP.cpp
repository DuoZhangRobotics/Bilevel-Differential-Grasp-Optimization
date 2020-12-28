#ifdef ENVIRONMENT_SUPPORT
#include "MDP.h"
#include <Utils/Utils.h>
#include <Utils/CrossSpatialUtil.h>

PRJ_BEGIN

//ExternalWrench
template <typename T>
ExternalWrench<T>::ExternalWrench()
{
  _DBDq.clear();
  _DBDX.assign(3,std::vector<Mat6XT,Eigen::aligned_allocator<Mat6XT>>(4,Mat6XT::Zero(6,0)));
  _jointId=-1;
  _B.resize(6,0);
  _w.resize(0);
}
//DebugWrenchConstructor
template <typename T>
DebugWrenchConstructor<T>::DebugWrenchConstructor(const std::vector<ExternalWrench<T>>& ref):EnvWrenchConstructor<T>(std::shared_ptr<Environment<T>>()),_refWrench(ref) {}
template <typename T>
void DebugWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,const NEArticulatedGradientInfo<T>& info,bool)
{
  externalWrench=_refWrench;
  for(sizeType i=0; i<(sizeType)externalWrench.size(); i++) {
    ExternalWrench<T>& wi=externalWrench[i];
    if(!wi._DBDq.empty()) {
      for(sizeType d=0; d<(sizeType)wi._DBDq.size(); d++) {
        const Mat6XT& DBDq=wi._DBDq[d].second;
        sizeType qid=wi._DBDq[d].first;
        wi._B+=DBDq*info._qM[qid];
      }
    } else if(wi._DBDX[0][0].size()>0) {
      Mat3X4T TJ=info.NEArticulatedGradientInfoMap<T>::getTrans(wi._jointId);
      for(sizeType r=0; r<3; r++)
        for(sizeType c=0; c<4; c++)
          wi._B+=wi._DBDX[r][c]*TJ(r,c);
    }
  }
}
template <typename T>
NEArticulatedGradientInfo<T> DebugWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,const ArticulatedBody& body,VecCM qM)
{
  NEArticulatedGradientInfo<T> info(body,qM);
  operator()(externalWrench,info,true);
  return info;
}
//MDP
template <typename T>
MDP<T>::MDP(const ArticulatedBody& body,Options& ops):_solLQP(ops),_body(body) {}
template <typename T>
MDP<T>::MDP(const MDP<T>& other):_solLQP(other._solLQP),_body(other._body),_info(other._info) {}
template <typename T>
bool MDP<T>::solveMDPQPF(MatTM DdqHatDq,MatTM DdqHatDdq,bool profileMDPError)
{
  if(DdqHatDq.data())
    DdqHatDq.setZero();
  if(DdqHatDdq.data())
    DdqHatDdq.setIdentity();
  _deltaDq.setZero(_info._dqM.size());
  _dqMDP=_info._dqM;
  if(_externalWrench.empty())
    return false;
  _joints.clear();
  _solLQP._foot.resize((sizeType)_externalWrench.size());
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++) {
    _joints.push_back(_externalWrench[i]._jointId);
    _solLQP._foot[i]=_externalWrench[i]._B.cols();
  }
  if(_solLQP._foot.sum()==0)
    return false;
  _info.calcHInvH(_body,true);
  assembleHcQP();
  //solve W
  bool succ;
  _lastW=_solLQP.solve(succ,_lastW.size()==_solLQP._c.size()?&_lastW:NULL);
  if(!succ && profileMDPError) {
    WARNING("MultiPrecisionQP failed!")
    _solLQP.writeProb("error.dat");
    //exit(-1);
  }
  //apply
  Bf(mapCV(_lastW));
  _deltaDq=_ScInvH.transpose()*_Bf;
  _dqMDP=_info._dqM+_deltaDq;
  //differentiable
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  if(DdqHatDq.data() && DdqHatDdq.data()) {
    sizeType N=_info._qM.size();
    //DdqHatDdq
    BT(_Sc,_BTSc);
    _BTScInvH=_BTSc*_info._invHM;
    _solLQP.computeFGH(_solLQP._muFinal,_lastW,NULL,&_Sigma);
    SolveNewton<T>::inverse(_Sigma,_invHLQP);
    _Sigma=_BTScInvH.transpose()*_invHLQP;
    DdqHatDdq=MatT::Identity(N,N)-_Sigma*_BTSc;
    //DdqHatDq: part-1
    _DScDq.resize(6*jointsM.size(),N);
    DBDqT(_Sc*_dqMDP,_DBDqT);
    _info.DScDq(_body,jointsM,mapM(_DScDq),mapCV(_dqMDP));
    BT(_DScDq,_BTDScDq);
    DdqHatDq=-_Sigma*(_BTDScDq+_DBDqT);
    //DdqHatDq: part-2
    _DHDq.resize(N,N);
    _DScDqT.resize(N,N);
    _info.DHDq(_body,mapM(_DHDq),mapCV(_deltaDq));
    _info.DScDqT(_body,jointsM,mapM(_DScDqT),mapCV(_Bf));
    DBDq(_lastW,_DBDq);
    DdqHatDq+=DdqHatDdq*(_info._invHM*(-_DHDq+_DScDqT+_Sc.transpose()*_DBDq));
  }
  return true;
}
template <typename T>
bool MDP<T>::solveMDPQPI(MatTM DdqHatDq,MatTM DdqHatDdq,bool profileMDPError)
{
  if(DdqHatDq.data())
    DdqHatDq.setZero();
  if(DdqHatDdq.data())
    DdqHatDdq.setZero();
  _deltaDq.setZero(_info._dqM.size());
  _dqMDP.setZero(0);
  if(_externalWrench.empty())
    return false;
  _joints.clear();
  _solLQP._foot.resize((sizeType)_externalWrench.size());
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++) {
    _joints.push_back(_externalWrench[i]._jointId);
    _solLQP._foot[i]=_externalWrench[i]._B.cols();
  }
  if(_solLQP._foot.sum()==0)
    return false;
  _info.calcHInvH(_body,true);
  assembleHcQP();
  //solve W
  bool succ;
  _lastW=_solLQP.solve(succ,_lastW.size()==_solLQP._c.size()?&_lastW:NULL);
  if(!succ && profileMDPError) {
    WARNING("MultiPrecisionQP failed!")
    _solLQP.writeProb("error.dat");
    //exit(-1);
  }
  //apply
  Bf(mapCV(_lastW));
  _deltaDq=_Sc.transpose()*_Bf;
  //differentiable
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  if(DdqHatDq.data() && DdqHatDdq.data()) {
    _dqMDP=_info._dqM+_ScInvH.transpose()*_Bf;
    sizeType N=_info._qM.size();
    //DdqHatDdq
    BT(_Sc,_BTSc);
    _solLQP.computeFGH(_solLQP._muFinal,_lastW,NULL,&_Sigma);
    SolveNewton<T>::inverse(_Sigma,_invHLQP);
    _Sigma=_BTSc.transpose()*_invHLQP;
    DdqHatDdq=-_Sigma*_BTSc;
    //DdqHatDq: part-1
    _DScDq.resize(6*jointsM.size(),N);
    DBDqT(_Sc*_dqMDP,_DBDqT);
    _info.DScDq(_body,jointsM,mapM(_DScDq),mapCV(_dqMDP));
    BT(_DScDq,_BTDScDq);
    DdqHatDq=-_Sigma*(_BTDScDq+_DBDqT);
    //DdqHatDq: part-2
    _DHDq.resize(N,N);
    _DScDqT.resize(N,N);
    _dqMDP-=_info._dqM;
    _info.DHDq(_body,mapM(_DHDq),mapCV(_dqMDP));
    _info.DScDqT(_body,jointsM,mapM(_DScDqT),mapCV(_Bf));
    DBDq(_lastW,_DBDq);
    _Sigma=_DScDqT+_Sc.transpose()*_DBDq; //temporary reuse Sigma
    DdqHatDq+=_Sigma+DdqHatDdq*_info._invHM*(-_DHDq+_Sigma);
    _dqMDP.setZero(0);  //flag that _dqMDP is not used
  }
  return true;
}
template <typename T>
bool MDP<T>::solveMDPLPF(MatTM DddqDq,MatTM DddqDdq,bool profileMDPError)
{
  if(DddqDq.data())
    DddqDq.setZero();
  if(DddqDdq.data())
    DddqDdq.setZero();
  _deltaDq.setZero(_info._dqM.size());
  _dqMDP.setZero(0);
  if(_externalWrench.empty())
    return false;
  _joints.clear();
  _solLQP._foot.resize((sizeType)_externalWrench.size());
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++) {
    _joints.push_back(_externalWrench[i]._jointId);
    _solLQP._foot[i]=_externalWrench[i]._B.cols();
  }
  if(_solLQP._foot.sum()==0)
    return false;
  _info.calcHInvH(_body,true);
  assembleHcLP();
  //solve W
  bool succ;
  _lastW=_solLQP.solve(succ,_lastW.size()==_solLQP._c.size()?&_lastW:NULL);
  if(!succ && profileMDPError) {
    WARNING("MultiPrecisionLP failed!")
    _solLQP.writeProb("error.dat");
    //exit(-1);
  }
  //apply
  Bf(mapCV(_lastW));
  _deltaDq=_ScInvH.transpose()*_Bf;
  //differentiable
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  if(DddqDq.data() && DddqDdq.data()) {
    sizeType N=_info._qM.size();
    //DddqDdq
    BT(_Sc,_BTSc);
    _BTScInvH=_BTSc*_info._invHM;
    _solLQP.computeFGH(_solLQP._muFinal,_lastW,NULL,&_Sigma);
    SolveNewtonLP<T>::inverse(_solLQP._muFinal,_lastW,_invHLQP,_solLQP._foot);
    _Sigma=_BTScInvH.transpose()*_invHLQP;
    DddqDdq=-_Sigma*_BTSc;
    //DddqDq: part-1
    _DScDq.resize(6*jointsM.size(),N);
    DBDqT(_Sc*_info._dqM,_DBDqT);
    _info.DScDq(_body,jointsM,mapM(_DScDq),mapV2CV(_info._dqM));
    BT(_DScDq,_BTDScDq);
    DddqDq=-_Sigma*(_BTDScDq+_DBDqT);
    //DddqDq: part-2
    _DHDq.resize(N,N);
    _DScDqT.resize(N,N);
    _info.DHDq(_body,mapM(_DHDq),mapCV(_deltaDq));
    _info.DScDqT(_body,jointsM,mapM(_DScDqT),mapCV(_Bf));
    DBDq(_lastW,_DBDq);
    DddqDq+=_info._invHM*(-_DHDq+_DScDqT+_Sc.transpose()*_DBDq);
  }
  return true;
}
template <typename T>
bool MDP<T>::solveMDPLPI(MatTM DddqDq,MatTM DddqDdq,bool profileMDPError)
{
  if(DddqDq.data())
    DddqDq.setZero();
  if(DddqDdq.data())
    DddqDdq.setZero();
  _deltaDq.setZero(_info._dqM.size());
  _dqMDP.setZero(0);
  if(_externalWrench.empty())
    return false;
  _joints.clear();
  _solLQP._foot.resize((sizeType)_externalWrench.size());
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++) {
    _joints.push_back(_externalWrench[i]._jointId);
    _solLQP._foot[i]=_externalWrench[i]._B.cols();
  }
  if(_solLQP._foot.sum()==0)
    return false;
  _info.calcH(_body);
  assembleHcLP();
  //solve W
  bool succ;
  _lastW=_solLQP.solve(succ,_lastW.size()==_solLQP._c.size()?&_lastW:NULL);
  if(!succ && profileMDPError) {
    WARNING("MultiPrecisionLPI failed!")
    _solLQP.writeProb("error.dat");
    //exit(-1);
  }
  //apply
  Bf(mapCV(_lastW));
  _deltaDq=_Sc.transpose()*_Bf;
  //differentiable
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  if(DddqDq.data() && DddqDdq.data()) {
    sizeType N=_info._qM.size();
    //DddqDdq
    BT(_Sc,_BTSc);
    _solLQP.computeFGH(_solLQP._muFinal,_lastW,NULL,&_Sigma);
    SolveNewtonLP<T>::inverse(_solLQP._muFinal,_lastW,_invHLQP,_solLQP._foot);
    _Sigma=_BTSc.transpose()*_invHLQP;
    DddqDdq=-_Sigma*_BTSc;
    //DddqDq: part-1
    _DScDq.resize(6*jointsM.size(),N);
    DBDqT(_Sc*_info._dqM,_DBDqT);
    _info.DScDq(_body,jointsM,mapM(_DScDq),mapV2CV(_info._dqM));
    BT(_DScDq,_BTDScDq);
    DddqDq=-_Sigma*(_BTDScDq+_DBDqT);
    //DddqDq: part-2
    _DScDqT.resize(N,N);
    _info.DScDqT(_body,jointsM,mapM(_DScDqT),mapCV(_Bf));
    DBDq(_lastW,_DBDq);
    DddqDq+=_DScDqT+_Sc.transpose()*_DBDq;
  }
  return true;
}
template <typename T>
MDP<T>& MDP<T>::operator=(const MDP<T>& other)
{
  _solLQP=other._solLQP;
  _info=other._info;
  return *this;
}
template <typename T>
void MDP<T>::reset(Options& ops)
{
  _solLQP.reset(ops);
}
template <typename T>
void MDP<T>::randomize(T scale,bool QP)
{
  std::vector<sizeType> joints;
  sizeType nrJ=_body.nrJ();
  while(joints.size()<5) {
    joints.clear();
    for(sizeType j=0; j<nrJ; j++)
      if(RandEngine::randR01()>0.5)
        joints.push_back(j);
  }
  ColiCM jointsM(&joints[0],(sizeType)joints.size());
  std::ostringstream oss;
  oss << jointsM.transpose();
  INFOV("Testing random joints: %s!",oss.str().c_str())

  _externalWrench.resize(jointsM.size());
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++) {
    sizeType nrB=RandEngine::randI(1,10);
    _externalWrench[i]._jointId=jointsM[i];
    _externalWrench[i]._B=MatT::Random(6,nrB)*scale;
    _externalWrench[i]._DBDq.clear();
    _externalWrench[i]._DBDX[0][0].resize(6,0);
    if(RandEngine::randR01()>0.5) {
      for(sizeType d=0; d<_body.nrDOF(); d++)
        if(RandEngine::randR01()>0.5)
          _externalWrench[i]._DBDq.push_back(std::make_pair(d,MatT::Random(6,nrB)*scale));
    } else {
      for(sizeType r=0; r<3; r++)
        for(sizeType c=0; c<4; c++)
          _externalWrench[i]._DBDX[r][c]=MatT::Random(6,nrB)*scale;
    }
  }
  _joints.clear();
  _solLQP._foot.resize((sizeType)_externalWrench.size());
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++) {
    _joints.push_back(_externalWrench[i]._jointId);
    _solLQP._foot[i]=_externalWrench[i]._B.cols();
  }
  _info.reset(_body,Vec::Random(_body.nrDOF()),Vec::Random(_body.nrDOF()));
  _info.calcHInvH(_body,true);
  if(QP)
    assembleHcQP();
  else assembleHcLP();
}
template <typename T>
void MDP<T>::readAndTestProb()
{
  if(!exists("error.dat"))
    return;
  _solLQP.readAndTestProb("error.dat");
  exit(-1);
}
template <typename T>
void MDP<T>::assembleHcQP()
{
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  //this is H
  _Sc.resize(6*jointsM.size(),_info._HM.cols());
  _ScInvH.resize(6*jointsM.size(),_info._HM.cols());
  _ScInvHScT.resize(6*jointsM.size(),6*jointsM.size());
  _info.ScInvHScT(_body,jointsM,mapM<MatT>(_Sc),mapM<MatT>(_ScInvH),mapM<MatT>(_ScInvHScT));
  //solQP: H,c
  Vec Scdq=_Sc*_info._dqM;
  _solLQP._c.resize(_solLQP._foot.sum());
  _solLQP._H.resize(_solLQP._foot.sum(),_solLQP._foot.sum());
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); offi+=_externalWrench[i]._B.cols(),i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    _solLQP._c.segment(offi,_externalWrench[i]._B.cols())=wi._B.transpose()*Scdq.template segment<6>(i*6);
    //Hii
    Eigen::Block<MatT> Hii=_solLQP._H.block(offi,offi,wi._B.cols(),wi._B.cols());
    Hii=wi._B.transpose()*_ScInvHScT.template block<6,6>(i*6,i*6)*wi._B;
    //Hij/Hji
    for(sizeType j=0,offj=0; j<i; offj+=_externalWrench[j]._B.cols(),j++) {
      const ExternalWrench<T>& wj=_externalWrench[j];
      Eigen::Block<MatT> Hij=_solLQP._H.block(offi,offj,wi._B.cols(),wj._B.cols());
      Eigen::Block<MatT> Hji=_solLQP._H.block(offj,offi,wj._B.cols(),wi._B.cols());
      Hij=wi._B.transpose()*_ScInvHScT.template block<6,6>(i*6,j*6)*wj._B;
      Hji=Hij.transpose();
    }
  }
}
template <typename T>
void MDP<T>::assembleHcLP()
{
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  //this is H
  _Sc.resize(6*jointsM.size(),_info._HM.cols());
  _ScInvH.resize(6*jointsM.size(),_info._HM.cols());
  _ScInvHScT.resize(0,0);
  _info.ScInvH(_body,jointsM,mapM<MatT>(_Sc),mapM<MatT>(_ScInvH));
  //solQP: H,c
  Vec Scdq=_Sc*_info._dqM;
  _solLQP._c.resize(_solLQP._foot.sum());
  _solLQP._H.resize(0,0);
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); offi+=_externalWrench[i]._B.cols(),i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    _solLQP._c.segment(offi,_externalWrench[i]._B.cols())=wi._B.transpose()*Scdq.template segment<6>(i*6);
  }
}
template <typename T>
void MDP<T>::assembleHcLPI()
{
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  //this is H
  _Sc.resize(6*jointsM.size(),_info._HM.cols());
  _ScInvH.resize(0,0);
  _ScInvHScT.resize(0,0);
  _info.Sc(_body,jointsM,mapM<MatT>(_Sc));
  //solQP: H,c
  Vec Scdq=_Sc*_info._dqM;
  _solLQP._c.resize(_solLQP._foot.sum());
  _solLQP._H.resize(0,0);
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); offi+=_externalWrench[i]._B.cols(),i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    _solLQP._c.segment(offi,_externalWrench[i]._B.cols())=wi._B.transpose()*Scdq.template segment<6>(i*6);
  }
}
template <typename T>
sizeType MDP<T>::nW() const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)_externalWrench.size(); i++)
    ret+=_externalWrench[i]._B.cols();
  return ret;
}
template <typename T>
typename MDP<T>::MatT MDP<T>::B() const
{
  MatT ret;
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  ret.setZero(jointsM.size()*6,_solLQP._foot.sum());
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); offi+=_externalWrench[i]._B.cols(),i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    ret.block(i*6,offi,wi._B.rows(),wi._B.cols())=wi._B;
  }
  return ret;
}
template <typename T>
void MDP<T>::Bf(const Vec& f)
{
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  _Bf=Vec::Zero(6*jointsM.size());
  for(sizeType i=0,off=0; i<(sizeType)_externalWrench.size(); off+=_externalWrench[i]._B.cols(),i++) {
    ExternalWrench<T>& wi=_externalWrench[i];
    wi._w=f.segment(off,wi._B.cols());
    _Bf.template segment<6>(i*6)=wi._B*wi._w;
  }
}
template <typename T>
void MDP<T>::BT(const MatT& b,MatT& BTb) const
{
  BTb.resize(_solLQP._foot.sum(),b.cols());
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); offi+=_externalWrench[i]._B.cols(),i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    BTb.block(offi,0,wi._B.cols(),b.cols())=wi._B.transpose()*b.block(i*6,0,6,b.cols());
  }
}
template <typename T>
void MDP<T>::DBDq(const Vec& b,MatT& DBDqb) const
{
  Mat3X4T DBbDX;
  ColiCM jointsM(&_joints[0],(sizeType)_joints.size());
  DBDqb.setZero(jointsM.size()*6,_info._qM.size());
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); offi+=_externalWrench[i]._B.cols(),i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    VecCM bBlk(&(b.coeffRef(offi)),wi._B.cols());
    if(!wi._DBDq.empty()) {
      for(sizeType d=0; d<(sizeType)wi._DBDq.size(); d++) {
        const Mat6XT& DBDq=wi._DBDq[d].second;
        sizeType qid=wi._DBDq[d].first;
        DBDqb.template block<6,1>(i*6,qid)+=DBDq*bBlk;
      }
    } else if(wi._DBDX[0][0].size()>0) {
      for(sizeType d=0; d<6; d++) {
        for(sizeType r=0; r<3; r++)
          for(sizeType c=0; c<4; c++)
            DBbDX(r,c)=wi._DBDX[r][c].row(d).dot(bBlk);
        DXDq(DBDqb,i*6+d,wi._jointId,DBbDX);
      }
    }
  }
}
template <typename T>
void MDP<T>::DBDqT(const Vec& b,MatT& DBDqTb) const
{
  Mat3X4T DBTbDX;
  DBDqTb.setZero(_solLQP._foot.sum(),_info._qM.size());
  for(sizeType i=0,offi=0; i<(sizeType)_externalWrench.size(); offi+=_externalWrench[i]._B.cols(),i++) {
    const ExternalWrench<T>& wi=_externalWrench[i];
    Vec6TCM bBlk(&(b.coeffRef(i*6)),6);
    if(!wi._DBDq.empty()) {
      for(sizeType d=0; d<(sizeType)wi._DBDq.size(); d++) {
        const Mat6XT& DBDq=wi._DBDq[d].second;
        sizeType qid=wi._DBDq[d].first;
        DBDqTb.block(offi,qid,wi._B.cols(),1)+=DBDq.transpose()*bBlk;
      }
    } else if(wi._DBDX[0][0].size()>0) {
      for(sizeType d=0; d<wi._B.cols(); d++) {
        for(sizeType r=0; r<3; r++)
          for(sizeType c=0; c<4; c++)
            DBTbDX(r,c)=wi._DBDX[r][c].col(d).dot(bBlk);
        DXDq(DBDqTb,offi+d,wi._jointId,DBTbDX);
      }
    }
  }
}
template <typename T>
void MDP<T>::DXDq(MatT& DBDq,sizeType r,sizeType JID,const Mat3X4T& DBDX) const
{
  Vec6T V;
  Mat6T TJ;
  Mat3X4T RtJ=_info.NEArticulatedGradientInfoMap<T>::getTrans(JID);
  V.template segment<3>(3)=CTR(DBDX);
  V.template segment<3>(0)=invCrossMatTrace<T>(ROT(RtJ)*ROT(DBDX).transpose()+CTR(RtJ)*CTR(DBDX).transpose());
  while(JID>=0) {
    const Joint& J=_body.joint(JID);
    sizeType nrDOF=J.nrDOF();
    TJ=spatialInv<T>(TRANSI6(_info._invTM,JID));
    for(sizeType d=0; d<nrDOF; d++)
      DBDq(r,J._offDOF+d)+=V.dot(TJ*_info._SM.col(J._offDOF+d));
    JID=J._parent;
  }
}
template <typename T>
T MDP<T>::kineticEnergyInfo(const Vec& w)
{
  Bf(w);
  _deltaDq=_ScInvH.transpose()*_Bf;
  _dqMDP=_info._dqM+_deltaDq;
  return _dqMDP.dot(_info._HM*_dqMDP)/2;
}
template <typename T>
T MDP<T>::kineticEnergyAssembled(const Vec& w) const
{
  T E=(_solLQP._c+_solLQP._H*w/2).dot(w);
  E+=_info._dqM.dot(_info._HM*_info._dqM)/2;
  return E;
}
template <typename T>
void MDP<T>::debugWrenchConstructor(const ArticulatedBody& body,EnvWrenchConstructor<T>& debugWrench,sizeType nrIter)
{
  Options ops;
  MDP<T> mdp(body,ops);

  sizeType N=mdp._body.nrDOF();
  Vec qdq,qdq2,deltaq,f,fT,Bf,Bf2,BTf,BTf2;
  MatT DBDqb,DBDqTb;
  DEFINE_NUMERIC_DELTA_T(T)
  for(sizeType i=0; i<nrIter; i++) {
    INFO("-------------------------------------------------------------debugWrenchConstructor")
    qdq.setRandom(N);
    deltaq.setRandom(N);
    qdq2=qdq+deltaq*DELTA;
    mdp._info.reset(body,qdq);
    debugWrench(mdp._externalWrench,mdp._info,true);
    mdp._joints.clear();
    mdp._solLQP._foot.resize((sizeType)mdp._externalWrench.size());
    for(sizeType i=0; i<(sizeType)mdp._externalWrench.size(); i++) {
      mdp._joints.push_back(mdp._externalWrench[i]._jointId);
      mdp._solLQP._foot[i]=mdp._externalWrench[i]._B.cols();
    }
    f.setRandom(6*(sizeType)mdp._externalWrench.size());
    fT.setRandom(mdp._solLQP._foot.sum());
    Bf=mdp.B()*fT;
    BTf=mdp.B().transpose()*f;
    mdp.DBDq(fT,DBDqb);
    mdp.DBDqT(f,DBDqTb);

    mdp._info.reset(body,qdq2);
    debugWrench(mdp._externalWrench,mdp._info,false);
    mdp._joints.clear();
    mdp._solLQP._foot.resize((sizeType)mdp._externalWrench.size());
    for(sizeType i=0; i<(sizeType)mdp._externalWrench.size(); i++) {
      mdp._joints.push_back(mdp._externalWrench[i]._jointId);
      mdp._solLQP._foot[i]=mdp._externalWrench[i]._B.cols();
    }
    Bf2=mdp.B()*fT;
    BTf2=mdp.B().transpose()*f;
    DEBUG_GRADIENT("DBDqb",std::sqrt((DBDqb*deltaq).squaredNorm()),std::sqrt((DBDqb*deltaq-(Bf2-Bf)/DELTA).squaredNorm()))
    DEBUG_GRADIENT("DBDqTb",std::sqrt((DBDqTb*deltaq).squaredNorm()),std::sqrt((DBDqTb*deltaq-(BTf2-BTf)/DELTA).squaredNorm()))
  }
}
template <typename T>
void MDP<T>::debug(const ArticulatedBody& body,sizeType nrIter,T scale)
{
  DEFINE_NUMERIC_DELTA_T(T)
  T DELTA_MULT=std::max<T>(DELTA,1e-20f);

  Options ops;
  MDP<T> mdp(body,ops);
  ops.setOptions<MultiPrecisionLQP<T>,bool>("callback",false);
  mdp._solLQP.reset(ops);

  Mat3X4T DBDX;
  sizeType N=mdp._body.nrDOF();
  T K0,K1,BX,BX2;
  Vec w,qdq,qdq2,deltaq,f,fT,Bf,Bf2,BTf,BTf2;
  MatT DdqHatDq,DdqHatDdq,dqMDP,dqMDP2,DBDqb,DBDqTb,DBDqRow;
  for(sizeType i=0; i<nrIter; i++) {
    INFO("-------------------------------------------------------------DebugMDP")
    mdp.randomize(scale,false);
    mdp.randomize(scale,true);
    mdp._solLQP.debugGradient();
    DebugWrenchConstructor<T> debugWrench(mdp._externalWrench);
    w=Vec::Random(mdp._solLQP._c.size());
    T KI=mdp.kineticEnergyInfo(mapCV(w));
    T KA=mdp.kineticEnergyAssembled(mapCV(w));
    DEBUG_GRADIENT("MDPSystem",KI,KI-KA)

    bool succ;
    mdp._lastW.setZero(mdp._solLQP._c.size());
    K0=mdp.kineticEnergyInfo(mapCV(mdp._lastW));
    mdp._lastW=mdp._solLQP.solve(succ);//ASSERT(succ)
    K1=mdp.kineticEnergyInfo(mapCV(mdp._lastW));
    std::ostringstream oss2;
    oss2 << "MDPSolve kinetic energy before=" << K0 << " after=" << K1;
    INFO(oss2.str().c_str())

    qdq.resize(N);
    DBDX.setRandom();
    DBDqRow.setZero(1,N);
    deltaq.setRandom(N);
    qdq2=qdq+deltaq*DELTA;
    f.setRandom(6*(sizeType)mdp._externalWrench.size());
    fT.setRandom(mdp._solLQP._foot.sum());
    mdp._info=debugWrench(mdp._externalWrench,body,mapCV(qdq));
    Bf=mdp.B()*fT;
    BTf=mdp.B().transpose()*f;
    BX=(mdp._info.NEArticulatedGradientInfoMap<T>::getTrans(mdp._externalWrench.back()._jointId)*DBDX.transpose()).trace();
    mdp.DXDq(DBDqRow,0,mdp._externalWrench.back()._jointId,DBDX);
    mdp.DBDq(fT,DBDqb);
    mdp.DBDqT(f,DBDqTb);
    mdp._info=debugWrench(mdp._externalWrench,body,mapCV(qdq2));
    Bf2=mdp.B()*fT;
    BTf2=mdp.B().transpose()*f;
    BX2=(mdp._info.NEArticulatedGradientInfoMap<T>::getTrans(mdp._externalWrench.back()._jointId)*DBDX.transpose()).trace();
    DEBUG_GRADIENT("DXDq",DBDqRow.row(0).dot(deltaq),DBDqRow.row(0).dot(deltaq)-(BX2-BX)/DELTA)
    DEBUG_GRADIENT("DBDqb",std::sqrt((DBDqb*deltaq).squaredNorm()),std::sqrt((DBDqb*deltaq-(Bf2-Bf)/DELTA).squaredNorm()))
    DEBUG_GRADIENT("DBDqTb",std::sqrt((DBDqTb*deltaq).squaredNorm()),std::sqrt((DBDqTb*deltaq-(BTf2-BTf)/DELTA).squaredNorm()))

    mdp._lastW.setZero(0);
    qdq.setRandom(N*2);
    deltaq.setRandom(N);
    DdqHatDq.resize(N,N);
    DdqHatDdq.resize(N,N);

    std::vector<std::string> str;
    str.push_back("solveMDPQPF");
    str.push_back("solveMDPQPI");
    str.push_back("solveMDPLPF");
    str.push_back("solveMDPLPI");
    typedef bool(MDP<T>::*proc)(MatTM DdqHatDq,MatTM DdqHatDdq,bool profileMDPError);
    proc mainProc[4]= {&MDP<T>::solveMDPQPF,&MDP<T>::solveMDPQPI,&MDP<T>::solveMDPLPF,&MDP<T>::solveMDPLPI};
    for(sizeType pass=0; pass<4; pass++) {
      proc p=mainProc[pass];
      //mdp
      qdq2=qdq;
      qdq2.segment(N,N)+=deltaq*DELTA_MULT;
      mdp._info.reset(mdp._body,qdq.segment(0,N),qdq.segment(N,N));
      debugWrench(mdp._externalWrench,mdp._info);
      (mdp.*p)(mapM(DdqHatDq),mapM(DdqHatDdq),false);
      dqMDP=endsWith(str[pass],"solveMDPQPF")?mdp._dqMDP:mdp._deltaDq;
      //mdp2
      mdp._info.reset(mdp._body,qdq2.segment(0,N),qdq2.segment(N,N));
      debugWrench(mdp._externalWrench,mdp._info);
      (mdp.*p)(mapM<MatT>((MatT*)NULL),mapM<MatT>((MatT*)NULL),false);
      dqMDP2=endsWith(str[pass],"solveMDPQPF")?mdp._dqMDP:mdp._deltaDq;
      DEBUG_GRADIENT(str[pass]+"-DdqHatDdq",std::sqrt((DdqHatDdq*deltaq).squaredNorm()),std::sqrt((DdqHatDdq*deltaq-(dqMDP2-dqMDP)/DELTA_MULT).squaredNorm()))
      //mdp2
      qdq2=qdq;
      qdq2.segment(0,N)+=deltaq*DELTA_MULT;
      mdp._info.reset(mdp._body,qdq2.segment(0,N),qdq2.segment(N,N));
      debugWrench(mdp._externalWrench,mdp._info);
      (mdp.*p)(mapM<MatT>((MatT*)NULL),mapM<MatT>((MatT*)NULL),false);
      dqMDP2=endsWith(str[pass],"solveMDPQPF")?mdp._dqMDP:mdp._deltaDq;
      DEBUG_GRADIENT(str[pass]+"-DdqHatDq",std::sqrt((DdqHatDq*deltaq).squaredNorm()),std::sqrt((DdqHatDq*deltaq-(dqMDP2-dqMDP)/DELTA_MULT).squaredNorm()))
    }
  }
}
//instance
#define INSTANTIATE(T)    \
template class ExternalWrench<T>;   \
template class DebugWrenchConstructor<T>;   \
template class MDP<T>;
INSTANTIATE(double)
#ifdef ALL_TYPES
INSTANTIATE(__float128)
INSTANTIATE(mpfr::mpreal)
#endif
PRJ_END
#endif
