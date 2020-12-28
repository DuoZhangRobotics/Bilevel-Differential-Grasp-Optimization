#include <Utils/SpatialRotationUtil.h>
#include "JointFunc.h"

USE_PRJ_NAMESPACE

//JointFunc
//PBD
template <typename T>
typename JointFunc<T>::Mat3X4T JointFunc<T>::TDTDDT(const Joint& J,VecCM x,Mat3XTM DT,Mat3XTM DDT)
{
  Mat3X4T JLT;
  JLT.setIdentity();
  if(J._typeJoint == Joint::TRANS_3D) {
    CTR(JLT)=x.template segment<3>(J._offDOF);
    DT.col(J._offDOF+0)=Vec3T::UnitX();
    DT.col(J._offDOF+1)=Vec3T::UnitY();
    DT.col(J._offDOF+2)=Vec3T::UnitZ();
  } else if(J._typeJoint == Joint::TRANS_2D) {
    CTR(JLT)=Vec3T(x[J._offDOF+0],x[J._offDOF+1],0);
    DT.col(J._offDOF+0)=Vec3T::UnitX();
    DT.col(J._offDOF+1)=Vec3T::UnitY();
  } else if(J._typeJoint == Joint::TRANS_1D) {
    CTR(JLT)=Vec3T(x[J._offDOF+0],0,0);
    DT.col(J._offDOF+0)=Vec3T::UnitX();
  } else if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    Eigen::Map<Vec3T> DJTmp[3]= {
      Eigen::Map<Vec3T>(&(DT.coeffRef(0,J._offDOF+0))),
      Eigen::Map<Vec3T>(&(DT.coeffRef(0,J._offDOF+1))),
      Eigen::Map<Vec3T>(&(DT.coeffRef(0,J._offDOF+2))),
    };
    Eigen::Map<Vec3T> DDTTmp[9]= {
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+0))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+1))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+2))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+3))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+4))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+5))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+6))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+7))),
      Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT+8))),
    };
    if(J._typeJoint == Joint::ROT_3D_EXP)
      ROT(JLT)=expWGradV<T,Eigen::Map<Vec3T> >(x.template segment<3>(J._offDOF),DJTmp,DDTTmp);
    else ROT(JLT)=eulerX1Y3Z2<T,Eigen::Map<Vec3T> >(x.template segment<3>(J._offDOF),DJTmp,DDTTmp);
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    Eigen::Map<Vec3T> DJTmp[2]= {
      Eigen::Map<Vec3T>(&(DT.coeffRef(0,J._offDOF+0))),
      Eigen::Map<Vec3T>(&(DT.coeffRef(0,J._offDOF+1))),
    };
    Eigen::Map<Vec3T> DDTTmp=Eigen::Map<Vec3T>(&(DDT.coeffRef(0,J._offDDT)));
    ROT(JLT)=expWYZ<T,Eigen::Map<Vec3T> >(x[J._offDOF+0],x[J._offDOF+1],DJTmp+0,DJTmp+1,&DDTTmp);
  } else if(J._typeJoint == Joint::HINGE_JOINT) {
    Eigen::Map<Vec3T> DJTmp=Eigen::Map<Vec3T>(&(DT.coeffRef(0,J._offDOF)));
    ROT(JLT)=expWZ<T,Eigen::Map<Vec3T> >(x[J._offDOF+0],&DJTmp);
  } else if(J._typeJoint == Joint::FIX_JOINT) {
  } else {
    ASSERT_MSGV(false,"Unknown joint type in %s, J._typeJoint=%d!",__FUNCTION__,J._typeJoint)
  }
  return JLT;
}
template <typename T>
void JointFunc<T>::DDT(const Joint& J,std::function<void(sizeType,sizeType,T)> hess,Mat3XTCM DDT,const Vec3T& coefRss)
{
  if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    for(int d=0; d<3; d++)
      for(int d2=0; d2<3; d2++)
        hess(J._offDOF+d,J._offDOF+d2,DDT.col(J._offDDT+d+d2*3).dot(coefRss));
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    hess(J._offDOF+1,J._offDOF,DDT.col(J._offDDT).dot(coefRss));
  }
}
template <typename T>
void JointFunc<T>::DDT(const Joint& J,MatTM hess,Mat3XTCM DDT,const Vec3T& coefRss)
{
  JointFunc<T>::DDT(J,[&](sizeType row,sizeType col,T val) {
    hess(row,col)+=val;
  },DDT,coefRss);
}
template <typename T>
typename JointFunc<T>::Vec3T JointFunc<T>::DDT(const Joint& J,Mat3XTCM DDT,sizeType R,sizeType C)
{
  if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    return DDT.col(J._offDDT+(R-J._offDOF)+(C-J._offDOF)*3);
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    if(R-J._offDOF==1 && C-J._offDOF==0)
      return DDT.col(J._offDDT);
    else return Vec3T::Zero();
  } else return Vec3T::Zero();
}
template <typename T>
void JointFunc<T>::DDDTLambda(const Joint& J,VecCM x,VecCM lambda,Mat3XTM DDDTLambda)
{
  if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    Eigen::Map<Vec3T> DDDTLambdaTmp[9]= {
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+0))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+1))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+2))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+3))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+4))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+5))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+6))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+7))),
      Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT+8))),
    };
    if(J._typeJoint == Joint::ROT_3D_EXP)
      expWGradVLambda<T,Eigen::Map<Vec3T> >(x.template segment<3>(J._offDOF),lambda.template segment<3>(J._offDOF),DDDTLambdaTmp);
    else eulerX1Y3Z2Lambda<T,Eigen::Map<Vec3T> >(x.template segment<3>(J._offDOF),lambda.template segment<3>(J._offDOF),DDDTLambdaTmp);
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    Eigen::Map<Vec3T> DDDTLambdaTmp=Eigen::Map<Vec3T>(&(DDDTLambda.coeffRef(0,J._offDDT)));
    expWYZLambda<T,Eigen::Map<Vec3T> >(x[J._offDOF+1],lambda[J._offDOF+1],DDDTLambdaTmp);
  }
}
template <typename T>
void JointFunc<T>::DDTLambda(const Joint& J,Mat3XTM hess,Mat3XTCM DDT,const VecCM& lambda)
{
  if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    for(int d=0; d<3; d++)
      hess.col(J._offDOF+d)=DDT.template block<3,3>(0,J._offDDT+d*3)*lambda.template segment<3>(J._offDOF);
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    hess.col(J._offDOF+0)=DDT.col(J._offDDT)*lambda[J._offDOF+1];
    hess.col(J._offDOF+1).setZero();
  } else hess.block(0,J._offDOF,3,J.nrDOF()).setZero();
}
//NE
template <typename T>
void JointFunc<T>::JCALC(const Joint& J,const VecCM q,const VecCM dq,const VecCM ddq,
                         Mat6TM XJ,Mat6XTM S,Vec6TM vJ,Vec6TM dvJ,Mat6XTM DvJDq,Mat6XTM DdvJDq,Mat6XTM DdvJDdq)
{
  if(J._typeJoint == Joint::TRANS_3D) {
    if(vJ.data()) {
      vJ.setZero();
      vJ.template segment<3>(3)=dq.template segment<3>(J._offDOF);
    }
    if(dvJ.data()) {
      dvJ.setZero();
      if(ddq.data())
        dvJ.template segment<3>(3)=ddq.template segment<3>(J._offDOF);
    }
    if(XJ.data()) {
      XJ.setIdentity();
      SPATIAL_BLK10(XJ)=cross<T>(Vec3T(q[J._offDOF+0],q[J._offDOF+1],q[J._offDOF+2]));
    }
    if(S.data()) {
      S.template block<3,3>(0,J._offDOF).setZero();
      S.template block<3,3>(3,J._offDOF).setIdentity();
    }
    Mat6XTM DvJDqM(DvJDq.data()?&(DvJDq.coeffRef(0,J._offDOF)):NULL,6,3,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDqM(DdvJDq.data()?&(DdvJDq.coeffRef(0,J._offDOF)):NULL,6,3,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDdqM(DdvJDdq.data()?&(DdvJDdq.coeffRef(0,J._offDOF)):NULL,6,3,Eigen::OuterStride<>(S.outerStride()));
    if(DvJDq.data())
      DvJDqM.setZero();
    if(DdvJDq.data())
      DdvJDqM.setZero();
    if(DdvJDdq.data())
      DdvJDdqM.setZero();
  } else if(J._typeJoint == Joint::TRANS_2D) {
    if(vJ.data()) {
      vJ.setZero();
      vJ.template segment<2>(3)=dq.template segment<2>(J._offDOF);
    }
    if(dvJ.data()) {
      dvJ.setZero();
      if(ddq.data())
        dvJ.template segment<2>(3)=ddq.template segment<2>(J._offDOF);
    }
    if(XJ.data()) {
      XJ.setIdentity();
      SPATIAL_BLK10(XJ)=cross<T>(Vec3T(q[J._offDOF+0],q[J._offDOF+1],0));
    }
    if(S.data()) {
      S.template block<3,2>(0,J._offDOF).setZero();
      S.template block<3,2>(3,J._offDOF).setIdentity();
    }
    Mat6XTM DvJDqM(DvJDq.data()?&(DvJDq.coeffRef(0,J._offDOF)):NULL,6,2,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDqM(DdvJDq.data()?&(DdvJDq.coeffRef(0,J._offDOF)):NULL,6,2,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDdqM(DdvJDdq.data()?&(DdvJDdq.coeffRef(0,J._offDOF)):NULL,6,2,Eigen::OuterStride<>(S.outerStride()));
    if(DvJDq.data())
      DvJDqM.setZero();
    if(DdvJDq.data())
      DdvJDqM.setZero();
    if(DdvJDdq.data())
      DdvJDdqM.setZero();
  } else if(J._typeJoint == Joint::TRANS_1D) {
    if(vJ.data()) {
      vJ.setZero();
      vJ.template segment<1>(3)=dq.template segment<1>(J._offDOF);
    }
    if(dvJ.data()) {
      dvJ.setZero();
      if(ddq.data())
        dvJ.template segment<1>(3)=ddq.template segment<1>(J._offDOF);
    }
    if(XJ.data()) {
      XJ.setIdentity();
      SPATIAL_BLK10(XJ)=cross<T>(Vec3T(q[J._offDOF+0],0,0));
    }
    if(S.data()) {
      S.template block<3,1>(0,J._offDOF).setZero();
      S.template block<3,1>(3,J._offDOF).setIdentity();
    }
    Mat6XTM DvJDqM(DvJDq.data()?&(DvJDq.coeffRef(0,J._offDOF)):NULL,6,1,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDqM(DdvJDq.data()?&(DdvJDq.coeffRef(0,J._offDOF)):NULL,6,1,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDdqM(DdvJDdq.data()?&(DdvJDdq.coeffRef(0,J._offDOF)):NULL,6,1,Eigen::OuterStride<>(S.outerStride()));
    if(DvJDq.data())
      DvJDqM.setZero();
    if(DdvJDq.data())
      DdvJDqM.setZero();
    if(DdvJDdq.data())
      DdvJDdqM.setZero();
  } else if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    VecCM qM(&(q.coeffRef(J._offDOF+0)),3);
    VecCM dqM(dq.data()?&(dq.coeffRef(J._offDOF+0)):NULL,3);
    VecCM ddqM(ddq.data()?&(ddq.coeffRef(J._offDOF+0)):NULL,3);
    Mat6XTM SM(S.data()?&(S.coeffRef(0,J._offDOF)):NULL,6,3,S.outerStride());
    Mat6XTM DvJDqM(DvJDq.data()?&(DvJDq.coeffRef(0,J._offDOF)):NULL,6,3,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDqM(DdvJDq.data()?&(DdvJDq.coeffRef(0,J._offDOF)):NULL,6,3,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDdqM(DdvJDdq.data()?&(DdvJDdq.coeffRef(0,J._offDOF)):NULL,6,3,Eigen::OuterStride<>(S.outerStride()));
    if(J._typeJoint==Joint::ROT_3D_EXP)
      spatialExpW<T,VecCM,Vec6TM,Mat6TM,Mat6XTM,Mat6XTM>
      (qM,dq.data()?&dqM:NULL,ddq.data()?&ddqM:NULL,
       XJ.data()?&XJ:NULL,
       S.data()?&SM:NULL,
       vJ.data()?&vJ:NULL,
       dvJ.data()?&dvJ:NULL,
       DvJDqM.data()?&DvJDqM:NULL,
       DdvJDqM.data()?&DdvJDqM:NULL,
       DdvJDdqM.data()?&DdvJDdqM:NULL);
    else
      spatialEulerX1Y3Z2<T,VecCM,Vec6TM,Mat6TM,Mat6XTM,Mat6XTM>
      (qM,dq.data()?&dqM:NULL,ddq.data()?&ddqM:NULL,
       XJ.data()?&XJ:NULL,
       S.data()?&SM:NULL,
       vJ.data()?&vJ:NULL,
       dvJ.data()?&dvJ:NULL,
       DvJDqM.data()?&DvJDqM:NULL,
       DdvJDqM.data()?&DdvJDqM:NULL,
       DdvJDdqM.data()?&DdvJDdqM:NULL);
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    VecCM qM(&(q.coeffRef(J._offDOF+0)),2);
    VecCM dqM(dq.data()?&(dq.coeffRef(J._offDOF+0)):NULL,2);
    VecCM ddqM(ddq.data()?&(ddq.coeffRef(J._offDOF+0)):NULL,2);
    Mat6XTM SM(S.data()?&(S.coeffRef(0,J._offDOF)):NULL,6,2,S.outerStride());
    Mat6XTM DvJDqM(DvJDq.data()?&(DvJDq.coeffRef(0,J._offDOF)):NULL,6,2,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDqM(DdvJDq.data()?&(DdvJDq.coeffRef(0,J._offDOF)):NULL,6,2,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDdqM(DdvJDdq.data()?&(DdvJDdq.coeffRef(0,J._offDOF)):NULL,6,2,Eigen::OuterStride<>(S.outerStride()));
    spatialExpWYZ<T,VecCM,Vec6TM,Mat6TM,Mat6XTM,Mat6XTM>
    (qM,dq.data()?&dqM:NULL,ddq.data()?&ddqM:NULL,
     XJ.data()?&XJ:NULL,
     S.data()?&SM:NULL,
     vJ.data()?&vJ:NULL,
     dvJ.data()?&dvJ:NULL,
     DvJDqM.data()?&DvJDqM:NULL,
     DdvJDqM.data()?&DdvJDqM:NULL,
     DdvJDdqM.data()?&DdvJDdqM:NULL);
  } else if(J._typeJoint == Joint::HINGE_JOINT) {
    const T* ddqM=ddq.data()?&(ddq.coeffRef(J._offDOF+0)):NULL;
    Mat6XTM SM(S.data()?&(S.coeffRef(0,J._offDOF)):NULL,6,1,S.outerStride());
    Mat6XTM DvJDqM(DvJDq.data()?&(DvJDq.coeffRef(0,J._offDOF)):NULL,6,1,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDqM(DdvJDq.data()?&(DdvJDq.coeffRef(0,J._offDOF)):NULL,6,1,Eigen::OuterStride<>(S.outerStride()));
    Mat6XTM DdvJDdqM(DdvJDdq.data()?&(DdvJDdq.coeffRef(0,J._offDOF)):NULL,6,1,Eigen::OuterStride<>(S.outerStride()));
    spatialExpWZ<T,Vec6TM,Mat6TM,Mat6XTM,Mat6XTM>
    (q[J._offDOF+0],dq.data()?(dq.data()+J._offDOF):NULL,ddqM,
     XJ.data()?&XJ:NULL,
     S.data()?&SM:NULL,
     vJ.data()?&vJ:NULL,
     dvJ.data()?&dvJ:NULL,
     DvJDqM.data()?&DvJDqM:NULL,
     DdvJDqM.data()?&DdvJDqM:NULL,
     DdvJDdqM.data()?&DdvJDdqM:NULL);
  } else if(J._typeJoint == Joint::FIX_JOINT) {
    if(vJ.data())
      vJ.setZero();
    if(dvJ.data())
      dvJ.setZero();
    if(XJ.data())
      XJ.setIdentity();
  } else {
    ASSERT_MSGV(false,"Unknown joint type in %s, J._typeJoint=%d!",__FUNCTION__,J._typeJoint)
  }
}
template <typename T>
void JointFunc<T>::DSTDqf(const Joint& J,const VecCM q,MatTM DtauDq,const Vec6T& f)
{
  sizeType nrDOF=J.nrDOF();
  MatTM DtauDqM(&(DtauDq.coeffRef(J._offDOF,J._offDOF)),nrDOF,nrDOF,DtauDq.outerStride());
  if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    VecCM qM(&(q.coeffRef(J._offDOF+0)),3);
    if(J._typeJoint == Joint::ROT_3D_EXP)
      spatialExpWDSTDqf<T,VecCM,MatTM,Vec6T>(qM,DtauDqM,f);
    else spatialEulerX1Y3Z2DSTDqf<T,VecCM,MatTM,Vec6T>(qM,DtauDqM,f);
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    VecCM qM(&(q.coeffRef(J._offDOF+0)),2);
    spatialExpWYZDSTDqf<T,VecCM,MatTM,Vec6T>(qM,DtauDqM,f);
  }
}
template <typename T>
void JointFunc<T>::DSDqf(const Joint& J,const VecCM q,MatTM DtauDq,VecCM f)
{
  if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    VecCM qM(&(q.coeffRef(J._offDOF+0)),3);
    if(J._typeJoint == Joint::ROT_3D_EXP)
      spatialExpWDSDqf<T,VecCM,MatTM,VecCM>(qM,DtauDq,f);
    else spatialEulerX1Y3Z2DSDqf<T,VecCM,MatTM,VecCM>(qM,DtauDq,f);
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    VecCM qM(&(q.coeffRef(J._offDOF+0)),2);
    spatialExpWYZDSDqf<T,VecCM,MatTM,VecCM>(qM,DtauDq,f);
  }
}

//instance
PRJ_BEGIN
template struct JointFunc<double>;
#ifdef ALL_TYPES
template struct JointFunc<__float128>;
template struct JointFunc<mpfr::mpreal>;
#endif
PRJ_END
