#ifndef CAMERA_MODEL_H
#define CAMERA_MODEL_H

#include "MathBasic.h"

PRJ_BEGIN

//camera model
template<typename T>
bool getProjectionMatrixFrustum
(const typename ScalarUtil<T>::ScalarMat4& prj,
 T& left,T& right,T& bottom,T& top,T& zNear,T& zFar)
{
  if (prj(3,0)!=0.0f || prj(3,1)!=0.0f || prj(3,2)!=-1.0f || prj(3,3)!=0.0f)
    return false;

  // note: near and far must be used inside this method instead of zNear and zFar
  // because zNear and zFar are references and they may point to the same variable.
  T temp_near = prj(2,3) / (prj(2,2)-1.0f);
  T temp_far = prj(2,3) / (1.0f+prj(2,2));

  left = temp_near * (prj(0,2)-1.0f) / prj(0,0);
  right = temp_near * (1.0f+prj(0,2)) / prj(0,0);

  top = temp_near * (1.0f+prj(1,2)) / prj(1,1);
  bottom = temp_near * (prj(1,2)-1.0f) / prj(1,1);

  zNear = temp_near;
  zFar = temp_far;

  return true;
}
template<typename T>
void getViewMatrixFrame
(const typename ScalarUtil<T>::ScalarMat4& m,
 typename ScalarUtil<T>::ScalarVec3& c,
 typename ScalarUtil<T>::ScalarVec3& x,
 typename ScalarUtil<T>::ScalarVec3& y,
 typename ScalarUtil<T>::ScalarVec3& z)
{
  typename ScalarUtil<T>::ScalarMat4 mInv=m.inverse();

  typename ScalarUtil<T>::ScalarVec4 cH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,0.0f,0.0f,1.0f);
  c=(typename ScalarUtil<T>::ScalarVec3(cH.x(),cH.y(),cH.z())*(1.0f/cH.w()));

  typename ScalarUtil<T>::ScalarVec4 xH=mInv*typename ScalarUtil<T>::ScalarVec4(1.0f,0.0f,0.0f,1.0f);
  x=(typename ScalarUtil<T>::ScalarVec3(xH.x(),xH.y(),xH.z())*(1.0f/xH.w())-c).normalized();

  typename ScalarUtil<T>::ScalarVec4 yH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,1.0f,0.0f,1.0f);
  y=(typename ScalarUtil<T>::ScalarVec3(yH.x(),yH.y(),yH.z())*(1.0f/yH.w())-c).normalized();

  typename ScalarUtil<T>::ScalarVec4 zH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,0.0f,-1.0f,1.0f);
  z=(typename ScalarUtil<T>::ScalarVec3(zH.x(),zH.y(),zH.z())*(1.0f/zH.w())-c).normalized();
}
template<typename T>
typename ScalarUtil<T>::ScalarMat4 getProjection(T fov,T aspect,T near,T far)
{
  T D2R = M_PI/180.0;
  T yScale = 1.0/std::tan(D2R*fov/2);
  T xScale = yScale/aspect;
  T nearmfar = near-far;
  T m[] = {
    xScale, 0,      0,                    0,
    0,      yScale, 0,                    0,
    0,      0,      (far+near)/nearmfar, -1,
    0,      0,      2*far*near/nearmfar,  0
  };
  return Eigen::Map<typename ScalarUtil<T>::ScalarMat4>(m);
}
template<typename T>
typename ScalarUtil<T>::ScalarMat4 lookAt
(const typename ScalarUtil<T>::ScalarVec3& eye,
 const typename ScalarUtil<T>::ScalarVec3& center,
 const typename ScalarUtil<T>::ScalarVec3& up)
{
  typename ScalarUtil<T>::ScalarMat4 m;
  m.setIdentity();
  m.template block<3,1>(0,3)=eye;
  m.template block<3,1>(0,0)=up.cross(eye-center).normalized();
  m.template block<3,1>(0,2)=(eye-center).normalized();
  m.template block<3,1>(0,1)=m.template block<3,1>(0,2).cross(m.template block<3,1>(0,0));
  return m.inverse();
}
template<typename T>
bool getViewMatrixFrame(const typename ScalarUtil<T>::ScalarMat4& m,
                        const typename ScalarUtil<T>::ScalarMat4& prj,
                        typename ScalarUtil<T>::ScalarVec3& c,
                        typename ScalarUtil<T>::ScalarVec3& x,
                        typename ScalarUtil<T>::ScalarVec3& y,
                        typename ScalarUtil<T>::ScalarVec3& z)
{
  T left,right,bottom,top,zNear,zFar;
  if(!getProjectionMatrixFrustum<T>(prj,left,right,bottom,top,zNear,zFar))
    return false;

  typename ScalarUtil<T>::ScalarMat4 mInv=m.inverse();
  typename ScalarUtil<T>::ScalarMat3 mRotInv=(m.block(0,0,3,3)).inverse();

  typename ScalarUtil<T>::ScalarVec4 cH=mInv*typename ScalarUtil<T>::ScalarVec4(0.0f,0.0f,0.0f,1.0f);
  c=(typename ScalarUtil<T>::ScalarVec3(cH.x(),cH.y(),cH.z())*(1.0f/cH.w()));

  x=mRotInv*(typename ScalarUtil<T>::ScalarVec3(1.0f,0.0f,0.0f));
  x=x.normalized()*right;

  y=mRotInv*(typename ScalarUtil<T>::ScalarVec3(0.0f,1.0f,0.0f));
  y=y.normalized()*top;

  z=mRotInv*(typename ScalarUtil<T>::ScalarVec3(0.0f,0.0f,-1.0f));
  z=z.normalized()*zNear;

  return true;
}
template<typename T>
typename ScalarUtil<T>::ScalarVec3
transformHomo(const typename ScalarUtil<T>::ScalarMat4& m,
              const typename ScalarUtil<T>::ScalarVec3& p)
{
  return m.template block<3,3>(0,0)*p+m.template block<3,1>(0,3);
}
template<typename T>
typename ScalarUtil<T>::ScalarVec3
transformHomo(const typename ScalarUtil<T>::ScalarMat3X4& m,
              const typename ScalarUtil<T>::ScalarVec3& p)
{
  return m.template block<3,3>(0,0)*p+m.template block<3,1>(0,3);
}
template<typename T>
typename ScalarUtil<T>::ScalarMat3X4
transformMul(const typename ScalarUtil<T>::ScalarMat3X4& m,
             const typename ScalarUtil<T>::ScalarMat3X4& m2)
{
  typename ScalarUtil<T>::ScalarMat3X4 ret;
  ret.template block<3,3>(0,0)=m.template block<3,3>(0,0)*m2.template block<3,3>(0,0);
  ret.template block<3,1>(0,3)=m.template block<3,3>(0,0)*m2.template block<3,1>(0,3)+m.template block<3,1>(0,3);
  return ret;
}
template<typename T>
typename ScalarUtil<T>::ScalarVec3
transformHomoInv(const typename ScalarUtil<T>::ScalarMat4& m,
                 const typename ScalarUtil<T>::ScalarVec3& p)
{
  return m.template block<3,3>(0,0).inverse()*(p-m.template block<3,1>(0,3));
}
template<typename T>
typename ScalarUtil<T>::ScalarVec3
transformHomoInv(const typename ScalarUtil<T>::ScalarMat3X4& m,
                 const typename ScalarUtil<T>::ScalarVec3& p)
{
  return m.template block<3,3>(0,0).inverse()*(p-m.template block<3,1>(0,3));
}
template<typename T>
typename ScalarUtil<T>::ScalarVec3
transformHomoInvT(const typename ScalarUtil<T>::ScalarMat4& m,
                  const typename ScalarUtil<T>::ScalarVec3& p)
{
  return m.template block<3,3>(0,0).transpose()*(p-m.template block<3,1>(0,3));
}
template<typename T>
typename ScalarUtil<T>::ScalarVec3
transformHomoInvT(const typename ScalarUtil<T>::ScalarMat3X4& m,
                  const typename ScalarUtil<T>::ScalarVec3& p)
{
  return m.template block<3,3>(0,0).transpose()*(p-m.template block<3,1>(0,3));
}
template<typename T>
typename ScalarUtil<T>::ScalarMat4
getOrtho2D(const typename ScalarUtil<T>::ScalarVec3& minC,
           const typename ScalarUtil<T>::ScalarVec3& maxC)
{
  typename ScalarUtil<T>::ScalarMat4 mt=typename ScalarUtil<T>::ScalarMat4::Zero();
  mt(0,0)=2.0f/(maxC.x()-minC.x());
  mt(0,3)=(minC.x()+maxC.x())/(minC.x()-maxC.x());

  mt(1,1)=2.0f/(maxC.y()-minC.y());
  mt(1,3)=(minC.y()+maxC.y())/(minC.y()-maxC.y());

  mt(2,2)=1.0f;
  mt(3,3)=1.0f;
  return mt;
}

PRJ_END

#endif
