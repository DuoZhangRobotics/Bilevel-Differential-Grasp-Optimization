#ifndef CROSS_SPATIAL_UTIL_H
#define CROSS_SPATIAL_UTIL_H

#include <CommonFile/MathBasic.h>
#include <Utils/DebugGradient.h>
#ifdef B2
#undef B2
#endif
#include <Utils/Scalar.h>
#include <iostream>

PRJ_BEGIN

//cross
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat3 cross(const typename ScalarUtil<T>::ScalarVec3& v)
{
  typename ScalarUtil<T>::ScalarMat3 ret;
  ret <<
      0,-v[2],v[1],
      v[2],0,-v[0],
      -v[1],v[0],0;
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 invCross(const typename ScalarUtil<T>::ScalarMat3& wCross)
{
  return typename ScalarUtil<T>::ScalarVec3(wCross(2,1),wCross(0,2),wCross(1,0));
}
template <typename T>   //trace([w]*m) = w.dot(invCrossMatTrace(m))
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 invCrossMatTrace(const typename ScalarUtil<T>::ScalarMat3& m)
{
  return typename ScalarUtil<T>::ScalarVec3(m(1,2)-m(2,1),m(2,0)-m(0,2),m(0,1)-m(1,0));
}
template <typename T>   //trace([wA]*m*[wB]) = wA.dot(invDoubleCrossMatTrace(m)*wB)
FORCE_INLINE typename ScalarUtil<T>::ScalarMat3 invDoubleCrossMatTrace(const typename ScalarUtil<T>::ScalarMat3& m)
{
  typename ScalarUtil<T>::ScalarMat3 ret;
  ret <<
      -m(1,1)-m(2,2),m(1,0),m(2,0),
      m(0,1),-m(0,0)-m(2,2),m(2,1),
      m(0,2),m(1,2),-m(0,0)-m(1,1);
  return ret;
}
template <typename T>   //trace([l]*[wA]*m*[wB]) = wA.dot(invDoubleCrossMatTrace(l,m)*wB)
FORCE_INLINE typename ScalarUtil<T>::ScalarMat3 invDoubleCrossMatTrace(const typename ScalarUtil<T>::ScalarVec3& l,const typename ScalarUtil<T>::ScalarMat3& m)
{
  typename ScalarUtil<T>::ScalarMat3 ret;
  ret<<
     l[0]*m(2,1)-l[0]*m(1,2),
     -l[2]*m(2,2)-l[0]*m(2,0)-l[1]*m(1,2),
     l[2]*m(2,1)+l[1]*m(1,1)+l[0]*m(1,0),

     l[2]*m(2,2)+l[1]*m(2,1)+l[0]*m(0,2),
     m(0,2)*l[1]-l[1]*m(2,0),
     -l[2]*m(2,0)-m(0,1)*l[1]-l[0]*m(0,0),

     -m(1,2)*l[2]-l[1]*m(1,1)-l[0]*m(0,1),
     m(0,2)*l[2]+l[1]*m(1,0)+l[0]*m(0,0),
     m(1,0)*l[2]-m(0,1)*l[2];
  return ret;
}
template <typename T>
FORCE_INLINE void debugCross()
{
  typedef typename ScalarUtil<T>::ScalarVec3 Vec3T;
  typedef typename ScalarUtil<T>::ScalarMat3 Mat3T;
  DEFINE_NUMERIC_DELTA_T(T)
  INFO("-------------------------------------------------------------DebugCross")
  Vec3T l=Vec3T::Random();
  Vec3T w=Vec3T::Random();
  Vec3T wA=Vec3T::Random();
  Vec3T wB=Vec3T::Random();
  Mat3T m=Mat3T::Random();
  DEBUG_GRADIENT("cross",std::sqrt(w.squaredNorm()),std::sqrt((w-invCross<T>(cross<T>(w))).squaredNorm()))
  DEBUG_GRADIENT("invCrossMatTrace",std::sqrt((cross<T>(w)*m).trace()),std::sqrt((cross<T>(w)*m).trace()-w.dot(invCrossMatTrace<T>(m))))
  DEBUG_GRADIENT("invDoubleCrossMatTrace",(cross<T>(wA)*m*cross<T>(wB)).trace(),(cross<T>(wA)*m*cross<T>(wB)).trace()-wA.dot(invDoubleCrossMatTrace<T>(m)*wB))
  DEBUG_GRADIENT("invDoubleCrossMatTrace2",(cross<T>(l)*cross<T>(wA)*m*cross<T>(wB)).trace(),(cross<T>(l)*cross<T>(wA)*m*cross<T>(wB)).trace()-wA.dot(invDoubleCrossMatTrace<T>(l,m)*wB))
}
//spatial
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat6 toSpatial(const typename ScalarUtil<T>::ScalarMat3X4 t)
{
  typename ScalarUtil<T>::ScalarMat6 ret;
  ret.setZero();
  ret.template block<3,3>(0,0)=t.template block<3,3>(0,0);
  ret.template block<3,3>(3,3)=t.template block<3,3>(0,0);
  ret.template block<3,3>(3,0)=cross<T>(t.template block<3,1>(0,3))*t.template block<3,3>(0,0);
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat3X4 fromSpatial(const typename ScalarUtil<T>::ScalarMat6 t)
{
  typename ScalarUtil<T>::ScalarMat3X4 ret;
  ret.template block<3,3>(0,0)=t.template block<3,3>(0,0);
  ret.template block<3,1>(0,3)=invCross<T>(t.template block<3,3>(3,0)*t.template block<3,3>(0,0).transpose());
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat6 spatialCross(const typename ScalarUtil<T>::ScalarVec6& v)
{
  typename ScalarUtil<T>::ScalarMat6 ret;
  ret.setZero();
  ret.template block<3,3>(3,0)=cross<T>(v.template segment<3>(3));
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=cross<T>(v.template segment<3>(0));
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec6 spatialCross(const typename ScalarUtil<T>::ScalarVec6& v,const typename ScalarUtil<T>::ScalarVec6& w)
{
  typename ScalarUtil<T>::ScalarVec6 ret;
  ret.template segment<3>(0)=v.template segment<3>(0).cross(w.template segment<3>(0));
  ret.template segment<3>(3)=v.template segment<3>(0).cross(w.template segment<3>(3))+v.template segment<3>(3).cross(w.template segment<3>(0));
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat6 spatialCrossStar(const typename ScalarUtil<T>::ScalarVec6& v)
{
  typename ScalarUtil<T>::ScalarMat6 ret;
  ret.setZero();
  ret.template block<3,3>(0,3)=cross<T>(v.template segment<3>(3));
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=cross<T>(v.template segment<3>(0));
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec6 spatialCrossStar(const typename ScalarUtil<T>::ScalarVec6& v,const typename ScalarUtil<T>::ScalarVec6& w)
{
  typename ScalarUtil<T>::ScalarVec6 ret;
  ret.template segment<3>(0)=v.template segment<3>(0).cross(w.template segment<3>(0))+v.template segment<3>(3).cross(w.template segment<3>(3));
  ret.template segment<3>(3)=v.template segment<3>(0).cross(w.template segment<3>(3));
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat6 spatialXStar(const typename ScalarUtil<T>::ScalarMat6& X)
{
  typename ScalarUtil<T>::ScalarMat6 ret;
  ret.setZero();
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=X.template block<3,3>(0,0);
  ret.template block<3,3>(0,3)=-X.template block<3,3>(0,0)*X.template block<3,3>(3,0).transpose()*X.template block<3,3>(0,0);
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarMat6 spatialInv(const typename ScalarUtil<T>::ScalarMat6& X)
{
  typename ScalarUtil<T>::ScalarMat6 ret;
  ret.setZero();
  ret.template block<3,3>(0,0)=ret.template block<3,3>(3,3)=X.template block<3,3>(0,0).transpose();
  ret.template block<3,3>(3,0)=-X.template block<3,3>(0,0).transpose()*X.template block<3,3>(3,0)*X.template block<3,3>(0,0).transpose();
  return ret;
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 spatialVel(const typename ScalarUtil<T>::ScalarVec6& V,const typename ScalarUtil<T>::ScalarVec3& v)
{
  return V.template segment<3>(0).cross(v)+V.template segment<3>(3);
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 applyTrans(const typename ScalarUtil<T>::ScalarMat3X4& X,const typename ScalarUtil<T>::ScalarVec3& v)
{
  return X.template block<3,3>(0,0)*v+X.template block<3,1>(0,3);
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 applyTransInv(const typename ScalarUtil<T>::ScalarMat3X4& X,const typename ScalarUtil<T>::ScalarVec3& v)
{
  return X.template block<3,3>(0,0).transpose()*(v-X.template block<3,1>(0,3));
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 spatialApplyTrans(const typename ScalarUtil<T>::ScalarMat6& X,const typename ScalarUtil<T>::ScalarVec3& v)
{
  return applyTrans<T>(fromSpatial<T>(X),v);
}
template <typename T>
FORCE_INLINE typename ScalarUtil<T>::ScalarVec3 spatialApplyTransInv(const typename ScalarUtil<T>::ScalarMat6& X,const typename ScalarUtil<T>::ScalarVec3& v)
{
  return applyTrans<T>(fromSpatial<T>(spatialInv<T>(X)),v);
}
template <typename T>
FORCE_INLINE void debugSpatial()
{
  typedef typename ScalarUtil<T>::ScalarMat3X4 Mat3X4T;
  typedef typename ScalarUtil<T>::ScalarVec3 Vec3T;
  typedef typename ScalarUtil<T>::ScalarVec4 Vec4T;
  typedef typename ScalarUtil<T>::ScalarMat6 Mat6T;
  typedef typename ScalarUtil<T>::ScalarVec6 Vec6T;
  typedef typename ScalarUtil<T>::ScalarQuat QuatT;
  DEFINE_NUMERIC_DELTA_T(T)
  INFO("-------------------------------------------------------------DebugSpatial")
  Vec4T q;
  Vec6T a,a2;
  Mat3X4T t,t2;
  q.setRandom();
  q/=std::sqrt(q.squaredNorm());
  a.setRandom();
  a2.setRandom();
  ROT(t)=QuatT(q[0],q[1],q[2],q[3]).toRotationMatrix();
  CTR(t)=Vec3T::Random();
  Mat6T ST=toSpatial<T>(t);
  t2=fromSpatial<T>(ST);
  Mat6T invST=spatialInv<T>(ST);
  Mat6T invTST=spatialXStar<T>(ST);
  DEBUG_GRADIENT("toFromSpatial",std::sqrt(t.squaredNorm()),std::sqrt((t-t2).squaredNorm()))
  DEBUG_GRADIENT("invSpatial",std::sqrt(ST.squaredNorm()),std::sqrt((ST*invST-Mat6T::Identity()).squaredNorm()))
  DEBUG_GRADIENT("spatialXStar",std::sqrt(invST.squaredNorm()),std::sqrt((invST-invTST.transpose()).squaredNorm()))
  DEBUG_GRADIENT("spatialCross",std::sqrt(spatialCross<T>(a,a2).squaredNorm()),std::sqrt((spatialCross<T>(a,a2)-spatialCross<T>(a)*a2).squaredNorm()))
  DEBUG_GRADIENT("spatialCrossStar",std::sqrt(spatialCrossStar<T>(a,a2).squaredNorm()),std::sqrt((spatialCrossStar<T>(a,a2)-spatialCrossStar<T>(a)*a2).squaredNorm()))

}

PRJ_END

#endif
