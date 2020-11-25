#ifndef ROTATION_UTIL_H
#define ROTATOIN_UTIL_H

#include <CommonFile/MathBasic.h>
#ifdef B2
#undef B2
#endif
#include <Eigen/Eigenvalues>
#include <iostream>
#include <math.h>

PRJ_BEGIN

template <typename T>
struct SinCosTraits;
template <>
struct SinCosTraits<scalarD>
{
  static FORCE_INLINE EIGEN_DEVICE_FUNC void sincos(scalarD theta,scalarD* a,scalarD* b)
  {
    ::sincos(theta,a,b);
  }
};
template <>
struct SinCosTraits<scalarF>
{
  static FORCE_INLINE EIGEN_DEVICE_FUNC void sincos(scalarF theta,scalarF* a,scalarF* b)
  {
    ::sincosf(theta,a,b);
  }
};
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 cross(const typename ScalarUtil<T>::ScalarVec3& v)
{
  typename ScalarUtil<T>::ScalarMat3 ret;
  ret <<
      0,-v[2],v[1],
      v[2],0,-v[0],
      -v[1],v[0],0;
  return ret;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarVec3 invCross(const typename ScalarUtil<T>::ScalarMat3& wCross)
{
  return typename ScalarUtil<T>::ScalarVec3(wCross(2,1),wCross(0,2),wCross(1,0));
}
template <typename T>   //trace([w]*m) = w.dot(invCrossMatTrace(m))
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarVec3 invCrossMatTrace(const typename ScalarUtil<T>::ScalarMat3& m)
{
  return typename ScalarUtil<T>::ScalarVec3(m(1,2)-m(2,1),m(2,0)-m(0,2),m(0,1)-m(1,0));
}
template <typename T>   //trace([wA]*m*[wB]) = wA.dot(invDoubleCrossMatTrace(m)*wB)
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 invDoubleCrossMatTrace(const typename ScalarUtil<T>::ScalarMat3& m)
{
  typename ScalarUtil<T>::ScalarMat3 ret;
  ret <<
      -m(1,1)-m(2,2),m(1,0),m(2,0),
      m(0,1),-m(0,0)-m(2,2),m(2,1),
      m(0,2),m(1,2),-m(0,0)-m(1,1);
  return ret;
}
template <typename T>   //trace([l]*[wA]*m*[wB]) = wA.dot(invDoubleCrossMatTrace(l,m)*wB)
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 invDoubleCrossMatTrace(const typename ScalarUtil<T>::ScalarVec3& l,const typename ScalarUtil<T>::ScalarMat3& m)
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
template <typename T,typename T2>   //hinge joint
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 expWZ(T z,T2* DRDZ)
{
  T sz,cz;
  sincos(z,&sz,&cz);

  typename ScalarUtil<T>::ScalarMat3 ret;
  ret <<
      cz,-sz,0,
      sz, cz,0,
      0 , 0 ,1;

  if(DRDZ)
    *DRDZ <<   0, 0,1;
  return ret;
}
template <typename T,typename T2>   //ball joint
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 expWYZ(T y,T z,T2* DRDY,T2* DRDZ,T2* DRDYDZ)
{
  T sy,cy,sz,cz;
  SinCosTraits<T>::sincos(y,&sy,&cy);
  SinCosTraits<T>::sincos(z,&sz,&cz);

  typename ScalarUtil<T>::ScalarMat3 ret;
  ret <<
      cy*cz,-sz,sy*cz,
      cy*sz, cz,sy*sz,
      -sy,  0,   cy;

  if(DRDY)
    *DRDY << -sz,cz,0;
  if(DRDZ)
    *DRDZ <<   0, 0,1;

  if(DRDYDZ)
    *DRDYDZ << -cz,-sz,0;
  return ret;
}
template <typename T,typename T2>
FORCE_INLINE EIGEN_DEVICE_FUNC void expWYZLambda(T z,T lambdaZ,T2& DRDYDZ)
{
  T sz,cz;
  sincos(z,&sz,&cz);
  DRDYDZ << sz*lambdaZ,-cz*lambdaZ,0;
}
template <typename T,typename T2>   //rotation 3D euler angle
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 eulerX1Y3Z2(const typename ScalarUtil<T>::ScalarVec3& w,T2 diffV[3],T2 ddiffV[9])
{
  //this is a special type of Euler-Angle in favor of our definition in RootInvariantMLP
  T sx,cx,sy,cy,sz,cz;
  sincos(w[0],&sx,&cx);
  sincos(w[1],&sy,&cy);
  sincos(w[2],&sz,&cz);

  typename ScalarUtil<T>::ScalarMat3 ret;
  ret(0,0)=cy*cz;
  ret(0,1)=sx*sy-cx*cy*sz;
  ret(0,2)=sx*cy*sz+cx*sy;
  ret(1,0)=sz;
  ret(1,1)=cx*cz;
  ret(1,2)=-sx*cz;
  ret(2,0)=-sy*cz;
  ret(2,1)=cx*sy*sz+sx*cy;
  ret(2,2)=cx*cy-sx*sy*sz;

  if(diffV) {
    diffV[0] << cy*cz,sz,-sy*cz;
    diffV[1] <<     0, 1,     0;
    diffV[2] <<    sy, 0,    cy;
  }

  if(ddiffV) {
    ddiffV[0] <<      0, 0,     0;
    ddiffV[1] << -sy*cz, 0,-cy*cz;
    ddiffV[2] << -cy*sz,cz, sy*sz;

    ddiffV[3] << 0,0,0;
    ddiffV[4] << 0,0,0;
    ddiffV[5] << 0,0,0;

    ddiffV[6] <<  0,0,  0;
    ddiffV[7] << cy,0,-sy;
    ddiffV[8] <<  0,0,  0;
  }
  return ret;
}
template <typename T,typename T2>
FORCE_INLINE EIGEN_DEVICE_FUNC void eulerX1Y3Z2Lambda(const typename ScalarUtil<T>::ScalarVec3& w,const typename ScalarUtil<T>::ScalarVec3& lambda,T2 ddiffV[9])
{
  T sx,cx,sy,cy,sz,cz;
  sincos(w[0],&sx,&cx);
  sincos(w[1],&sy,&cy);
  sincos(w[2],&sz,&cz);

  ddiffV[0] <<                                0,            0,                              0;
  ddiffV[1] << -cy*cz*lambda[1]+sy*sz*lambda[2],            0,sy*cz*lambda[1]+cy*sz*lambda[2];
  ddiffV[2] <<  sy*sz*lambda[1]-cy*cz*lambda[2],-sz*lambda[2],cy*sz*lambda[1]+sy*cz*lambda[2];

  ddiffV[3] << 0,0,0;
  ddiffV[4] << 0,0,0;
  ddiffV[5] << 0,0,0;

  ddiffV[6] <<             0,0,            0;
  ddiffV[7] << -sy*lambda[1],0,-cy*lambda[1];
  ddiffV[8] <<             0,0,            0;
}
template <typename T>   //rotation 3D rodrigues, helper for diffV
FORCE_INLINE EIGEN_DEVICE_FUNC T f_tMst_t3(T t,T st,T ct,T eps)
{
  T t2=t*t,t3=t*t2;
  if(std::abs(t)<eps)
    return (1037836800+t2*(t2*(1235520+t2*(t2*(156-t2)-17160))-51891840))/6227020800;
  else return (t-st)/t3;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T f_1Mct_t2(T t,T st,T ct,T eps)
{
  T t2=t*t;
  if(std::abs(t)<eps)
    return (239500800+t2*(t2*(665280+t2*(t2*(132-t2)-11880))-19958400))/479001600;
  else return (1-ct)/t2;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T f_st_t(T t,T st,T ct,T eps)
{
  T t2=t*t;
  if(std::abs(t)<eps)
    return (39916800+t2*(t2*(332640+t2*(t2*(110-t2)-7920))-6652800))/39916800;
  else return st/t;
}
template <typename T>   //helper for ddiffV
FORCE_INLINE EIGEN_DEVICE_FUNC T f_3stM2tMtct_t5(T t,T st,T ct,T eps)
{
  T t2=t*t,t5=t2*t2*t;
  if(std::abs(t)<eps)
    return (t2*(t2*(t2*(t2*(t2-175)+21840)-1801800)+86486400)-1816214400)/108972864000;
  else return (3*st-2*t-t*ct)/t5;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T f_tstP2ctM2_t4(T t,T st,T ct,T eps)
{
  T t2=t*t,t4=t2*t2;
  if(std::abs(t)<eps)
    return (t2*(t2*(t2*(t2*(3*t2-455)+48048)-3243240)+121080960)-1816214400)/21794572800;
  else return (-2+2*ct+t*st)/t4;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T f_stMtct_t3(T t,T st,T ct,T eps)
{
  T t2=t*t,t3=t*t2;
  if(std::abs(t)<eps)
    return (172972800+t2*(t2*(617760+t2*(t2*(130-t2)-11440))-17297280))/518918400;
  else return (st-t*ct)/t3;
}
template <typename T>   //helper for dddiffVLambda
FORCE_INLINE EIGEN_DEVICE_FUNC T fD_tMst_t3(T t,T st,T ct,T eps)
{
  T t2=t*t,t5=t2*t2*t;
  if(std::abs(t)<eps)
    return (t2*(t2*(t2*(t2*(t2-175)+21840)-1801800)+86486400)-1816214400)/108972864000;
  else return (-2*t-t*ct+3*st)/t5;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T fD_1Mct_t2(T t,T st,T ct,T eps)
{
  T t2=t*t,t4=t2*t2;
  if(std::abs(t)<eps)
    return (t2*(t2*(t2*(t2*(3*t2-455)+48048)-3243240)+121080960)-1816214400)/21794572800;
  else return (-2+2*ct+t*st)/t4;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T fD_3stM2tMtct_t5(T t,T st,T ct,T eps)
{
  T t2=t*t,t7=t2*t2*t2*t;
  if(std::abs(t)<eps)
    return (23524300800+t2*(t2*(17821440+t2*(t2*(1360-7*t2)-190400))-980179200))/14820309504000;
  else return (8*t+7*t*ct+(t2-15)*st)/t7;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T fD_tstP2ctM2_t4(T t,T st,T ct,T eps)
{
  T t2=t*t,t6=t2*t2*t2;
  if(std::abs(t)<eps)
    return (9686476800+t2*(t2*(11531520+t2*(t2*(1200-7*t2)-145600))-518918400))/871782912000;
  else return -(-8+(8-t2)*ct+5*t*st)/t6;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T fD_stMtct_t3(T t,T st,T ct,T eps)
{
  T t2=t*t,t5=t2*t2*t;
  if(std::abs(t)<eps)
    return (t2*(t2*(t2*(t2*(t2-150)+15600)-1029600)+37065600)-518918400)/7783776000;
  else return (3*t*ct+(t2-3)*st)/t5;
}
template <typename T,typename T2>
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 expWGradV(const typename ScalarUtil<T>::ScalarVec3& w,T2* diffV=NULL,T2* ddiffV=NULL)
{
  typedef typename ScalarUtil<T>::ScalarMat3 MAT3;
  typedef typename ScalarUtil<T>::ScalarVec3 VEC3;
  T epsT=1E-1f;
  T t=w.norm(),st,ct;
  SinCosTraits<T>::sincos(t,&st,&ct);
  MAT3 cw=cross<T>(w);
  //W
  T _tMst_t3=f_tMst_t3<T>(t,st,ct,epsT);
  T _1Mct_t2=f_1Mct_t2<T>(t,st,ct,epsT);
  T _st_t=f_st_t<T>(t,st,ct,epsT);
  if(diffV) {
#define DIFFV(I) diffV[I]=_tMst_t3*w[I]*w+_1Mct_t2*w.cross(VEC3::Unit(I))+_st_t*VEC3::Unit(I);
    DIFFV(0)
    DIFFV(1)
    DIFFV(2)
#undef DIFFV
  }
  //DW
  if(ddiffV) {
    T _3stM2tMtct_t5=f_3stM2tMtct_t5<T>(t,st,ct,epsT);
    T _tstP2ctM2_t4=f_tstP2ctM2_t4<T>(t,st,ct,epsT);
    T _stMtct_t3=f_stMtct_t3<T>(t,st,ct,epsT);
#define DDIFFV(I,J) ddiffV[I*3+J]=  \
    _tMst_t3*(VEC3::Unit(J)*w[I]+(I==J?1:0)*w)+ \
    _3stM2tMtct_t5*(w*w[I]*w[J])+   \
    _1Mct_t2*(VEC3::Unit(J)).cross(VEC3::Unit(I))+  \
    _tstP2ctM2_t4*w.cross(VEC3::Unit(I))*w[J]-  \
    _stMtct_t3*VEC3::Unit(I)*w[J];
    DDIFFV(0,0)
    DDIFFV(0,1)
    DDIFFV(0,2)
    DDIFFV(1,0)
    DDIFFV(1,1)
    DDIFFV(1,2)
    DDIFFV(2,0)
    DDIFFV(2,1)
    DDIFFV(2,2)
#undef DDIFFV
  }
  //R
  return MAT3::Identity()+cw*_st_t+cw*cw*_1Mct_t2;
}
template <typename T,typename T2>
FORCE_INLINE EIGEN_DEVICE_FUNC void expWGradVLambda(const typename ScalarUtil<T>::ScalarVec3& w,const typename ScalarUtil<T>::ScalarVec3& lambda,T2* ddiffV)
{
  typedef typename ScalarUtil<T>::ScalarVec3 VEC3;
  T epsT=1E-1f;
  T t=w.norm(),st,ct,dwl=w.dot(lambda);
  sincos(t,&st,&ct);
  //W
  T _tMst_t3=f_tMst_t3<T>(t,st,ct,epsT);
  T _3stM2tMtct_t5=f_3stM2tMtct_t5<T>(t,st,ct,epsT);
  T _tstP2ctM2_t4=f_tstP2ctM2_t4<T>(t,st,ct,epsT);
  T _stMtct_t3=f_stMtct_t3<T>(t,st,ct,epsT);

  T D_tMst_t3=fD_tMst_t3<T>(t,st,ct,epsT)*dwl;
  T D_1Mct_t2=fD_1Mct_t2<T>(t,st,ct,epsT)*dwl;
  T D_3stM2tMtct_t5=fD_3stM2tMtct_t5<T>(t,st,ct,epsT)*dwl;
  T D_tstP2ctM2_t4=fD_tstP2ctM2_t4<T>(t,st,ct,epsT)*dwl;
  T D_stMtct_t3=fD_stMtct_t3<T>(t,st,ct,epsT)*dwl;
#define DDIFFV(I,J) ddiffV[I*3+J]=  \
  _tMst_t3*(VEC3::Unit(J)*lambda[I]+(I==J?1:0)*lambda)+ \
  _3stM2tMtct_t5*(lambda*w[I]*w[J]+w*lambda[I]*w[J]+w*w[I]*lambda[J])+   \
  _tstP2ctM2_t4*(lambda.cross(VEC3::Unit(I))*w[J]+w.cross(VEC3::Unit(I))*lambda[J])-  \
  _stMtct_t3*VEC3::Unit(I)*lambda[J]+   \
    \
  D_tMst_t3*(VEC3::Unit(J)*w[I]+(I==J?1:0)*w)+ \
  D_3stM2tMtct_t5*(w*w[I]*w[J])+   \
  D_1Mct_t2*(VEC3::Unit(J)).cross(VEC3::Unit(I))+  \
  D_tstP2ctM2_t4*w.cross(VEC3::Unit(I))*w[J]-  \
  D_stMtct_t3*VEC3::Unit(I)*w[J];
  DDIFFV(0,0)
  DDIFFV(0,1)
  DDIFFV(0,2)
  DDIFFV(1,0)
  DDIFFV(1,1)
  DDIFFV(1,2)
  DDIFFV(2,0)
  DDIFFV(2,1)
  DDIFFV(2,2)
#undef DDIFFV
}
template <typename T,typename T2>
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarMat3 expWZY(T z,T y,T2* DRDZ,T2* DRDY,T2* DRDZDY)
{
  T sz=sin(z);
  T cz=cos(z);
  T sy=sin(y);
  T cy=cos(y);

  typename ScalarUtil<T>::ScalarMat3 ret;
  ret <<
      cy*cz,-cy*sz,sy,
      sz,    cz, 0,
      -sy*cz, sy*sz,cy;
  if(DRDZ)
    *DRDZ << sy,0,cy;
  if(DRDY)
    *DRDY <<  0,1, 0;
  if(DRDZDY)
    *DRDZDY << cy,0,-sy;
  return ret;
}
template <typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC typename ScalarUtil<T>::ScalarVec3 invExpW(const typename ScalarUtil<T>::ScalarMat3& R)
{
  T cosTheta=(R.trace()-1)/2,coef;
  if(cosTheta < ScalarUtil<T>::scalar_eps()-1) {
    sizeType a,b,c;
    R.diagonal().maxCoeff(&a);
    b=(a+1)%3;
    c=(a+2)%3;

    T s=std::sqrt(R(a,a)-R(b,b)-R(c,c)+1);
    typename ScalarUtil<T>::ScalarVec3 ret;
    ret[a]=s/2;
    ret[b]=(R(a,b)+R(b,a))/(2*s);
    ret[c]=(R(a,c)+R(c,a))/(2*s);
    return ret.normalized()*M_PI;
  } else if(cosTheta > 1-ScalarUtil<T>::scalar_eps()) {
    T theta=acos(std::min<T>(cosTheta,1));
    coef=0.5f+theta*theta/12.0f;  //using taylor expansion
  } else {
    T theta=acos(cosTheta),sinTheta=sin(theta);
    coef=theta/(2*sinTheta);
  }
  return typename ScalarUtil<T>::ScalarVec3(R(2,1)-R(1,2),R(0,2)-R(2,0),R(1,0)-R(0,1))*coef;
}
template<typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T getAngle3D(const typename ScalarUtil<T>::ScalarVec3& a,const typename ScalarUtil<T>::ScalarVec3& b)
{
  T denom=std::max<T>(a.norm()*b.norm(),ScalarUtil<T>::scalar_eps());
  return acos(std::max<T>(std::min<T>(a.dot(b)/denom,(T)(1.0f-ScalarUtil<T>::scalar_eps())),(T)(-1.0f+ScalarUtil<T>::scalar_eps())));
}
template<typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T getAngle2D(const typename ScalarUtil<T>::ScalarVec2& a,const typename ScalarUtil<T>::ScalarVec2& b)
{
  T denom=std::max<T>(a.norm()*b.norm(),ScalarUtil<T>::scalar_eps());
  return acos(std::max<T>(std::min<T>(a.dot(b)/denom,(T)(1.0f-ScalarUtil<T>::scalar_eps())),(T)(-1.0f+ScalarUtil<T>::scalar_eps())));
}

PRJ_END

#endif
