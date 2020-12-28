#include "MultiPrecisionSeparatingPlane.h"
#include <Utils/RotationUtil.h>

USE_PRJ_NAMESPACE

template <typename T>
MultiPrecisionSeparatingPlane<T>::MultiPrecisionSeparatingPlane(Options& ops,Vec4T& plane,T d0)
  :MultiPrecisionLQP<T>(ops),_plane(plane),_d0(d0)
{
  MultiPrecisionLQP<T>::_H.setZero(1,1);
  MultiPrecisionLQP<T>::_c.setZero(4);
}
template <typename T>
void MultiPrecisionSeparatingPlane<T>::clearPoints()
{
  _pssL.clear();
  _pssR.clear();
}
template <typename T>
void MultiPrecisionSeparatingPlane<T>::resetPoints(const PSS& pss)
{
  _plane.template segment<3>(0)/=std::sqrt(_plane.template segment<3>(0).squaredNorm());
  for(const Vec3T& p:pss)
    if(_plane.template segment<3>(0).dot(p)+_plane[3]<0)
      _pssR.push_back(p);
    else _pssL.push_back(p);
}
template <typename T>
void MultiPrecisionSeparatingPlane<T>::resetPoints(const Mat3XT& pssL,const Mat3XT& pssR)
{
  if(_plane.template segment<3>(0).dot(pssL.col(0))+_plane[3]<0) {
    for(sizeType i=0; i<pssL.cols(); i++)
      _pssR.push_back(pssL.col(i));
    for(sizeType i=0; i<pssR.cols(); i++)
      _pssL.push_back(pssR.col(i));
  } else {
    for(sizeType i=0; i<pssR.cols(); i++)
      _pssR.push_back(pssR.col(i));
    for(sizeType i=0; i<pssL.cols(); i++)
      _pssL.push_back(pssL.col(i));
  }
}
template <typename T>
T MultiPrecisionSeparatingPlane<T>::computeFGH(T mu,const Vec& w,Vec* g,DMat* h) const
{
  T ret=0;
  Vec4T dE;
  Vec3T dR[3],ddR[9],dN;
  if(g)
    g->setZero(4);
  if(h)
    h->setZero(4,4);
  Vec3T n=eulerX1Y3Z2<T,Vec3T>(w.template segment<3>(0),g?dR:NULL,h?ddR:NULL)*_plane.template segment<3>(0);
  for(sizeType pass=0; pass<2; pass++)
    for(const Vec3T& pL:(pass==0?_pssL:_pssR)) {
      T sgn=pass==0?1:-1;
      T tmp=(pL.dot(n)+w[3]+_plane[3])*sgn;
      T D,DD,E=clog(tmp,&D,&DD,_d0,mu);
      if(!std::isfinite(E))
        return E;
      else if(E==0)
        continue;
      ret+=E;
      if(g) {
        dE=Vec4T(dR[0].cross(n).dot(pL),dR[1].cross(n).dot(pL),dR[2].cross(n).dot(pL),1)*sgn;
        *g+=D*dE;
      }
      if(h) {
        *h+=DD*dE*dE.transpose();
        //we avoid this block to make positive definite
        //dN=D*pL*sgn;
        //for(sizeType r=0; r<3; r++)
        //  for(sizeType c=0; c<3; c++)
        //    h->coeffRef(r,c)+=((cross<T>(dR[r])*cross<T>(dR[c])+cross<T>(ddR[r*3+c]))*n).dot(dN);
      }
    }
  return ret;
}
template <typename T>
typename MultiPrecisionSeparatingPlane<T>::Vec MultiPrecisionSeparatingPlane<T>::sampleValidW() const
{
  Vec w=Vec4T::Random();
  while(true)
    if(std::isfinite(computeFGH(1,w)))
      break;
    else w*=0.5f;
  return w;
}
template <typename T>
typename MultiPrecisionSeparatingPlane<T>::Vec4T MultiPrecisionSeparatingPlane<T>::solve(bool& succ)
{
  Vec ret=MultiPrecisionLQP<T>::solve(succ);
  Vec3T n=eulerX1Y3Z2<T,Vec3T>(ret.template segment<3>(0),NULL,NULL)*_plane.template segment<3>(0);
  n/=std::sqrt(n.squaredNorm());
  return _plane=Vec4T(n[0],n[1],n[2],_plane[3]+ret[3]);
}
template <typename T>
T MultiPrecisionSeparatingPlane<T>::clog(T d,T* D,T* DD,T d0,T coef)
{
  if(d<=0)
    return ScalarUtil<T>::scalar_nanq();
  else if(d>d0) {
    if(D)
      *D=0;
    if(DD)
      *DD=0;
    return 0;
  }
  T valLog=std::log(d/d0);
  T valLogC=valLog*(d-d0);
  T relD=(d-d0)/d;
  if(D)
    *D=-(2*valLogC+(d-d0)*relD)*coef;
  if(DD)
    *DD=-(4*relD-relD*relD+2*valLog)*coef;
  return -valLogC*(d-d0)*coef;
}
template <typename T>
void MultiPrecisionSeparatingPlane<T>::initialGuess(Vec& w) const
{
  w.setZero(4);
}
template <typename T>
void MultiPrecisionSeparatingPlane<T>::limitAlpha(T& alpha,const Vec&,const Vec&) const
{
  alpha=std::min<T>(alpha,1);
}
//instance
PRJ_BEGIN
template class MultiPrecisionSeparatingPlane<double>;
#ifdef ALL_TYPES
template class MultiPrecisionSeparatingPlane<__float128>;
template class MultiPrecisionSeparatingPlane<mpfr::mpreal>;
#endif
PRJ_END
