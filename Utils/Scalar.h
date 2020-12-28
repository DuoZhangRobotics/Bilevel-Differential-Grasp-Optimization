#ifndef SCALAR_H
#define SCALAR_H

#include "ScalarMpfr.h"
#include "ScalarQ.h"

namespace std
{
double to_double(float a);
double to_double(double a);
double to_double(sizeType a);
}

PRJ_BEGIN

template <typename MAT>
typename MAT::Scalar maxCoeff(const MAT& m,sizeType* index=NULL)
{
  typename MAT::Scalar ret;
  for(sizeType r=0; r<m.size(); r++)
    if(r==0 || ret<m.data()[r]) {
      ret=m.data()[r];
      if(index)
        *index=r;
    }
  return ret;
}
template <typename MAT>
typename MAT::Scalar minCoeff(const MAT& m,sizeType* index=NULL)
{
  typename MAT::Scalar ret;
  for(sizeType r=0; r<m.size(); r++)
    if(r==0 || ret>m.data()[r]) {
      ret=m.data()[r];
      if(index)
        *index=r;
    }
  return ret;
}
template <typename MAT>
typename MAT::Scalar LInfNorm(const MAT& m)
{
  typename MAT::Scalar ret=0;
  for(sizeType r=0; r<m.rows(); r++)
    for(sizeType c=0; c<m.cols(); c++)
      ret=std::max(ret,std::abs(m(r,c)));
  return ret;
}
template <typename MAT>
typename MAT::Scalar L1Norm(const MAT& m)
{
  typename MAT::Scalar ret=0;
  for(sizeType r=0; r<m.rows(); r++)
    for(sizeType c=0; c<m.cols(); c++)
      ret+=std::abs(m(r,c));
  return ret;
}

PRJ_END

#endif
