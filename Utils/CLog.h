#ifndef CLOG_H
#define CLOG_H

#include <CommonFile/MathBasic.h>

PRJ_BEGIN

template <typename T>
T clog(T d,T* D,T* DD,T d0,T coef)
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

PRJ_END

#endif
