#ifndef SCALAR_Q_H
#define SCALAR_Q_H

#include <iostream>
#include <quadmath.h>
typedef __float128 scalarQ;
namespace std
{
#define STDQUAD1(NAME) scalarQ NAME(scalarQ a);
#define STDQUAD2(NAME) scalarQ NAME(scalarQ a,scalarQ b);
STDQUAD1(acos)
STDQUAD1(acosh)
STDQUAD1(asin)
STDQUAD1(asinh)
STDQUAD1(atan)
STDQUAD1(atanh)
STDQUAD2(atan2)
STDQUAD1(cbrt)
STDQUAD1(ceil)
STDQUAD1(cosh)
STDQUAD1(cos)
STDQUAD1(erf)
STDQUAD1(erfc)
STDQUAD1(exp)
STDQUAD1(fabs)
STDQUAD1(floor)
STDQUAD2(fmax)
STDQUAD2(fmin)
STDQUAD2(fmod)
STDQUAD1(isinf)
STDQUAD1(isnan)
STDQUAD1(round)
STDQUAD1(log)
STDQUAD1(log10)
STDQUAD1(log2)
STDQUAD2(pow)
STDQUAD1(sinh)
STDQUAD1(sin)
STDQUAD1(sqrt)
STDQUAD1(tanh)
STDQUAD1(tan)
#undef STDQUAD1
#undef STDQUAD2
scalarQ abs(scalarQ a);
bool isfinite(scalarQ a);
scalarQ frexp(scalarQ a,int* exp);
scalarQ ldexp(scalarQ a,int exp);
double to_double(scalarQ a);
void convert_scalar(double a,double& to);
void convert_scalar(double a,__float128& to);
void convert_scalar(__float128 a,double& to);
void convert_scalar(__float128 a,__float128& to);
std::string to_string(scalarQ a);
istream& operator>>(istream& input,scalarQ& x);
ostream& operator<<(ostream& output,scalarQ x);
}

#include <CommonFile/MathBasic.h>

PRJ_BEGIN

#include <CommonFile/BeginAllEigen.h>
#undef NAME_EIGEN_ALLTYPES
#define NAME_EIGEN_ALLTYPES(NAME,size1,size2) NAME_EIGEN(T,NAME,size1,size2)
template <>
struct ScalarUtil<__float128> {
  typedef __float128 T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max() {
    return std::numeric_limits<double>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1E-20f;
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_inf() {
    return strtoflt128("INF",NULL);
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_nanq() {
    return strtoflt128("NAN",NULL);
  }
};
#include <CommonFile/EndAllEigen.h>

PRJ_END

#endif
