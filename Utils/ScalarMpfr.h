#ifndef SCALAR_MPFR_H
#define SCALAR_MPFR_H

#define MPFR_WANT_FLOAT128
#include <mpreal.h>
#include <quadmath.h>
#include <unsupported/Eigen/MPRealSupport>
#include <CommonFile/MathBasic.h>
#include <CommonFile/IOFwd.h>
namespace std
{
mpfr::mpreal min(mpfr::mpreal a,mpfr::mpreal b);
mpfr::mpreal max(mpfr::mpreal a,mpfr::mpreal b);
mpfr::mpreal pow(mpfr::mpreal a,mpfr::mpreal p);
sizeType floor(mpfr::mpreal a);
sizeType ceil(mpfr::mpreal a);
mpfr::mpreal sqrt(mpfr::mpreal a);
sizeType round(mpfr::mpreal a);
mpfr::mpreal log(mpfr::mpreal a);
mpfr::mpreal exp(mpfr::mpreal a);
mpfr::mpreal abs(mpfr::mpreal a);
mpfr::mpreal sin(mpfr::mpreal a);
mpfr::mpreal cos(mpfr::mpreal a);
mpfr::mpreal asin(mpfr::mpreal a);
mpfr::mpreal acos(mpfr::mpreal a);
double to_double(mpfr::mpreal a);
void convert_scalar(double a,mpfr::mpreal& to);
void convert_scalar(__float128 a,mpfr::mpreal& to);
void convert_scalar(mpfr::mpreal a,mpfr::mpreal& to);
void convert_scalar(mpfr::mpreal a,double& to);
void convert_scalar(mpfr::mpreal a,__float128& to);
std::string to_string(mpfr::mpreal a);
bool isfinite(mpfr::mpreal a);
}

PRJ_BEGIN

#include <CommonFile/BeginAllEigen.h>
#undef NAME_EIGEN_ALLTYPES
#define NAME_EIGEN_ALLTYPES(NAME,size1,size2) NAME_EIGEN(T,NAME,size1,size2)
template <>
struct ScalarUtil<mpfr::mpreal> {
  typedef mpfr::mpreal T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE static T scalar_max() {
    return std::numeric_limits<double>::max();
  }
  FORCE_INLINE static T scalar_eps() {
    return 1E-20f;
  }
  FORCE_INLINE static T scalar_inf() {
    T ret;
    mpfr_set_float128(ret.mpfr_ptr(),strtoflt128("INF",NULL),MPFR_RNDN);
    return ret;
  }
  FORCE_INLINE static T scalar_nanq() {
    T ret;
    mpfr_set_float128(ret.mpfr_ptr(),strtoflt128("NAN",NULL),MPFR_RNDN);
    return ret;
  }
};
#include <CommonFile/EndAllEigen.h>

std::ostream& writeBinaryData(mpfr::mpreal val,std::ostream& os,IOData*);
std::istream& readBinaryData(mpfr::mpreal& val,std::istream& is,IOData*);

PRJ_END

#endif
