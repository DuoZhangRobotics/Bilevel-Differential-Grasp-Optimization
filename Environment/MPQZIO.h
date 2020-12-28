#ifndef MPQZ_IO_H
#define MPQZ_IO_H

#include <Utils/Scalar.h>
#include <CommonFile/IOFwd.h>
#include <gmpxx.h>

PRJ_BEGIN

extern void castRational(const mpq_class& q,float& to);
extern void castRational(const mpq_class& q,double& to);
extern void castRational(const mpq_class& q,__float128& to);
extern void castRational(const mpq_class& q,mpfr::mpreal& to);
extern void castRational(const float& q,mpq_class& to);
extern void castRational(const double& q,mpq_class& to);
extern void castRational(const __float128& q,mpq_class& to);
extern void castRational(const mpfr::mpreal& q,mpq_class& to);
template <typename VEC,typename VECIN>
VEC castRational(const Eigen::MatrixBase<VECIN>& val)
{
  VEC ret;
  ret.resize(val.rows(),val.cols());
  for(sizeType r=0; r<val.rows(); r++)
    for(sizeType c=0; c<val.cols(); c++)
      castRational(val(r,c),ret(r,c));
  return ret;
}
std::istream& readBinaryData(mpz_class& v,std::istream& is,IOData* dat=NULL);
std::ostream& writeBinaryData(const mpz_class& v,std::ostream& os,IOData* dat=NULL);
std::istream& readBinaryData(mpq_class& v,std::istream& is,IOData* dat=NULL);
std::ostream& writeBinaryData(const mpq_class& v,std::ostream& os,IOData* dat=NULL);

#define NO_CONFLICT
#define FIXED_ONLY
#include "BeginAllEigen.h"
//redefine atomic operation
#undef NAME_EIGEN
#define NAME_EIGEN(type,NAME,size1,size2) \
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat=NULL);
//realize
NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE()
NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE()
#include "EndAllEigen.h"

PRJ_END

#endif
