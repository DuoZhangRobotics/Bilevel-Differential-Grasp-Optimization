#ifndef MPQZ_IO_H
#define MPQZ_IO_H

#include "CommonFile/IOBasic.h"
#include "gmpxx.h"
#include "gmp.h"

PRJ_BEGIN

template <typename VEC,typename VECIN>
VEC castRational(const Eigen::MatrixBase<VECIN>& val)
{
  VEC ret;
  ret.resize(val.rows(),val.cols());
  for(sizeType r=0;r<val.rows();r++)
    for(sizeType c=0;c<val.cols();c++)
      ret(r,c)=val(r,c).get_d();
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
