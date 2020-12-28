#include "MPQZIO.h"
#include <CommonFile/IOBasic.h>

PRJ_BEGIN

void castRational(const mpq_class& q,float& to)
{
  to=q.get_d();
}
void castRational(const mpq_class& q,double& to)
{
  to=q.get_d();
}
void castRational(const mpq_class& q,__float128& to)
{
  to=mpfr_get_float128(mpfr::mpreal(q.get_mpq_t()).mpfr_ptr(),GMP_RNDN);
}
void castRational(const mpq_class& q,mpfr::mpreal& to)
{
  to=mpfr::mpreal(q.get_mpq_t());
}
void castRational(const float& q,mpq_class& to)
{
  to=mpq_class(q);
}
void castRational(const double& q,mpq_class& to)
{
  to=mpq_class(q);
}
void castRational(const __float128& q,mpq_class& to)
{
  mpfr::mpreal qq;
  mpfr_set_float128(qq.mpfr_ptr(),q,MPFR_RNDN);
  mpfr_get_q(to.get_mpq_t(),qq.mpfr_srcptr());
}
void castRational(const mpfr::mpreal& q,mpq_class& to)
{
  mpfr_get_q(to.get_mpq_t(),q.mpfr_srcptr());
}
std::istream& readBinaryData(mpz_class& v,std::istream& is,IOData*)
{
  std::string str;
  readBinaryData(str,is);
  v.set_str(str,10);
  return is;
}
std::ostream& writeBinaryData(const mpz_class& v,std::ostream& os,IOData*)
{
  writeBinaryData(v.get_str(),os);
  return os;
}
std::istream& readBinaryData(mpq_class& v,std::istream& is,IOData*)
{
  readBinaryData(v.get_num(),is);
  readBinaryData(v.get_den(),is);
  return is;
}
std::ostream& writeBinaryData(const mpq_class& v,std::ostream& os,IOData*)
{
  writeBinaryData(v.get_num(),os);
  writeBinaryData(v.get_den(),os);
  return os;
}

//io for fixed matrix
#define NO_CONFLICT
#define FIXED_ONLY
#include "BeginAllEigen.h"
//redefine atomic operation
#undef NAME_EIGEN
#define NAME_EIGEN(type,NAME,size1,size2) \
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData*) \
{ \
  typedef Eigen::Matrix<type,size1,size2> TYPE; \
  sizeType d0=TYPE::RowsAtCompileTime;os.write((char*)&d0,sizeof(sizeType)); \
  sizeType d1=TYPE::ColsAtCompileTime;os.write((char*)&d1,sizeof(sizeType)); \
  for(sizeType r=0;r<TYPE::RowsAtCompileTime;r++) \
  for(sizeType c=0;c<TYPE::ColsAtCompileTime;c++) \
    writeBinaryData(v(r,c),os); \
  return os; \
} \
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData*) \
{ \
  typedef Eigen::Matrix<type,size1,size2> TYPE; \
  sizeType d0;is.read((char*)&d0,sizeof(sizeType));ASSERT(d0 == TYPE::RowsAtCompileTime) \
  sizeType d1;is.read((char*)&d1,sizeof(sizeType));ASSERT(d1 == TYPE::ColsAtCompileTime) \
  for(sizeType r=0;r<TYPE::RowsAtCompileTime;r++) \
  for(sizeType c=0;c<TYPE::ColsAtCompileTime;c++) \
    readBinaryData(v(r,c),is);  \
  return is; \
}
//realize
NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE()
NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE()
#include "EndAllEigen.h"

PRJ_END
