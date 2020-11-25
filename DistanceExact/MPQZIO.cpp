#include "MPQZIO.h"

PRJ_BEGIN

std::istream& readBinaryData(mpz_class& v,std::istream& is,IOData* dat)
{
  std::string str;
  readBinaryData(str,is);
  v.set_str(str,10);
  return is;
}
std::ostream& writeBinaryData(const mpz_class& v,std::ostream& os,IOData* dat)
{
  writeBinaryData(v.get_str(),os);
  return os;
}
std::istream& readBinaryData(mpq_class& v,std::istream& is,IOData* dat)
{
  readBinaryData(v.get_num(),is);
  readBinaryData(v.get_den(),is);
  return is;
}
std::ostream& writeBinaryData(const mpq_class& v,std::ostream& os,IOData* dat)
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
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat) \
{ \
  typedef Eigen::Matrix<type,size1,size2> TYPE; \
  sizeType d0=TYPE::RowsAtCompileTime;os.write((char*)&d0,sizeof(sizeType)); \
  sizeType d1=TYPE::ColsAtCompileTime;os.write((char*)&d1,sizeof(sizeType)); \
  for(sizeType r=0;r<TYPE::RowsAtCompileTime;r++) \
  for(sizeType c=0;c<TYPE::ColsAtCompileTime;c++) \
    writeBinaryData(v(r,c),os); \
  return os; \
} \
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat) \
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
