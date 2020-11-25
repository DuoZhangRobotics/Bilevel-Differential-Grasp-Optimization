#ifndef MATH_BASIC_H
#define MATH_BASIC_H

#include "Config.h"
#include <vector>
#ifndef __APPLE__
#include <omp.h>
#else
#define NO_OPENMP
#endif
#include <cmath>
#include <float.h>
#include <Eigen/Core>
#include <Eigen/Geometry>
#define EIGEN_DEVICE_FUNC
#define DEVICE_ONLY_FUNC

//convert function for IO.h
namespace std
{
template <typename T>
struct convert {
  template <typename T2>
  FORCE_INLINE T operator()(T2 d) const {
    return (T)d;
  }
  template <typename T2,int r,int c,int o,int mr,int mc>
  FORCE_INLINE T operator()(Eigen::Matrix<T2,r,c,o,mr,mc> d) const {
    ASSERT(false) //this is just a trait for compilation, it do nothing
    return T();
  }
};
}

#ifdef M_PI
#undef M_PI
#endif

PRJ_BEGIN

#ifdef QUADMATH_SUPPORT
typedef scalarD scalar;
#define M_PI scalar(3.14159265358979323846)
#elif defined(DOUBLE_PRECISION)
typedef scalarD scalar;
#define M_PI scalar(3.14159265358979323846)
#else
typedef scalarF scalar;
#define M_PI scalar(3.14159265358979323846f)
#endif

//define all eigen special sized matrices
#include "BeginAllEigen.h"
NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE()
NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE()
NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(,f,scalarF)
NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(,d,scalarD)
#undef NAME_EIGEN_ALLTYPES
#define NAME_EIGEN_ALLTYPES(NAME,size1,size2) NAME_EIGEN(T,NAME,size1,size2)
template <typename T>
struct ScalarUtil;
template <>
struct ScalarUtil<float> {
  typedef float T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max()
  {
    return std::numeric_limits<float>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1E-5f;
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_inf()
  {
    return std::numeric_limits<float>::infinity();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_nanq()
  {
    return std::numeric_limits<float>::quiet_NaN();
  }
};
template <>
struct ScalarUtil<double> {
  typedef double T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max()
  {
    return std::numeric_limits<double>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1E-10f;
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_inf()
  {
    return std::numeric_limits<double>::infinity();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_nanq()
  {
    return std::numeric_limits<double>::quiet_NaN();
  }
};
#if defined(QUADMATH_SUPPORT) && defined(__GNUC__)
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
#endif
template <>
struct ScalarUtil<char> {
  typedef char T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max()
  {
    return std::numeric_limits<char>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1;
  }
};
template <>
struct ScalarUtil<unsigned char> {
  typedef unsigned char T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max()
  {
    return std::numeric_limits<unsigned char>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1;
  }
};
template <>
struct ScalarUtil<int> {
  typedef int T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max()
  {
    return std::numeric_limits<int>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1;
  }
};
template <>
struct ScalarUtil<unsigned int> {
  typedef unsigned int T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max()
  {
    return std::numeric_limits<unsigned int>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1;
  }
};
template <>
struct ScalarUtil<int64_t> {
  typedef int64_t T;
  NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE(Scalar)
  NAME_EIGEN_SPECIAL_ALLTYPES_SPECIALSIZE(Scalar,,T)
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_max()
  {
    return std::numeric_limits<int64_t>::max();
  }
  FORCE_INLINE DEVICE_ONLY_FUNC static T scalar_eps() {
    return 1;
  }
};
template<typename T> sizeType getEigenStride()
{
  std::vector<T,Eigen::aligned_allocator<T> > tester(2);
  return (sizeType)(((char*)&(tester[1]))-((char*)&(tester[0])));
}
#include "EndAllEigen.h"

//basic elementwise operation
#define DECL_FLOAT_TRAITS(T)  \
bool compL(T a,T b);  \
bool compLE(T a,T b); \
bool compG(T a,T b);  \
bool compGE(T a,T b); \
T compMin(T a,T b); \
T compMax(T a,T b); \
T floorV(T a);  \
T ceilV(T a);
DECL_FLOAT_TRAITS(char)
DECL_FLOAT_TRAITS(unsigned char)
DECL_FLOAT_TRAITS(scalarF)
DECL_FLOAT_TRAITS(scalarD)
DECL_FLOAT_TRAITS(sizeType)
#undef DECL_FLOAT_TRAITS
template <typename Derived>
EIGEN_DEVICE_FUNC bool isFinite(const Eigen::MatrixBase<Derived>& x)
{
  //for(sizeType i=0; i<x.rows(); i++)
  //  for(sizeType j=0; j<x.cols(); j++)
  //    if(isnan(x(i,j)) || isinf(x(i,j)))
  //      return false;
  //return true;
  return ((x-x).array()==(x-x).array()).all();
}
template <typename Derived,typename Derived2>
EIGEN_DEVICE_FUNC bool compL(const Eigen::MatrixBase<Derived>& a,const Eigen::MatrixBase<Derived2>& b)
{
  return a.array().binaryExpr(b.array(),std::less<typename Derived::Scalar>()).all();
}
template <typename Derived,typename Derived2>
EIGEN_DEVICE_FUNC bool compLE(const Eigen::MatrixBase<Derived>& a,const Eigen::MatrixBase<Derived2>& b)
{
  return a.array().binaryExpr(b.array(),std::less_equal<typename Derived::Scalar>()).all();
}
template <typename Derived,typename Derived2>
EIGEN_DEVICE_FUNC bool compG(const Eigen::MatrixBase<Derived>& a,const Eigen::MatrixBase<Derived2>& b)
{
  return a.array().binaryExpr(b.array(),std::greater<typename Derived::Scalar>()).all();
}
template <typename Derived,typename Derived2>
EIGEN_DEVICE_FUNC bool compGE(const Eigen::MatrixBase<Derived>& a,const Eigen::MatrixBase<Derived2>& b)
{
  return a.array().binaryExpr(b.array(),std::greater_equal<typename Derived::Scalar>()).all();
}
template <typename Derived,typename Derived2>
EIGEN_DEVICE_FUNC Eigen::Matrix<typename Derived::Scalar,Derived::RowsAtCompileTime,Derived::ColsAtCompileTime>
compMax(const Eigen::MatrixBase<Derived>& a,const Eigen::MatrixBase<Derived2>& b)
{
  return a.cwiseMax(b);
}
template <typename Derived,typename Derived2>
EIGEN_DEVICE_FUNC Eigen::Matrix<typename Derived::Scalar,Derived::RowsAtCompileTime,Derived::ColsAtCompileTime>
compMin(const Eigen::MatrixBase<Derived>& a,const Eigen::MatrixBase<Derived2>& b)
{
  return a.cwiseMin(b);
}
template <typename Derived>
EIGEN_DEVICE_FUNC Eigen::Matrix<sizeType,Derived::RowsAtCompileTime,Derived::ColsAtCompileTime>
ceilV(const Eigen::MatrixBase<Derived>& a)
{
  return a.array().unaryExpr([](typename Derived::Scalar elem) {
    return std::ceil(elem);
  }).matrix().template cast<sizeType>();
}
template <typename Derived>
EIGEN_DEVICE_FUNC Eigen::Matrix<sizeType,Derived::RowsAtCompileTime,Derived::ColsAtCompileTime>
floorV(const Eigen::MatrixBase<Derived>& a)
{
  return a.array().unaryExpr([](typename Derived::Scalar elem) {
    return std::floor(elem);
  }).matrix().template cast<sizeType>();
}
template <typename Derived>
EIGEN_DEVICE_FUNC Derived roundV(const Eigen::MatrixBase<Derived>& a)
{
  return floorV((a.array()+0.5f).matrix());
}
//warping
template <typename T>
EIGEN_DEVICE_FUNC T handleScalarWarp(T& x,T minC,T maxC,sizeType* warpTimeRet)
{
  T warpRange=maxC-minC;
  sizeType warpTime=(sizeType)-std::floor((x-minC)/warpRange);
  if(warpTimeRet)
    *warpTimeRet=warpTime;
  x+=(T)warpTime*warpRange;
  return warpRange;
}
template <typename T>
EIGEN_DEVICE_FUNC T handleScalarWarpDiff(T& x,T warpRange,sizeType* warpTimeRet)
{
  sizeType warpTime=(sizeType)-std::floor((x+warpRange*0.5f)/warpRange);
  if(warpTimeRet)
    *warpTimeRet=warpTime;
  x+=(T)warpTime*warpRange;
  return warpRange;
}
template <typename T>
EIGEN_DEVICE_FUNC T handleScalarWarpDiff(T& x,T minC,T maxC,sizeType* warpTimeRet)
{
  return handleScalarWarpDiff<T>(x,maxC-minC,warpTimeRet);
}
//fetch sign
template <typename T>
EIGEN_DEVICE_FUNC T sgn(T val)
{
  return val < 0.0f ? -1 : val > 0.0f ? 1 : 0;
}
template <typename Derived>
EIGEN_DEVICE_FUNC Derived sgn(const Eigen::MatrixBase<Derived>& a)
{
  return a.array().unaryExpr([](typename Derived::Scalar elem) {
    return sgn(elem);
  }).matrix();
}
//align power of two
template <typename T>
EIGEN_DEVICE_FUNC T makePOT(T val,T decimate)
{
  ASSERT(decimate > 1)
  T ret=1;
  while(ret < val)
    ret*=decimate;
  return ret;
}
template <typename Derived>
EIGEN_DEVICE_FUNC Derived makePOTV(const Eigen::MatrixBase<Derived>& a,sizeType decimate)
{
  return a.array().uniaryExpr(makePOT).matrix();
}
//minmod clamp function
template <typename T>
EIGEN_DEVICE_FUNC T minmod(T a,T b,T c)
{
  if(a > 0.0f && b > 0.0f && c > 0.0f)
    return std::min<T>(std::min<T>(a,b),c);
  else if(a < 0.0f && b < 0.0f && c < 0.0f)
    return std::max<T>(std::max<T>(a,b),c);
  else return 0.0f;
}
template <typename Derived>
EIGEN_DEVICE_FUNC Eigen::Matrix<typename Derived::Scalar,Derived::RowsAtCompileTime,Derived::ColsAtCompileTime>
minmodV(const Eigen::MatrixBase<Derived>& a,const Eigen::MatrixBase<Derived>& b,const Eigen::MatrixBase<Derived>& c)
{
  Eigen::Matrix<typename Derived::Scalar,Derived::RowsAtCompileTime,Derived::ColsAtCompileTime> ret;
  ret.resize(a.rows(),a.cols());
  for(sizeType R=0; R<a.rows(); R++)
    for(sizeType C=0; C<a.cols(); C++)
      ret(R,C)=minmod(a(R,C),b(R,C),c(R,C));
  return ret;
}
//RandEngine
class RandEngine
{
public:
  //real
  static scalar randR01();
  static scalar randR(scalar l,scalar u);
  static scalar randSR(sizeType seed,scalar l,scalar u);
  //int
  static sizeType randSI(sizeType seed);
  static sizeType randI(sizeType l,sizeType u);
  static sizeType randSI(sizeType seed,sizeType l,sizeType u);
  //normal
  static scalar normal();
  static scalar normal(scalar mean,scalar cov);
  //settings
  static void seedTime();
  static void seed(sizeType i);
  static void useDeterministic();
  static void useNonDeterministic();
  //common random number
  static const std::vector<int>& getRecordedRandom(sizeType id);
  static void resetAllRecordRandom();
  static void resetRecordRandom(sizeType id);
  static void beginRecordRandom(sizeType id);
  static void endRecordRandom();
  static void clearAllRecord();
  static void clearRecord(sizeType id);
  static std::vector<sizeType> getRecordId();
  static void printRecordId();
};
//matlab funcs
template <typename T>
typename ScalarUtil<T>::ScalarVec3 randBary()
{
  typename ScalarUtil<T>::ScalarVec3 ret;

  ret.x()=(T)RandEngine::randR(0,1.0f);
  ret.y()=(T)RandEngine::randR(0,1.0f);
  ret.y()*=(1.0f-ret.x());
  ret.z()=1.0f-ret.x()-ret.y();

  return ret;
}
template <typename VECA,typename VECB>
EIGEN_DEVICE_FUNC Eigen::Matrix<typename VECA::Scalar,-1,1> concat(const VECA& A,const VECB& B)
{
  Eigen::Matrix<typename VECA::Scalar,-1,1> ret;
  ret.resize(A.size()+B.size());
  ret << A,B;
  return ret;
}
template <typename MATA,typename MATB>
EIGEN_DEVICE_FUNC Eigen::Matrix<typename MATA::Scalar,-1,-1> concatRow(const MATA& A,const MATB& B)
{
  Eigen::Matrix<typename MATA::Scalar,-1,-1> ret;
  ret.resize(A.rows()+B.rows(),A.cols());
  ret << A,B;
  return ret;
}
template <typename MATA,typename MATB>
EIGEN_DEVICE_FUNC Eigen::Matrix<typename MATA::Scalar,-1,-1> concatCol(const MATA& A,const MATB& B)
{
  Eigen::Matrix<typename MATA::Scalar,-1,-1> ret;
  ret.resize(A.rows(),A.cols()+B.cols());
  ret << A,B;
  return ret;
}
//bit treaking
template<typename RET>
EIGEN_DEVICE_FUNC RET countBits(const unsigned char& v)
{
  static const unsigned char BitsSetTable256[256] = {
#   define B2(n) n,     n+1,     n+1,     n+2
#   define B4(n) B2(n), B2(n+1), B2(n+1), B2(n+2)
#   define B6(n) B4(n), B4(n+1), B4(n+1), B4(n+2)
    B6(0), B6(1), B6(1), B6(2)
  };
  return BitsSetTable256[v];
#undef B2
#undef B4
#undef B6
}
template<typename RET>
EIGEN_DEVICE_FUNC RET countBits(const sizeType& v)
{
  RET ret=0;
  unsigned char* vc=(unsigned char*)&v;
  for(sizeType i=0; i<sizeof(sizeType); i++)
    ret+=countBits<RET>(vc[i]);
  return ret;
}
//eigen traits
template <typename T2>
struct EigenTraits {
  typedef T2 ScalarType;
  static EIGEN_DEVICE_FUNC T2 value() {
    return (T2)0.0f;
  }
  static EIGEN_DEVICE_FUNC T2 constant(T2 val) {
    return val;
  }
  static EIGEN_DEVICE_FUNC T2 max() {
    return ScalarUtil<T2>::scalar_max();
  }
  static EIGEN_DEVICE_FUNC T2 eps() {
    return ScalarUtil<T2>::scalar_eps();
  }
  static EIGEN_DEVICE_FUNC T2 add(T2 a,T2 b) {
    return a+b;
  }
  static EIGEN_DEVICE_FUNC T2 sub(T2 a,T2 b) {
    return a-b;
  }
  static EIGEN_DEVICE_FUNC T2 div(T2 a,T2 b) {
    return a/b;
  }
  static EIGEN_DEVICE_FUNC T2 mul(T2 a,T2 b) {
    return a*b;
  }
  static EIGEN_DEVICE_FUNC T2 abs(T2 a) {
    return std::abs(a);
  }
  static EIGEN_DEVICE_FUNC T2 sqrt(T2 a) {
    return std::sqrt(a);
  }
};
template <typename T2,int r,int c,int o,int mr,int mc>
struct EigenTraits<Eigen::Matrix<T2,r,c,o,mr,mc> > {
  typedef Eigen::Matrix<T2,r,c,o,mr,mc> ValueType;
  typedef T2 ScalarType;
  static EIGEN_DEVICE_FUNC ValueType value() {
    ValueType ret;
    ret.setZero();
    return ret;
  }
  static EIGEN_DEVICE_FUNC ValueType constant(T2 val) {
    return ValueType::Constant(val);
  }
  static EIGEN_DEVICE_FUNC ValueType max() {
    return ValueType::Constant(ScalarUtil<T2>::scalar_max());
  }
  static EIGEN_DEVICE_FUNC ValueType eps() {
    return ValueType::Constant(ScalarUtil<T2>::scalar_eps());
  }
  static EIGEN_DEVICE_FUNC ValueType add(ValueType a,ValueType b) {
    return (a.array()+b.array()).matrix();
  }
  static EIGEN_DEVICE_FUNC ValueType sub(ValueType a,ValueType b) {
    return (a.array()-b.array()).matrix();
  }
  template <typename T3>
  static EIGEN_DEVICE_FUNC ValueType add(ValueType a,T3 b) {
    return (a.array()+(ScalarType)b).matrix();
  }
  template <typename T3>
  static EIGEN_DEVICE_FUNC ValueType sub(ValueType a,T3 b) {
    return (a.array()-(ScalarType)b).matrix();
  }
  static EIGEN_DEVICE_FUNC ValueType div(ValueType a,ValueType b) {
    return (a.array()/b.array()).matrix();
  }
  static EIGEN_DEVICE_FUNC ValueType mul(ValueType a,ValueType b) {
    return (a.array()*b.array()).matrix();
  }
  static EIGEN_DEVICE_FUNC ValueType abs(ValueType a) {
    return a.cwiseAbs();
  }
  static EIGEN_DEVICE_FUNC ValueType sqrt(ValueType a) {
    return a.array().sqrt().matrix();
  }
};

//openmp macro prepare
#ifdef _MSC_VER
#define STRINGIFY_OMP(X) #X
#define PRAGMA __pragma
#else
#define STRINGIFY_OMP(X) #X
#define PRAGMA _Pragma
#endif
#ifndef NO_OPENMP
//openmp convenient functions
#define OMP_PARALLEL_FOR_ PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(OmpSettings::getOmpSettings().nrThreads()) schedule(static)))
#define OMP_PARALLEL_FOR_I(...) PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(OmpSettings::getOmpSettings().nrThreads()) schedule(static) __VA_ARGS__))
#define OMP_PARALLEL_FOR_X(X) PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(X) schedule(static)))
#define OMP_PARALLEL_FOR_XI(X,...) PRAGMA(STRINGIFY_OMP(omp parallel for num_threads(X) schedule(static) __VA_ARGS__))
#define OMP_ADD(...) reduction(+: __VA_ARGS__)
#define OMP_PRI(...) private(__VA_ARGS__)
#define OMP_FPRI(...) firstprivate(__VA_ARGS__)
#define OMP_ATOMIC_	PRAGMA(STRINGIFY_OMP(omp atomic))
#define OMP_ATOMIC_CAPTURE_	PRAGMA(STRINGIFY_OMP(omp atomic capture))
#define OMP_CRITICAL_ PRAGMA(STRINGIFY_OMP(omp critical))
#define OMP_FLUSH_(X) PRAGMA(STRINGIFY_OMP(omp flush(X)))
#else
//openmp convenient functions
#define OMP_PARALLEL_FOR_
#define OMP_PARALLEL_FOR_I(...)
#define OMP_PARALLEL_FOR_X(X)
#define OMP_PARALLEL_FOR_XI(X,...)
#define OMP_ADD(...)
#define OMP_PRI(...)
#define OMP_FPRI(...)
#define OMP_ATOMIC_
#define OMP_ATOMIC_CAPTURE_
#define OMP_CRITICAL_
#define OMP_FLUSH_(X)
#endif
struct OmpSettings {
public:
  static const OmpSettings& getOmpSettings();
  static OmpSettings& getOmpSettingsNonConst();
  int nrThreads() const;
  int threadId() const;
  void setNrThreads(int nr);
  void useAllThreads();
private:
  OmpSettings();
  int _nrThreads;
  static OmpSettings _ompSettings;
};

PRJ_END

#endif
