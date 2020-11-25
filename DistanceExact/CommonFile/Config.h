#ifndef CONFIG_H
#define CONFIG_H

#define NAMESPACE COMMON
#define PRJ_BEGIN namespace NAMESPACE {
#define PRJ_END	}
#define USE_PRJ_NAMESPACE using namespace NAMESPACE;

//#define _DEBUG
#ifdef _DEBUG
#define CUSTOM_DEBUG
#endif
#define PROFILE

#ifdef _MSC_VER
#define ALIGN_8 __declspec(align(8))
#define ALIGN_16 __declspec(align(16))
#define FORCE_INLINE __forceinline
#else
#define ALIGN_8 __attribute__((aligned (8)))
#define ALIGN_16 __attribute__((aligned (16)))
#define FORCE_INLINE inline
#endif

#include <string>
void coutPrintf(const std::string fmt,...);
//#define NDEBUG
#ifndef NDEBUG

#include <assert.h>
#define ASSERT(x) do{assert((x));}while(0);
#define ASSERT_MSG(x,msg) do{if(!(x)){coutPrintf("[ERROR] %s \n",msg);fflush(stdout);assert(false);}}while(0);
#define ASSERT_MSGV(x,fmt,...) do{if(!(x)){coutPrintf("[ERROR] " fmt " \n",__VA_ARGS__);fflush(stdout);assert(false);}}while(0);

#else

#ifdef _MSC_VER
#pragma warning(disable:4552)
#pragma warning(disable:4553)
#endif
#define ASSERT(x) do{if(!(x)){exit(EXIT_FAILURE);}}while(0);
#define ASSERT_MSG(x,msg) do{if(!(x)){WARNING(msg);exit(EXIT_FAILURE);}}while(0);
#define ASSERT_MSGV(x,fmt,...) do{if(!(x)){WARNINGV(fmt,__VA_ARGS__);exit(EXIT_FAILURE);}}while(0);

#endif

#define WARNING(msg) do{coutPrintf("[WARNING] \x1B[31m %s \x1B[0m\n",msg);fflush(stdout);}while(0);
#define WARNINGV(fmt,...) do{coutPrintf("[WARNING] \x1B[31m " fmt " \x1B[0m\n",__VA_ARGS__);fflush(stdout);}while(0);
#define INFO(msg) do{coutPrintf("[INFO] %s \n",msg);fflush(stdout);}while(0);
#define INFOV(fmt,...) do{coutPrintf("[INFO] " fmt " \n",__VA_ARGS__);fflush(stdout);}while(0);
#define NOTIFY_MSG(msg) do{coutPrintf("[NOTIFY] %s \n",msg);fflush(stdout);}while(0); getchar();
#define NOTIFY_MSGV(fmt,...) do{coutPrintf("[NOTIFY] " fmt " \n",__VA_ARGS__);fflush(stdout);}while(0); getchar();

//OpenMP only support signed variable as index
#include <stdint.h>
#define USE_QUAD_SIZE
#ifdef USE_QUAD_SIZE
typedef int64_t sizeType;
#else
typedef int sizeType;
#endif

#if defined(QUADMATH_SUPPORT) && defined(__GNUC__)
#include <limits>
#include <iostream>
#include <quadmath.h>
#define DOUBLE_PRECISION
typedef __float128 scalarD;
typedef float scalarF;
namespace std
{
#define STDQUAD1(NAME) scalarD NAME(scalarD a);
#define STDQUAD2(NAME) scalarD NAME(scalarD a,scalarD b);
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
scalarD abs(scalarD a);
bool isfinite(scalarD a);
scalarD frexp(scalarD a,int* exp);
scalarD ldexp(scalarD a,int exp);
istream& operator>>(istream& input,scalarD& x);
ostream& operator<<(ostream& output,scalarD x);
}
void sincos(scalarD a,scalarD* s,scalarD* c);
#include <Eigen/Core>
namespace Eigen
{
template<> struct GenericNumTraits<__float128>
{
  enum {
    IsInteger = false,
    IsSigned = true,
    IsComplex = 0,
    RequireInitialization = 0,
    ReadCost = 1,
    AddCost = 1,
    MulCost = 1
  };

  typedef __float128 Real;
  typedef __float128 Literal;
  typedef __float128 Nested;

  static inline Real epsilon()
  {
    return std::numeric_limits<Real>::epsilon();
  }
  static inline Real dummy_precision()
  {
    return Real(0);
  }
  static inline __float128 highest()
  {
    return std::numeric_limits<Real>::max();
  }
  static inline __float128 lowest()
  {
    return std::numeric_limits<Real>::lowest();
  }
  static inline int digits10()
  {
    return std::numeric_limits<Real>::digits10;
  }
};
}
#else
typedef double scalarD;
typedef float scalarF;
#endif

typedef int vtkSizeType;
#define INVALID sizeType(-1)

#define FUNCTION_NOT_IMPLEMENTED ASSERT_MSGV(false,"Function \"%s\" not implemented!",__FUNCTION__)

#endif
