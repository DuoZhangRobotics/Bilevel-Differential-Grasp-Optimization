#ifndef DEBUG_GRADIENT_H
#define DEBUG_GRADIENT_H

#include <iostream>
#include <Utils/Scalar.h>

//numeric delta
#define DEFINE_NUMERIC_DELTA_T(T) T DELTA=DELTAPrecision<T>::delta();
template <typename T>
struct DELTAPrecision;
template <>
struct DELTAPrecision<float>
{
  static float delta()
  {
    return 1e-5f;
  }
};
template <>
struct DELTAPrecision<double>
{
  static double delta()
  {
    return 1e-9f;
  }
};
template <>
struct DELTAPrecision<scalarQ>
{
  static scalarQ delta()
  {
    return 1e-15f;
  }
};
template <>
struct DELTAPrecision<mpfr::mpreal>
{
  static mpfr::mpreal delta()
  {
    mpfr_prec_t p=mpfr_get_default_prec();
    return mpfr::mpreal((std::string("1e-")+std::to_string(long(p/5))).c_str());
  }
};
#define DEFINE_NUMERIC_DELTA DEFINE_NUMERIC_DELTA_T(scalar)

//gradient debug
#define DEBUG_GRADIENT(NAME,A,B) \
if(std::abs(B) > std::sqrt(DELTA)) { \
  std::ostringstream oss;   \
  oss << "\033[31m" << NAME << ": " << A << " Err: " << B << "\033[30m"; \
  INFO(oss.str().c_str())   \
} else {  \
  std::ostringstream oss;   \
  oss << NAME << ": " << A << " Err: " << B;  \
  INFO(oss.str().c_str())   \
}

#define DEBUG_GRADIENT_REL(NAME,A,B) \
if(std::abs(B) > std::sqrt(DELTA)*std::abs(A)) { \
  std::ostringstream oss;   \
  oss << "\033[31m" << NAME << ": " << A << " Err: " << B << "\033[30m"; \
  INFO(oss.str().c_str())   \
} else {  \
  std::ostringstream oss;   \
  oss << NAME << ": " << A << " Err: " << B;  \
  INFO(oss.str().c_str())   \
}

#define DEBUG_GRADIENT_ASSERT(NAME,A,B) \
if(std::abs(B) > std::sqrt(DELTA)) { \
  std::ostringstream oss;   \
  oss << "\033[31m" << NAME << ": " << A << " Err: " << B << "\033[30m"; \
  INFO(oss.str().c_str())   \
  ASSERT(false) \
} else {  \
  std::ostringstream oss;   \
  oss << NAME << ": " << A << " Err: " << B;  \
  INFO(oss.str().c_str())   \
}

#define DEBUG_GRADIENT_REL_ASSERT(NAME,A,B) \
if(std::abs(B) > std::sqrt(DELTA)*std::abs(A)) { \
  std::ostringstream oss;   \
  oss << "\033[31m" << NAME << ": " << A << " Err: " << B << "\033[30m"; \
  INFO(oss.str().c_str())   \
  ASSERT(false) \
} else {  \
  std::ostringstream oss;   \
  oss << NAME << ": " << A << " Err: " << B;  \
  INFO(oss.str().c_str())   \
}

#endif
