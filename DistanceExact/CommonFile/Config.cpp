#include "Config.h"
#include <cstdio>
#include <string>
#include <cstdarg>
#include <stdexcept>
#include <iostream>

#if defined(QUADMATH_SUPPORT) && defined(__GNUC__)
namespace std
{
#define STDQUAD1(NAME) scalarD NAME(scalarD a) {return NAME##q(a);}
#define STDQUAD2(NAME) scalarD NAME(scalarD a,scalarD b) {return NAME##q(a,b);}
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
scalarD abs(scalarD a)
{
  return fabs(a);
}
bool isfinite(scalarD a)
{
  return !isinfq(a) && !isnanq(a);
}
scalarD frexp(scalarD a,int* exp)
{
  return frexpq(a,exp);
}
scalarD ldexp(scalarD a,int exp)
{
  return ldexpq(a,exp);
}
istream& operator>>(istream& input,scalarD& x)
{
  double tmp;
  input >> tmp;
  x=tmp;
  return input;
}
ostream& operator<<(ostream& output,scalarD x)
{
  output << (double)x;
  return output;
}
}
void sincos(scalarD a,scalarD* s,scalarD* c)
{
  sincosq(a,s,c);
}

//coutPrintf
bool isInt(char c)
{
  return c == 'd';
}
bool isUInt(char c)
{
  return c == 'u';
}
bool isFloat(char c)
{
  return c == 'f' || c == 'F';
}
bool isString(char c)
{
  return c == 's';
}
bool isType(char c)
{
  return c == 'd' || c == 'i' || c == 'u' || c == 'o' ||
         c == 'x' || c == 'X' || c == 'f' || c == 'F' ||
         c == 'e' || c == 'E' || c == 'g' || c == 'G' ||
         c == 'a' || c == 'A' || c == 'c' || c == 's' ||
         c == 'p' || c == 'n';
}
void coutPrintf(std::string format,...)
{
  va_list args;
  va_start(args,format);
  while(format[0] != '[')
    format.erase(format.begin());
  while(format.back() != '\n')
    format.pop_back();
  for(size_t i=0; i<format.size();)
    if(format[i] == '%') {
      size_t j=i+1;
      while(j<format.size())
        if(isType(format[j])) {
          if(isFloat(format[j])) {
            std::cout << va_arg(args,scalarD);
            break;
          } else if(isString(format[j])) {
            std::cout << va_arg(args,const char*);
            break;
          } else if(isInt(format[j])) {
            std::cout << va_arg(args,sizeType);
            break;
          } else if(isUInt(format[j])) {
            std::cout << va_arg(args,unsigned int);
            break;
          } else {
            assert(false);
          }
        } else j++;
      assert(j<format.size());
      i=j+1;
    } else {
      std::cout << format[i];
      i++;
    }
  std::cout << std::flush;
  va_end(args);
}
#else
void coutPrintf(std::string format,...)
{
  va_list args,args_copy;
  va_start(args,format);
  va_copy(args_copy,args);
  const auto sz=std::vsnprintf(nullptr,0,format.c_str(),args)+1;
  try {
    std::string result(sz,' ');
    std::vsnprintf(&result.front(),sz,format.c_str(),args_copy);
    va_end(args_copy);
    va_end(args);
    // do whatever else with result
    while(!result.empty() && result.back() != '\n')
      result.pop_back();
    std::cout << result << std::flush;
  } catch(const std::bad_alloc&) {
    va_end(args_copy);
    va_end(args);
    throw;
  }
}
#endif
