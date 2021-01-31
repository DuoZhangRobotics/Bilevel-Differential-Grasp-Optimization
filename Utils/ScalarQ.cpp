#include "ScalarQ.h"

namespace std
{
#define STDQUAD1(NAME) __float128 NAME(__float128 a) {return NAME##q(a);}
#define STDQUAD2(NAME) __float128 NAME(__float128 a,__float128 b) {return NAME##q(a,b);}
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
__float128 abs(const __float128& a)
{
  return fabs(a);
}
bool isfinite(__float128 a)
{
  return !isinfq(a) && !isnanq(a);
}
__float128 frexp(__float128 a,int* exp)
{
  return frexpq(a,exp);
}
__float128 ldexp(__float128 a,int exp)
{
  return ldexpq(a,exp);
}
void convert_scalar(double a,double& to)
{
  to=a;
}
void convert_scalar(double a,__float128& to)
{
  to=a;
}
void convert_scalar(__float128 a,double& to)
{
  to=a;
}
void convert_scalar(__float128 a,__float128& to)
{
  to=a;
}
double to_double(__float128 a)
{
  return a;
}
std::string to_string(__float128 a)
{
  ostringstream oss;
  oss << a;
  return oss.str();
}
istream& operator>>(istream& input,__float128& x)
{
  double tmp;
  input >> tmp;
  x=tmp;
  return input;
}
ostream& operator<<(ostream& output,__float128 x)
{
  output << (double)x;
  return output;
}
}

PRJ_BEGIN

std::ostream& writeBinaryData(const __float128& val,std::ostream& os,IOData*)
{
  double valD=std::to_double(val);
  os.write((char*)&valD,sizeof(double));
  return os;
}
std::istream& readBinaryData(__float128& val,std::istream& is,IOData*)
{
  double valD;
  is.read((char*)&valD,sizeof(double));
  val=valD;
  return is;
}

PRJ_END
