#include "ScalarMpfr.h"

namespace std
{
mpfr::mpreal min(mpfr::mpreal a,mpfr::mpreal b) {
  return a<b?a:b;
}
mpfr::mpreal max(mpfr::mpreal a,mpfr::mpreal b) {
  return a>b?a:b;
}
mpfr::mpreal pow(mpfr::mpreal a,mpfr::mpreal p) {
  return mpfr::pow(a,p);
}
mpfr::mpreal sqrt(mpfr::mpreal a) {
  return mpfr::sqrt(a);
}
sizeType floor(mpfr::mpreal a) {
  return mpfr::floor(a).toLong();
}
sizeType ceil(mpfr::mpreal a) {
  return mpfr::ceil(a).toLong();
}
sizeType round(mpfr::mpreal a) {
  return mpfr::round(a).toLong();
}
mpfr::mpreal log(mpfr::mpreal a) {
  return mpfr::log(a);
}
mpfr::mpreal exp(mpfr::mpreal a) {
  return mpfr::exp(a);
}
mpfr::mpreal abs(mpfr::mpreal a) {
  return mpfr::abs(a);
}
mpfr::mpreal sin(mpfr::mpreal a) {
  return mpfr::sin(a);
}
mpfr::mpreal cos(mpfr::mpreal a) {
  return mpfr::cos(a);
}
mpfr::mpreal asin(mpfr::mpreal a) {
  return mpfr::asin(a);
}
mpfr::mpreal acos(mpfr::mpreal a) {
  return mpfr::acos(a);
}
double to_double(mpfr::mpreal a) {
  return a.toDouble();
}
void convert_scalar(double a,mpfr::mpreal& to) {
  to=mpfr::mpreal(a);
}
void convert_scalar(__float128 a,mpfr::mpreal& to) {
  mpfr_set_float128(to.mpfr_ptr(),a,MPFR_RNDN);
}
void convert_scalar(mpfr::mpreal a,mpfr::mpreal& to) {
  to=a;
}
void convert_scalar(mpfr::mpreal a,double& to) {
  to=std::to_double(a);
}
void convert_scalar(mpfr::mpreal a,__float128& to) {
  to=mpfr_get_float128(a.mpfr_ptr(),MPFR_RNDN);
}
std::string to_string(mpfr::mpreal a) {
  ostringstream oss;
  oss << a;
  return oss.str();
}
bool isfinite(mpfr::mpreal a) {
  return mpfr::isfinite(a);
}
}

PRJ_BEGIN

std::ostream& writeBinaryData(const mpfr::mpreal& val,std::ostream& os,IOData*)
{
  double valD=std::to_double(val);
  os.write((char*)&valD,sizeof(double));
  return os;
}
std::istream& readBinaryData(mpfr::mpreal& val,std::istream& is,IOData*)
{
  double valD;
  is.read((char*)&valD,sizeof(double));
  val=valD;
  return is;
}

PRJ_END
