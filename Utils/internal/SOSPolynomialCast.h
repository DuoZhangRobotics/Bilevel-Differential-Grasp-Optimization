//meta-template function to cast
template <typename T,typename T2>
struct Cast
{
  typedef T2 Type;
  static Type cast(const T& val) {
    return val;
  }
};
template <typename T,char LABEL,typename T2>
struct Cast<SOSTerm<T,LABEL>,T2>
{
  typedef SOSTerm<typename Cast<T,T2>::Type,LABEL> Type;
  static Type cast(const SOSTerm<T,LABEL>& val) {
    return val.template cast<T2>();
  }
};
template <typename T,char LABEL,typename T2>
struct Cast<SOSPolynomial<T,LABEL>,T2>
{
  typedef SOSPolynomial<typename Cast<T,T2>::Type,LABEL> Type;
  static Type cast(const SOSPolynomial<T,LABEL>& val) {
    return val.template cast<T2>();
  }
};
