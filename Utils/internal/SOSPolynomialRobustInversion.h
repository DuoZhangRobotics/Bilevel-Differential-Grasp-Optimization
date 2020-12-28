//meta-template function to deduce contracted poly type
template <typename T,char LABEL>
struct RobustInversion {
  typedef typename SOSPolynomial<T,LABEL>::MATD MATD;
  static MATD eval(MATD& m) {
    return m.inverse();
  }
};
template <char LABEL>
struct RobustInversion<__float128,LABEL> {
  typedef typename SOSPolynomial<__float128,LABEL>::MATD MATD;
  static MATD eval(MATD& m) {
    return m.template cast<double>().inverse().template cast<__float128>();
  }
};
