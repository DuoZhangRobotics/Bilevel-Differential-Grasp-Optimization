//meta-template function to deduce contracted poly type
template <typename T,char LABEL2>
struct Contract {
  typedef T Type;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  static Type eval(const T& val,const COLD& x) {
    return val;
  }
};
template <typename T,char LABEL,char LABEL2>
struct Contract<SOSPolynomial<T,LABEL>,LABEL2> {
  typedef typename Contract<T,LABEL2>::Type TypeCT;
  typedef SOSPolynomial<TypeCT,LABEL> Type;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  static Type eval(const SOSPolynomial<T,LABEL>& val,const COLD& x) {
    Type ret;
    for(sizeType i=0; i<(sizeType)val._terms.size(); i++) {
      SOSTerm<TypeCT,LABEL> t(Contract<T,LABEL2>::eval(val._terms[i]._coef,x));
      t._id=val._terms[i]._id;
      t._order=val._terms[i]._order;
      ret._terms.push_back(t);
    }
    return ret;
  }
};
template <typename T,char LABEL>
struct Contract<SOSPolynomial<T,LABEL>,LABEL> {
  typedef typename Contract<T,LABEL>::Type Type;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  static Type eval(const SOSPolynomial<T,LABEL>& val,const COLD& x) {
    return val.eval(x);
  }
};
