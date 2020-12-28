//meta-template function to check if a polynomial/term is zero
template <typename T>
struct IsZero {
  static bool isZero(const T& val,typename ScalarOfT<T>::Type eps) {
    return val<=eps && -val<=eps;
  }
  static T removeZero(const T& val,typename ScalarOfT<T>::Type eps) {
    return val;
  }
};
template <typename T,char LABEL>
struct IsZero<SOSPolynomial<T,LABEL>> {
  static bool isZero(const SOSPolynomial<T,LABEL>& val,typename ScalarOfT<T>::Type eps) {
    for(sizeType i=0;i<(sizeType)val._terms.size();i++)
      if(!IsZero<T>::isZero(val._terms[i]._coef,eps))
        return false;
    return true;
  }
  static SOSPolynomial<T,LABEL> removeZero(const SOSPolynomial<T,LABEL>& val,typename ScalarOfT<T>::Type eps) {
    SOSPolynomial<T,LABEL> ret;
    for(sizeType i=0;i<(sizeType)val._terms.size();i++) {
      SOSTerm<T,LABEL> t=val._terms[i];
      t._coef=IsZero<T>::removeZero(t._coef,eps);
      if(!IsZero<T>::isZero(t._coef,eps))
        ret+=t;
    }
    return ret;
  }
};
