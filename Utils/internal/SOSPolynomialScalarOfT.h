//meta-template function to convert scalar to poly-of-poly
template <typename T>
struct ScalarOfT {
  typedef T Type;
  template <typename T2>
  static T convert(T2 val) {
    return T(val);
  }
  static T convertVar(sizeType id,char LABEL2) {
    ASSERT(false)
    return T(0);
  }
  template <typename T2>
  static T2 convertBack(T val) {
    return T2(val);
  }
};
template <typename T,char LABEL>
struct ScalarOfT<SOSTerm<T,LABEL>> {
  typedef typename ScalarOfT<T>::Type Type;
  template <typename T2>
  static SOSTerm<T,LABEL> convert(T2 val) {
    return SOSTerm<T,LABEL>(ScalarOfT<T>::convert(val));
  }
  static SOSTerm<T,LABEL> convertVar(sizeType id,char LABEL2) {
    SOSTerm<T,LABEL> ret;
    if(LABEL==LABEL2) {
      ret._id.push_back(id);
      ret._order.push_back(1);
      ret._coef=ScalarOfT<T>::convert(1);
    } else ret._coef=ScalarOfT<T>::convertVar(id,LABEL2);
    return ret;
  }
  template <typename T2>
  static T2 convertBack(const SOSTerm<T,LABEL>& val) {
    ASSERT(val._id.empty())
    return ScalarOfT<T>::template convertBack<T2>(val._coef);
  }
};
template <typename T,char LABEL>
struct ScalarOfT<SOSPolynomial<T,LABEL>> {
  typedef typename ScalarOfT<T>::Type Type;
  template <typename T2>
  static SOSPolynomial<T,LABEL> convert(T2 val) {
    return SOSPolynomial<T,LABEL>(ScalarOfT<SOSTerm<T,LABEL>>::convert(val));
  }
  static SOSPolynomial<T,LABEL> convertVar(sizeType id,char LABEL2) {
    return SOSPolynomial<T,LABEL>(ScalarOfT<SOSTerm<T,LABEL>>::convertVar(id,LABEL2));
  }
  template<typename T2>
  static T2 convertBack(const SOSPolynomial<T,LABEL>& val) {
    ASSERT(val._terms.size()<=1)
    if(val._terms.empty())
      return T2(0);
    else return ScalarOfT<SOSTerm<T,LABEL>>::template convertBack<T2>(val._terms[0]);
  }
};
