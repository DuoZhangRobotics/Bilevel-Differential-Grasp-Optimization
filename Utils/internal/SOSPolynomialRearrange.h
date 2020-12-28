//meta-template function to rearrange SOSPolynomials
template <typename T,char LABEL,typename T2>
struct RearrangeAdd {
  static SOSTerm<T,LABEL> add(T2 coef,const std::unordered_map<char,SOSInfo>& m) {
    SOSTerm<T,LABEL> ret;
    typename std::unordered_map<char,SOSInfo>::const_iterator it=m.find(LABEL);
    if(it!=m.end()) {
      ret._id=it->second._id;
      ret._order=it->second._order;
    }
    ret._coef=coef;
    //ASSERT(it!=m.end())   //so that we can deal with label deficient case
    return ret;
  }
};
template <typename T3,char LABEL3,char LABEL,typename T2>
struct RearrangeAdd<SOSPolynomial<T3,LABEL3>,LABEL,T2> {
  static SOSTerm<SOSPolynomial<T3,LABEL3>,LABEL> add(T2 coef,const std::unordered_map<char,SOSInfo>& m) {
    SOSTerm<SOSPolynomial<T3,LABEL3>,LABEL> ret;
    std::unordered_map<char,SOSInfo>::const_iterator it=m.find(LABEL);
    if(it!=m.end()) {
      ret._id=it->second._id;
      ret._order=it->second._order;
    }
    ret._coef=RearrangeAdd<T3,LABEL3,T2>::add(coef,m);
    //ASSERT(it!=m.end())   //so that we can deal with label deficient case
    return ret;
  }
};
template <typename T,char LABEL,typename T2,char LABEL2>
struct Rearrange {
  typedef SOSPolynomial<T,LABEL> OUT;
  typedef SOSPolynomial<T2,LABEL2> IN;
  static void arrange(OUT& to,const IN& from,std::unordered_map<char,SOSInfo> m) {
    for(sizeType t=0; t<(sizeType)from._terms.size(); t++) {
      std::unordered_map<char,SOSInfo> mt=m;
      mt[LABEL2]=from._terms[t];
      to.add(RearrangeAdd<T,LABEL,T2>::add(from._terms[t]._coef,mt));
    }
  }
};
template <typename T,char LABEL,typename T3,char LABEL3,char LABEL2>
struct Rearrange<T,LABEL,SOSPolynomial<T3,LABEL3>,LABEL2> {
  typedef SOSPolynomial<T,LABEL> OUT;
  typedef SOSPolynomial<SOSPolynomial<T3,LABEL3>,LABEL2> IN;
  static void arrange(OUT& to,const IN& from,std::unordered_map<char,SOSInfo> m) {
    for(sizeType t=0; t<(sizeType)from._terms.size(); t++) {
      std::unordered_map<char,SOSInfo> mt=m;
      mt[LABEL2]=from._terms[t];
      Rearrange<T,LABEL,T3,LABEL3>::arrange(to,from._terms[t]._coef,mt);
    }
  }
};
