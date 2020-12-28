//meta-template function to count max number of variables
template <typename T,char LABEL,char LABEL2>
struct NrVar
{
  static sizeType nrVar(const SOSPolynomial<T,LABEL>& s) {
    if(LABEL==LABEL2) {
      sizeType nrVar=0;
      for(sizeType i=0; i<(sizeType)s._terms.size(); i++)
        nrVar=std::max<sizeType>(nrVar,s._terms[i].nrVar());
      return nrVar;
    } else return 0;
  }
  static sizeType order(const SOSPolynomial<T,LABEL>& s) {
    if(LABEL==LABEL2)
      return orderAll(s);
    else return 0;
  }
  static sizeType orderAll(const SOSPolynomial<T,LABEL>& s) {
    sizeType order=0;
    for(sizeType i=0; i<(sizeType)s._terms.size(); i++)
      order=std::max<sizeType>(order,s._terms[i].order());
    return order;
  }
};
template <typename T,char LABEL3,char LABEL,char LABEL2>
struct NrVar<SOSPolynomial<T,LABEL3>,LABEL,LABEL2>
{
  static sizeType nrVar(const SOSPolynomial<SOSPolynomial<T,LABEL3>,LABEL>& s) {
    if(LABEL==LABEL2) {
      sizeType nrVar=0;
      for(sizeType i=0; i<(sizeType)s._terms.size(); i++)
        nrVar=std::max<sizeType>(nrVar,s._terms[i].nrVar());
      return nrVar;
    } else {
      sizeType nrVar=0;
      for(sizeType i=0; i<(sizeType)s._terms.size(); i++)
        nrVar=std::max<sizeType>(nrVar,NrVar<T,LABEL3,LABEL2>::nrVar(s._terms[i]._coef));
      return nrVar;
    }
  }
  static sizeType order(const SOSPolynomial<SOSPolynomial<T,LABEL3>,LABEL>& s) {
    if(LABEL==LABEL2) {
      sizeType order=0;
      for(sizeType i=0; i<(sizeType)s._terms.size(); i++)
        order=std::max<sizeType>(order,s._terms[i].order());
      return order;
    } else {
      sizeType order=0;
      for(sizeType i=0; i<(sizeType)s._terms.size(); i++)
        order=std::max<sizeType>(order,NrVar<T,LABEL3,LABEL2>::order(s._terms[i]._coef));
      return order;
    }
  }
  static sizeType orderAll(const SOSPolynomial<SOSPolynomial<T,LABEL3>,LABEL>& s) {
    sizeType order=0;
    for(sizeType i=0; i<(sizeType)s._terms.size(); i++)
      order=std::max<sizeType>(order,s._terms[i].order()+NrVar<T,LABEL3,LABEL2>::orderAll(s._terms[i]._coef));
    return order;
  }
};
