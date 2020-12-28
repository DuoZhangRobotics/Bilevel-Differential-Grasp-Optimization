//meta-template function to replace index of variable with selected label
template <typename T,char LABEL2>
struct AffineTransXId {
  static T transXId(const T& in,sizeType coef,sizeType off) {
    return in;
  }
};
template <typename T,char LABEL,char LABEL2>
struct AffineTransXId<SOSPolynomial<T,LABEL>,LABEL2> {
  static SOSPolynomial<T,LABEL> transXId(const SOSPolynomial<T,LABEL>& in,sizeType coef,sizeType off) {
    SOSPolynomial<T,LABEL> out=in;
    for(sizeType i=0; i<(sizeType)out._terms.size(); i++) {
      if(LABEL==LABEL2)
        for(sizeType j=0; j<(sizeType)out._terms[i]._id.size(); j++)
          out._terms[i]._id[j]=out._terms[i]._id[j]*coef+off;
      out._terms[i]._coef=AffineTransXId<T,LABEL2>::transXId(out._terms[i]._coef,coef,off);
    }
    return out;
  }
};
