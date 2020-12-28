#ifndef SOS_POLYNOMIAL_H
#define SOS_POLYNOMIAL_H

#include <CommonFile/IO.h>
#include <unordered_map>
#include "ParallelVector.h"

PRJ_BEGIN

template <typename T,char LABEL>
class SOSTerm;
template <typename T,char LABEL>
class SOSPolynomial;
class SOSInfo
{
public:
  std::vector<sizeType> _id,_order;
};
template <typename T2>
struct Zero<SOSPolynomial<T2,'a'>> {
  static T2 value() {
    return SOSPolynomial<T2,'a'>();
  }
  static T2 value(const T2&) {
    return SOSPolynomial<T2,'a'>();
  }
};
#include "internal/SOSPolynomialScalarOfT.h"
#include "internal/SOSPolynomialAffineTransXId.h"
#include "internal/SOSPolynomialContract.h"
#include "internal/SOSPolynomialRearrange.h"
#include "internal/SOSPolynomialCast.h"
#include "internal/SOSPolynomialNrVar.h"
//SOSTerm
template <typename T,char LABEL>
class SOSTerm : public SOSInfo, public SerializableBase
{
public:
  typedef Eigen::Matrix<T,-1,1> VEC;
  typedef Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,sizeType>> STRIPS;
  typedef Eigen::Triplet<T,sizeType> STRIP;
  SOSTerm();
  SOSTerm(T val);
  SOSTerm(T val,sizeType id,sizeType order);
  SOSTerm(const std::string& str);
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  virtual std::string type() const;
  //cmp
  bool operator<(const SOSTerm& other) const;
  bool operator>(const SOSTerm& other) const;
  bool operator==(const SOSTerm& other) const;
  bool operator!=(const SOSTerm& other) const;
  //op
  sizeType nrVar() const;
  sizeType order() const;
  bool hasId(sizeType id) const;
  sizeType order(sizeType id) const;
  SOSTerm removeId(sizeType id) const;
  SOSTerm operator*(T other) const;
  SOSTerm& operator*=(T other);
  SOSTerm operator*(const SOSTerm& other) const;
  SOSTerm& operator*=(const SOSTerm& other);
  SOSTerm operator+(const SOSTerm& other) const;
  SOSTerm& operator+=(const SOSTerm& other);
  SOSTerm operator-(const SOSTerm& other) const;
  SOSTerm& operator-=(const SOSTerm& other);
  SOSPolynomial<T,LABEL> integrate(sizeType nrVar) const;
  void gradient(std::vector<SOSPolynomial<T,LABEL>>& grad) const;
  void gradient(sizeType row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,sizeType>>& grad) const;
  T integrate(const COLD& L,const COLD& U) const;
  template <typename HESS>
  T eval(const COLD& x,VEC* grad=NULL,HESS* hess=NULL) const;
  T evalGradDir(const COLD& x,const COLD& dir) const;
  template <typename JAC>
  T evalJacTpl(sizeType row,const COLD& x,JAC* jac=NULL) const;
  T evalJac(sizeType row,const COLD& x,STRIPS* jac=NULL) const;
  T evalJac(sizeType row,const COLD& x,MAT* jac=NULL) const;
  //misc
  template <typename T2>
  typename Cast<SOSTerm,T2>::Type cast() const {
    typename Cast<SOSTerm,T2>::Type ret;
    ret._coef=Cast<T,T2>::cast(_coef);
    ret._id=_id;
    ret._order=_order;
    return ret;
  }
  //io
  std::string toString(const std::unordered_map<sizeType,std::string>* varNames=NULL) const;
  std::string formattedString(const std::unordered_map<sizeType,std::string>* varNames=NULL) const;
  void readFormattedString(std::string str);
  void add(sizeType id,sizeType order);
  void consistencyCheck() const;
  //data
  T _coef;
};
//SOSPolynomial
template <typename T,char LABEL>
class SOSPolynomial : public SerializableBase
{
public:
  typedef Eigen::Matrix<T,-1,1> VEC;
  typedef Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::SparseMatrix<T,0,sizeType> SMAT;
  typedef Eigen::Matrix<SOSPolynomial<T,LABEL>,-1,1> VECP;
  typedef Eigen::Matrix<SOSPolynomial<T,LABEL>,-1,-1> MATP;
  typedef Eigen::SparseMatrix<SOSPolynomial<T,LABEL>,0,sizeType> SMATP;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,sizeType>> STRIPS;
  typedef Eigen::Triplet<T,sizeType> STRIP;
  SOSPolynomial();
  SOSPolynomial(int other);
  SOSPolynomial(sizeType other);
  SOSPolynomial(typename ScalarOfT<T>::Type other);
  SOSPolynomial(const SOSTerm<T,LABEL>& other);
  SOSPolynomial(const std::string& str);
  template <typename T2,char LABEL2>
  SOSPolynomial(const SOSPolynomial<T2,LABEL2>& other)
  {
    Rearrange<T,LABEL,T2,LABEL2>::arrange(*this,other,std::unordered_map<char,SOSInfo>());
  }
  virtual bool read(std::istream& is);
  virtual bool write(std::ostream& os) const;
  virtual std::string type() const;
  //cmp
  bool operator==(const SOSPolynomial& other) const;
  //op
  sizeType nrVar() const;
  template <char LABEL2>
  sizeType nrVar() const
  {
    return NrVar<T,LABEL,LABEL2>::nrVar(*this);
  }
  sizeType order() const;
  template <char LABEL2>
  sizeType order() const
  {
    return NrVar<T,LABEL,LABEL2>::order(*this);
  }
  sizeType orderAll() const;
  template <char LABEL2>
  sizeType orderAll() const
  {
    return NrVar<T,LABEL,LABEL2>::orderAll(*this);
  }
  SOSPolynomial& operator=(int other);
  SOSPolynomial& operator=(sizeType other);
  SOSPolynomial& operator=(typename ScalarOfT<T>::Type other);
  SOSPolynomial& operator=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial& operator=(const std::string& str);
  SOSPolynomial operator*(const T& other) const;
  SOSPolynomial& operator*=(const T& other);
  SOSPolynomial operator*(const SOSTerm<T,LABEL>& other) const;
  SOSPolynomial& operator*=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial operator*(const SOSPolynomial& other) const;
  SOSPolynomial& operator*=(const SOSPolynomial& other);
  SOSPolynomial operator+(const SOSTerm<T,LABEL>& other) const;
  SOSPolynomial& operator+=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial operator+(const SOSPolynomial& other) const;
  SOSPolynomial& operator+=(const SOSPolynomial& other);
  SOSPolynomial operator-(const SOSTerm<T,LABEL>& other) const;
  SOSPolynomial& operator-=(const SOSTerm<T,LABEL>& other);
  SOSPolynomial operator-(const SOSPolynomial& other) const;
  SOSPolynomial& operator-=(const SOSPolynomial& other);
  SOSPolynomial integrate() const;
  std::vector<SOSPolynomial> gradient() const;
  std::vector<std::vector<SOSPolynomial>> hessian() const;
  VECP gradientV() const;
  MATP hessianM() const;
  void gradientSparse(sizeType row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,sizeType>>& grad) const;
  SMATP hessianSparse() const;
  void sum(const std::vector<SOSPolynomial,Eigen::aligned_allocator<SOSPolynomial>>& polys);
  VEC gradientCoef() const;
  MAT hessianCoef() const;
  MAT JTJCoef() const;
  T integrate(const COLD& L,const COLD& U) const;
  T eval(const COLD& x,VEC* grad=NULL,MAT* hess=NULL) const;
  T evalTrips(const COLD& x,VEC* grad=NULL,STRIPS* hess=NULL) const;
  T evalGradDir(const COLD& x,const COLD& dir) const;
  template <typename JAC>
  T evalJac(sizeType row,const COLD& x,JAC* jac=NULL) const
  {
    T ret=ScalarOfT<T>::convert(0);
    for(sizeType i=0; i<(sizeType)_terms.size(); i++)
      ret+=_terms[i].evalJac(row,x,jac);
    return ret;
  }
  template <typename JAC>
  VEC evalJac(const COLD& x,JAC* jac=NULL) const
  {
    sizeType nrT=(sizeType)_terms.size();
    VEC vec=VEC::Constant(nrT,ScalarOfT<T>::convert(0));
    for(sizeType row=0; row<nrT; row++)
      vec[row]=_terms[row].evalJac(row,x,jac);
    return vec;
  }
  SOSPolynomial rename(const std::vector<sizeType>& ids) const;
  SOSPolynomial varRemap(const std::vector<sizeType>& ids) const;
  SOSPolynomial setAllCoef(T coef) const;
  SOSPolynomial removeZero(typename ScalarOfT<T>::Type eps) const;
  SOSPolynomial removeVariableId(sizeType id) const;
  SOSPolynomial linearConstraint(sizeType id,const SOSPolynomial& cons) const;
  SOSPolynomial linearTransform(const std::unordered_map<sizeType,typename ScalarOfT<T>::Type>& cons) const;
  SOSPolynomial linearTransform(const std::unordered_map<sizeType,SOSPolynomial>& cons,sizeType reportInterval=0) const;
  static COLD solve(const std::vector<SOSPolynomial>& LHS,const std::vector<SOSPolynomial>& RHS);
  //misc
  static SMAT hessianGradDir(const SMATP& h,const COLD& x,const COLD& dir)
  {
    STRIPS trips;
    for(sizeType k=0; k<h.outerSize(); ++k)
      for(typename SMATP::InnerIterator it(h,k); it; ++it)
        trips.push_back(STRIP(it.row(),it.col(),it.value().evalGradDir(x,dir)));

    SMAT ret;
    ret.resize(h.rows(),h.cols());
    ret.setFromTriplets(trips.begin(),trips.end());
    return ret;
  }
  template <char LABEL2>
  typename Contract<SOSPolynomial,LABEL2>::Type eval(const COLD& x) const
  {
    return Contract<SOSPolynomial,LABEL2>::eval(*this,x);
  }
  template <char LABEL2>
  SOSPolynomial affineTransXId(sizeType coef,sizeType off) const
  {
    return AffineTransXId<SOSPolynomial,LABEL2>::transXId(*this,coef,off);
  }
  template <typename T2>
  typename Cast<SOSPolynomial,T2>::Type cast() const {
    typename Cast<SOSPolynomial,T2>::Type ret;
    ret._terms.resize(_terms.size());
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_terms.size(); i++)
      ret._terms[i]=_terms[i].template cast<T2>();
    return ret;
  }
  //io
  static bool read(std::vector<VEC>& v,std::istream& is);
  static bool read(std::vector<MAT>& m,std::istream& is);
  static bool write(const std::vector<VEC>& v,std::ostream& os);
  static bool write(const std::vector<MAT>& m,std::ostream& os);
  static bool read(VEC& v,std::istream& is);
  static bool read(MAT& m,std::istream& is);
  static bool write(const VEC& v,std::ostream& os);
  static bool write(const MAT& m,std::ostream& os);
  operator T() const;
  //string
  std::string toString(const std::unordered_map<sizeType,std::string>* varNames=NULL) const;
  std::string formattedString(const std::unordered_map<sizeType,std::string>* varNames=NULL) const;
  void readFormattedString(std::string str);
  void operator<<(std::string str);
  //debug
  void debugIntegrate() const;
  void debugGradient() const;
  void debugEval() const;
  void add(const SOSTerm<T,LABEL>& other);
  void consistencyCheck() const;
  void makeConsistent();
  //data
  std::vector<SOSTerm<T,LABEL>> _terms;
};
//common type
#define DECL_POLYTYPES(TT,POSTFIX)  \
typedef SOSTerm<TT,'x'> TermX##POSTFIX;  \
typedef SOSPolynomial<TT,'x'> PolyX##POSTFIX;  \
typedef SOSTerm<TT,'t'> TermT##POSTFIX;  \
typedef SOSPolynomial<TT,'t'> PolyT##POSTFIX;  \
typedef SOSTerm<TT,'a'> TermA##POSTFIX;  \
typedef SOSPolynomial<TT,'a'> PolyA##POSTFIX;  \
typedef SOSTerm<PolyX##POSTFIX,'a'> TermXA##POSTFIX;  \
typedef SOSPolynomial<PolyX##POSTFIX,'a'> PolyXA##POSTFIX;  \
typedef SOSTerm<PolyA##POSTFIX,'x'> TermAX##POSTFIX;  \
typedef SOSPolynomial<PolyA##POSTFIX,'x'> PolyAX##POSTFIX;  \
typedef SOSTerm<PolyT##POSTFIX,'x'> TermTX##POSTFIX;  \
typedef SOSPolynomial<PolyT##POSTFIX,'x'> PolyTX##POSTFIX;  \
typedef SOSTerm<PolyTX##POSTFIX,'a'> TermTXA##POSTFIX;  \
typedef SOSPolynomial<PolyTX##POSTFIX,'a'> PolyTXA##POSTFIX;

PRJ_END

namespace Eigen
{
namespace internal
{
template<typename T>
struct scalar_product_traits<COMMON::SOSPolynomial<T,'a'>,T>
{
  enum {
    Defined = 1
  };
  typedef COMMON::SOSPolynomial<T,'a'> ReturnType;
};
template<typename T>
struct scalar_product_traits<T,COMMON::SOSPolynomial<T,'a'>>
{
  enum {
    Defined = 1
  };
  typedef COMMON::SOSPolynomial<T,'a'> ReturnType;
};
}
}

#endif
