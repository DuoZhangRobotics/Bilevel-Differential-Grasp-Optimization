#include "SOSPolynomial.h"
#include "SMatNumber.h"
#include "DebugGradient.h"
#include "SparseUtils.h"
#include "Utils.h"

USE_PRJ_NAMESPACE

#define NO_CONSISTENCY_CHECK
#include "internal/SOSPolynomialLexer.h"
#include "internal/SOSPolynomialSolve.h"
#include "internal/SOSPolynomialConvert.h"
#include "internal/SOSPolynomialEvaluate.h"
#include "internal/SOSPolynomialRobustInversion.h"
#include "internal/SOSPolynomialIsZero.h"
//SOSTerm
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm() {}
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm(T val):_coef(val) {}
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm(T val,sizeType id,sizeType order):_coef(val) {
  if(order>0) {
    _id.assign(1,id);
    _order.assign(1,order);
  }
  consistencyCheck();
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>::SOSTerm(const std::string& str) {
  readFormattedString(str);
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::read(std::istream& is) {
  readBinaryData(_id,is);
  readBinaryData(_order,is);
  scalarD coef;
  readBinaryData(coef,is);
  _coef=coef;
  return is.good();
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::write(std::ostream& os) const {
  writeBinaryData(_id,os);
  writeBinaryData(_order,os);
  scalarD coef=std::to_double(_coef);
  writeBinaryData(coef,os);
  return os.good();
}
template <typename T,char LABEL>
std::string SOSTerm<T,LABEL>::type() const {
  return typeid(SOSTerm).name();
}
//cmp
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator<(const SOSTerm& other) const {
  sizeType i=0;
  while(i<(sizeType)_id.size()&&i<(sizeType)other._id.size()) {
    if(_id[i]<other._id[i])
      return true;
    else if(_id[i]>other._id[i])
      return false;
    else {
      if(_order[i]<other._order[i])
        return true;
      else if(_order[i]>other._order[i])
        return false;
    }
    i++;
  }
  if(_id.size()<other._id.size())
    return true;
  else if(_id.size()>other._id.size())
    return false;
  else return false;
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator>(const SOSTerm& other) const {
  return other<*this;
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator==(const SOSTerm& other) const {
  return !(*this<other) && !(other<*this);
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::operator!=(const SOSTerm& other) const {
  return !(*this==other);
}
//op
template <typename T,char LABEL>
sizeType SOSTerm<T,LABEL>::nrVar() const
{
  if(_id.empty())
    return 0;
  else return _id.back()+1;
}
template <typename T,char LABEL>
sizeType SOSTerm<T,LABEL>::order() const
{
  if(_order.empty())
    return 0;
  else {
    sizeType order=0;
    for(sizeType i=0; i<(sizeType)_order.size(); i++)
      order+=_order[i];
    return order;
  }
}
template <typename T,char LABEL>
bool SOSTerm<T,LABEL>::hasId(sizeType id) const
{
  std::vector<sizeType>::const_iterator it=std::lower_bound(_id.begin(),_id.end(),id);
  return it!=_id.end() && *it==id;
}
template <typename T,char LABEL>
sizeType SOSTerm<T,LABEL>::order(sizeType id) const
{
  std::vector<sizeType>::const_iterator it=std::lower_bound(_id.begin(),_id.end(),id);
  sizeType off=it-_id.begin();
  return _order[off];
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::removeId(sizeType id) const
{
  SOSTerm ret=*this;
  std::vector<sizeType>::const_iterator it=std::lower_bound(ret._id.begin(),ret._id.end(),id);
  if(it!=ret._id.end() && *it==id) {
    sizeType off=it-ret._id.begin();
    ret._id.erase(ret._id.begin()+off);
    ret._order.erase(ret._order.begin()+off);
  }
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator*(T other) const {
  SOSTerm<T,LABEL> ret=*this;
  ret._coef*=other;
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator*=(T other) {
  return *this=*this*other;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator*(const SOSTerm& other) const {
  SOSTerm<T,LABEL> ret;
  sizeType i=0,j=0;
  while(i<(sizeType)_id.size() && j<(sizeType)other._id.size()) {
    if(_id[i]<other._id[j]) {
      ret.add(_id[i],_order[i]);
      i++;
    } else if(_id[i]>other._id[j]) {
      ret.add(other._id[j],other._order[j]);
      j++;
    } else {
      ret.add(_id[i],_order[i]+other._order[j]);
      i++;
      j++;
    }
  }
  while(i<(sizeType)_id.size()) {
    ret.add(_id[i],_order[i]);
    i++;
  }
  while(j<(sizeType)other._id.size()) {
    ret.add(other._id[j],other._order[j]);
    j++;
  }
  ret._coef=_coef*other._coef;
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator*=(const SOSTerm& other) {
  return *this=*this*other;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator+(const SOSTerm& other) const {
  ASSERT(*this==other)
  SOSTerm<T,LABEL> ret=*this;
  ret._coef+=other._coef;
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator+=(const SOSTerm& other) {
  return *this=*this+other;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL> SOSTerm<T,LABEL>::operator-(const SOSTerm& other) const {
  ASSERT(*this==other)
  SOSTerm<T,LABEL> ret=*this;
  ret._coef-=other._coef;
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSTerm<T,LABEL>& SOSTerm<T,LABEL>::operator-=(const SOSTerm& other) {
  return *this=*this-other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSTerm<T,LABEL>::integrate(sizeType nrVar) const
{
  sizeType j=0;
  SOSPolynomial<T,LABEL> ret,p;
  for(sizeType i=0; i<nrVar; i++) {
    sizeType order=0;
    if(j<(sizeType)_id.size()&&_id[j]==i)
      order=_order[j++];
    p =SOSTerm(ScalarOfT<T>::convert(1/scalarD(order+1)),i+nrVar,order+1);
    p-=SOSTerm(ScalarOfT<T>::convert(1/scalarD(order+1)),i      ,order+1);
    if(i==0)
      ret=p;
    else ret*=p;
  }
  ASSERT(j==(sizeType)_id.size())
  return ret*SOSTerm(_coef);
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::gradient(std::vector<SOSPolynomial<T,LABEL>>& grad) const
{
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    SOSTerm<T,LABEL> t=*this;
    if(_id.empty())
      continue;
    else if(_order[i]>1) {
      t._order[i]--;
      t._coef*=ScalarOfT<T>::convert(_order[i]);
    } else {
      t._id.erase(t._id.begin()+i);
      t._order.erase(t._order.begin()+i);
    }
    grad[_id[i]]+=t;
  }
  for(sizeType i=0; i<(sizeType)grad.size(); i++)
    grad[i].consistencyCheck();
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::gradient(sizeType row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,sizeType>>& grad) const
{
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    SOSTerm<T,LABEL> t=*this;
    if(_id.empty())
      continue;
    else if(_order[i]>1) {
      t._order[i]--;
      t._coef*=ScalarOfT<T>::convert(_order[i]);
    } else {
      t._id.erase(t._id.begin()+i);
      t._order.erase(t._order.begin()+i);
    }
    grad.push_back(Eigen::Triplet<SOSPolynomial<T,LABEL>,sizeType>(row,_id[i],t));
    //grad[_id[i]]+=t;
  }
  //for(sizeType i=0; i<(sizeType)grad.size(); i++)
  //  grad[i].consistencyCheck();
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::integrate(const COLD& L,const COLD& U) const
{
  T ret=ScalarOfT<T>::convert(1);
  sizeType j=0;
  for(sizeType i=0; i<L.size(); i++) {
    sizeType order=0;
    if(j<(sizeType)_id.size()&&_id[j]==i)
      order=_order[j++];
    ret*=ScalarOfT<T>::convert((std::pow(U[i],order+1)-std::pow(L[i],order+1))/(order+1));
  }
  ASSERT(j==(sizeType)_id.size())
  return _coef*ret;
}
template <typename T,char LABEL>
template <typename HESS>
T SOSTerm<T,LABEL>::eval(const COLD& x,VEC* grad,HESS* hess) const
{
  T ret=ScalarOfT<T>::convert(1),fgradX,fhessX,tmp;
  std::vector<T> xx0(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx1(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx2(_id.size(),ScalarOfT<T>::convert(1));
  //value
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    xx0[i]=ScalarOfT<T>::convert(std::pow(x[_id[i]],_order[i]));
    ret*=xx0[i];
    if(i>0)
      xx1[i]=xx1[i-1]*xx0[i-1];
  }
  //gradient
  if(grad)
    for(sizeType i=(sizeType)_id.size()-1; i>=0; i--) {
      fgradX=ScalarOfT<T>::convert(_order[i]*std::pow(x[_id[i]],_order[i]-1));
      ParallelEvaluate<T>::gradient(*grad,_id[i],_coef*fgradX*xx1[i]*xx2[i]);
      if(i>0)
        xx2[i-1]=xx2[i]*xx0[i];
    }
  //hessian
  if(hess)
    for(sizeType i=0; i<(sizeType)_id.size(); i++) {
      if(_order[i]>1) {
        fhessX=ScalarOfT<T>::convert(_order[i]*(_order[i]-1)*std::pow(x[_id[i]],_order[i]-2));
        ParallelEvaluate<T>::hessian(*hess,_id[i],_id[i],_coef*fhessX*xx1[i]*xx2[i]);
      }
      fgradX=ScalarOfT<T>::convert(_order[i]*std::pow(x[_id[i]],_order[i]-1));
      tmp=xx1[i]*fgradX;
      for(sizeType j=i+1; j<(sizeType)_id.size(); j++) {
        fgradX=ScalarOfT<T>::convert(_order[j]*std::pow(x[_id[j]],_order[j]-1));
        ParallelEvaluate<T>::hessianSym(*hess,_id[i],_id[j],_coef*tmp*fgradX*xx2[j]);
        tmp*=xx0[j];
      }
    }
  return _coef*ret;
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::evalGradDir(const COLD& x,const COLD& dir) const
{
  T ret=ScalarOfT<T>::convert(0),fgradX;
  std::vector<T> xx0(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx1(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx2(_id.size(),ScalarOfT<T>::convert(1));
  //value
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    xx0[i]=ScalarOfT<T>::convert(std::pow(x[_id[i]],_order[i]));
    ret*=xx0[i];
    if(i>0)
      xx1[i]=xx1[i-1]*xx0[i-1];
  }
  //gradient
  for(sizeType i=(sizeType)_id.size()-1; i>=0; i--) {
    fgradX=ScalarOfT<T>::convert(_order[i]*std::pow(x[_id[i]],_order[i]-1));
    ret+=ScalarOfT<T>::convert(dir[_id[i]])*_coef*fgradX*xx1[i]*xx2[i];
    if(i>0)
      xx2[i-1]=xx2[i]*xx0[i];
  }
  return ret;
}
template <typename T,char LABEL>
template <typename JAC>
T SOSTerm<T,LABEL>::evalJacTpl(sizeType row,const COLD& x,JAC* jac) const
{
  T ret=ScalarOfT<T>::convert(1),fgradX;
  std::vector<T> xx0(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx1(_id.size(),ScalarOfT<T>::convert(1));
  std::vector<T> xx2(_id.size(),ScalarOfT<T>::convert(1));
  //value
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    xx0[i]=ScalarOfT<T>::convert(std::pow(x[_id[i]],_order[i]));
    ret*=xx0[i];
    if(i>0)
      xx1[i]=xx1[i-1]*xx0[i-1];
  }
  //gradient
  if(jac)
    for(sizeType i=(sizeType)_id.size()-1; i>=0; i--) {
      fgradX=ScalarOfT<T>::convert(_order[i]*std::pow(x[_id[i]],_order[i]-1));
      ParallelEvaluate<T>::hessian(*jac,row,_id[i],_coef*fgradX*xx1[i]*xx2[i]);
      if(i>0)
        xx2[i-1]=xx2[i]*xx0[i];
    }
  return _coef*ret;
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::evalJac(sizeType row,const COLD& x,STRIPS* jac) const
{
  return evalJacTpl(row,x,jac);
}
template <typename T,char LABEL>
T SOSTerm<T,LABEL>::evalJac(sizeType row,const COLD& x,MAT* jac) const
{
  return evalJacTpl(row,x,jac);
}
//io
template <typename T,char LABEL>
std::string SOSTerm<T,LABEL>::toString(const std::unordered_map<sizeType,std::string>* varNames) const {
  bool bracketNeeded=false;
  std::string ret=convert(_coef,bracketNeeded);
  if(bracketNeeded)
    ret="("+ret+")";
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    if(varNames)
      ret+="*"+varNames->find(_id[i])->second;
    else ret+="*"+std::string(1,LABEL)+convert(_id[i],bracketNeeded);
    if(_order[i]>1)
      ret+="^"+convert(_order[i],bracketNeeded);
  }
  return ret;
}
template <typename T,char LABEL>
std::string SOSTerm<T,LABEL>::formattedString(const std::unordered_map<sizeType,std::string>* varNames) const {
  std::string ret="("+convertFormatted(_coef)+")";
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    if(varNames)
      ret+="*"+varNames->find(_id[i])->second;
    else ret+="*"+std::string(1,LABEL)+convertFormatted(_id[i]);
    ret+="^"+convertFormatted(_order[i]);
  }
  return ret;
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::readFormattedString(std::string str) {
  str.erase(std::remove(str.begin(),str.end(),' '),str.end());
  str.erase(std::remove(str.begin(),str.end(),'\n'),str.end());
  _id.clear();
  _order.clear();

  //coef part
  sizeType j=str.size();
  while(j>0 && str[j-1]!=')')
    j--;
  if(j>0)
    convertFormatted(str.substr(1,j-2),_coef);
  else {
    size_t pos=str.find_first_of('*');
    if(pos==std::string::npos)
      j=str.size();
    else j=(sizeType)pos;
    convertFormatted(str.substr(0,j),_coef);
  }
  if(j==(sizeType)str.size())
    return;

  //other part
  std::istringstream iss(str.substr(j));
  char m,c,o;
  sizeType id,order;
  while(!iss.eof()) {
    iss >> m >> c >> id >> o >> order;
    ASSERT(m=='*' && c==LABEL && o=='^')
    add(id,order);
  }
}
//helper
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::add(sizeType id,sizeType order) {
  std::vector<sizeType>::iterator it=std::lower_bound(_id.begin(),_id.end(),id);
  sizeType off=it-_id.begin();
  if(it==_id.end() || *it!=id) {
    _id.insert(off+_id.begin(),id);
    _order.insert(off+_order.begin(),order);
  } else {
    _order[off]+=order;
  }
}
template <typename T,char LABEL>
void SOSTerm<T,LABEL>::consistencyCheck() const {
#ifndef NO_CONSISTENCY_CHECK
  for(sizeType i=0; i<(sizeType)_id.size()-1; i++) {
    ASSERT(_id[i]<_id[i+1])
  }
  for(sizeType i=0; i<(sizeType)_id.size(); i++) {
    ASSERT(_order[i]>0)
  }
#endif
}

//SOSPolynomial
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial() {}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(int other)
{
  *this=(typename ScalarOfT<T>::Type)other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(sizeType other)
{
  *this=(typename ScalarOfT<T>::Type)other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(typename ScalarOfT<T>::Type other)
{
  *this=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(const SOSTerm<T,LABEL>& other):_terms(1,other) {}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::SOSPolynomial(const std::string& str) {
  *this<<str;
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::read(std::istream& is) {
  readBinaryData(_terms,is);
  return is.good();
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::write(std::ostream& os) const {
  writeBinaryData(_terms,os);
  return os.good();
}
template <typename T,char LABEL>
std::string SOSPolynomial<T,LABEL>::type() const {
  return typeid(SOSPolynomial).name();
}
//cmp
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::operator==(const SOSPolynomial& other) const
{
  return false;
}
//op
template <typename T,char LABEL>
sizeType SOSPolynomial<T,LABEL>::nrVar() const
{
  return nrVar<LABEL>();
}
template <typename T,char LABEL>
sizeType SOSPolynomial<T,LABEL>::order() const
{
  return order<LABEL>();
}
template <typename T,char LABEL>
sizeType SOSPolynomial<T,LABEL>::orderAll() const
{
  return orderAll<LABEL>();
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(int other)
{
  *this=(typename ScalarOfT<T>::Type)other;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(sizeType other)
{
  *this=(typename ScalarOfT<T>::Type)other;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(typename ScalarOfT<T>::Type other)
{
  *this=ScalarOfT<SOSPolynomial<T,LABEL>>::convert(other);
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(const SOSTerm<T,LABEL>& other)
{
  _terms.assign(1,other);
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator=(const std::string& str)
{
  *this<<str;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator*(const T& other) const {
  SOSPolynomial ret=*this;
  return ret*=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator*=(const T& other) {
  for(sizeType i=0; i<(sizeType)_terms.size(); i++)
    _terms[i]._coef*=other;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator*(const SOSTerm<T,LABEL>& other) const {
  SOSPolynomial ret=*this;
  return ret*=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator*=(const SOSTerm<T,LABEL>& other) {
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)_terms.size(); i++)
    _terms[i]*=other;
  std::sort(_terms.begin(),_terms.end());
  consistencyCheck();
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator*(const SOSPolynomial& other) const {
  SOSPolynomial ret=*this;
  return ret*=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator*=(const SOSPolynomial& other) {
  std::vector<SOSPolynomial,Eigen::aligned_allocator<SOSPolynomial>> otherByTerms(other._terms.size());
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)other._terms.size(); i++)
    otherByTerms[i]=*this*other._terms[i];
  sum(otherByTerms);
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator+(const SOSTerm<T,LABEL>& other) const {
  SOSPolynomial ret=*this;
  return ret+=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator+=(const SOSTerm<T,LABEL>& other) {
  typename std::vector<SOSTerm<T,LABEL>>::iterator it=std::lower_bound(_terms.begin(),_terms.end(),other);
  if(it==_terms.end() || *it!=other)
    _terms.insert(it,other);
  else it->_coef+=other._coef;
  return *this;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator+(const SOSPolynomial& other) const {
  //degenerate case
  if(_terms.empty())
    return other;
  else if(other._terms.empty())
    return *this;
  else if(other._terms.size()==1)
    return operator+(other._terms[0]);
  else if(_terms.size()==1)
    return other.operator+(_terms[0]);
  //merge sort
  sizeType i=0,j=0;
  SOSPolynomial ret;
  while(i<(sizeType)_terms.size() && j<(sizeType)other._terms.size()) {
    if(_terms[i]<other._terms[j]) {
      ret.add(_terms[i]);
      i++;
    } else if(_terms[i]>other._terms[j]) {
      ret.add(other._terms[j]);
      j++;
    } else {
      ret.add(_terms[i]+other._terms[j]);
      i++;
      j++;
    }
  }
  while(i<(sizeType)_terms.size()) {
    ret.add(_terms[i]);
    i++;
  }
  while(j<(sizeType)other._terms.size()) {
    ret.add(other._terms[j]);
    j++;
  }
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator+=(const SOSPolynomial& other) {
  return *this=*this+other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator-(const SOSTerm<T,LABEL>& other) const {
  SOSPolynomial ret=*this;
  return ret-=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator-=(const SOSTerm<T,LABEL>& other) {
  SOSTerm<T,LABEL> otherNeg=other;
  otherNeg._coef*=ScalarOfT<T>::convert(-1);
  return operator+=(otherNeg);
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::operator-(const SOSPolynomial& other) const {
  SOSPolynomial ret=*this;
  return ret-=other;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>& SOSPolynomial<T,LABEL>::operator-=(const SOSPolynomial& other) {
  return operator+=(other*ScalarOfT<T>::convert(-1));
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::integrate() const
{
  if(_terms.empty())
    return SOSPolynomial();
  SOSPolynomial ret=_terms[0].integrate(nrVar());
  for(sizeType i=1; i<(sizeType)_terms.size(); i++)
    ret+=_terms[i].integrate(nrVar());
  return ret;
}
template <typename T,char LABEL>
std::vector<SOSPolynomial<T,LABEL>> SOSPolynomial<T,LABEL>::gradient() const
{
  std::vector<SOSPolynomial<T,LABEL>> ret(nrVar());
  for(sizeType i=0; i<(sizeType)_terms.size(); i++)
    _terms[i].gradient(ret);
  return ret;
}
template <typename T,char LABEL>
std::vector<std::vector<SOSPolynomial<T,LABEL>>> SOSPolynomial<T,LABEL>::hessian() const
{
  std::vector<std::vector<SOSPolynomial>> ret;
  std::vector<SOSPolynomial<T,LABEL>> g=gradient();
  for(sizeType i=0; i<(sizeType)g.size(); i++)
    ret.push_back(g[i].gradient());
  return ret;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::VECP SOSPolynomial<T,LABEL>::gradientV() const
{
  VECP ret;
  ret.resize(nrVar());//=VECP::Zero(nrVar());
  std::vector<SOSPolynomial<T,LABEL>> retV=gradient();
  for(sizeType i=0; i<ret.size(); i++)
    ret[i]=retV[i];
  return ret;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::MATP SOSPolynomial<T,LABEL>::hessianM() const
{
  MATP ret;
  ret.resize(nrVar(),nrVar());//=MATP::Zero(nrVar(),nrVar());
  std::vector<std::vector<SOSPolynomial<T,LABEL>>> hessV=hessian();
  for(sizeType i=0; i<(sizeType)hessV.size(); i++)
    for(sizeType j=0; j<(sizeType)hessV[i].size(); j++)
      ret(i,j)=hessV[i][j];
  return ret;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::gradientSparse(sizeType row,ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,sizeType>>& grad) const
{
  for(sizeType i=0; i<(sizeType)_terms.size(); i++)
    _terms[i].gradient(row,grad);
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::SMATP SOSPolynomial<T,LABEL>::hessianSparse() const
{
  std::vector<SOSPolynomial<T,LABEL>> g=gradient();
  ParallelVector<Eigen::Triplet<SOSPolynomial<T,LABEL>,sizeType>> trips;
  for(sizeType i=0; i<(sizeType)g.size(); i++)
    g[i].gradientSparse(i,trips);

  SMATP ret;
  ret.resize(nrVar(),nrVar());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::sum(const std::vector<SOSPolynomial,Eigen::aligned_allocator<SOSPolynomial>>& polys)
{
  //merge
  _terms.clear();
  for(sizeType i=0; i<(sizeType)polys.size(); i++)
    _terms.insert(_terms.end(),polys[i]._terms.begin(),polys[i]._terms.end());
  if(_terms.empty())
    return;
  makeConsistent();
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::VEC SOSPolynomial<T,LABEL>::gradientCoef() const
{
  VEC gc=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
  std::vector<SOSPolynomial<T,LABEL>> g=gradient();
  for(sizeType i=0; i<nrVar(); i++)
    gc[i]=g.at(i).operator T();
  return gc;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::MAT SOSPolynomial<T,LABEL>::hessianCoef() const
{
  MAT hc=MAT::Constant(nrVar(),nrVar(),ScalarOfT<T>::convert(0));
  std::vector<std::vector<SOSPolynomial>> h=hessian();
  for(sizeType i=0; i<nrVar(); i++)
    for(sizeType j=0; j<(sizeType)h.at(i).size(); j++)
      hc(i,j)=h.at(i).at(j).operator T();
  return hc;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::MAT SOSPolynomial<T,LABEL>::JTJCoef() const
{
  sizeType nrT=(sizeType)_terms.size();
  MAT JTJ=MAT::Constant(nrT,nrT,ScalarOfT<T>::convert(0));
  OMP_PARALLEL_FOR_
  for(sizeType r=0; r<nrT; r++)
    for(sizeType c=0; c<nrT; c++)
      JTJ(r,c)=_terms[r]._coef*_terms[c]._coef;
  return JTJ;
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::integrate(const COLD& L,const COLD& U) const
{
  if(_terms.empty())
    return T();
  std::vector<T> termsI(_terms.size()),tmp;
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)_terms.size(); i++)
    termsI[i]=_terms[i].integrate(L,U);
  while(termsI.size()>1) {
    tmp.resize((termsI.size()+1)/2);
    OMP_PARALLEL_FOR_
    for(sizeType k=0; k<(sizeType)tmp.size(); k++)
      if(k*2==(sizeType)termsI.size()-1)
        tmp[k]=termsI[k*2];
      else tmp[k]=termsI[k*2]+termsI[k*2+1];
    termsI.swap(tmp);
  }
  return termsI[0];
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::eval(const COLD& x,VEC* grad,MAT* hess) const
{
  return ParallelEvaluate<T>::eval(*this,x,grad,hess);
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::evalTrips(const COLD& x,VEC* grad,STRIPS* hess) const
{
  return ParallelEvaluate<T>::eval(*this,x,grad,hess);
}
template <typename T,char LABEL>
T SOSPolynomial<T,LABEL>::evalGradDir(const COLD& x,const COLD& dir) const
{
  T ret=ScalarOfT<T>::convert(0);
  for(sizeType i=0; i<(sizeType)_terms.size(); i++)
    ret+=_terms[i].evalGradDir(x,dir);
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::rename(const std::vector<sizeType>& ids) const
{
  //collect ids
  std::unordered_map<sizeType,sizeType> idMap;
  for(sizeType i=0; i<(sizeType)ids.size(); i++)
    idMap[i]=ids[i];
  //remap
  SOSPolynomial ret=*this;
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)ret._terms.size(); i++) {
    const SOSTerm<T,LABEL>& t0=_terms[i];
    SOSTerm<T,LABEL>& t=ret._terms[i];
    t._id.clear();
    t._order.clear();
    for(sizeType d=0; d<(sizeType)t0._id.size(); d++)
      t.add(idMap.find(t0._id[d])->second,t0._order[d]);
  }
  std::sort(ret._terms.begin(),ret._terms.end());
  ret.consistencyCheck();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::varRemap(const std::vector<sizeType>& ids) const
{
  //remap
  ParallelVector<SOSTerm<T,LABEL>> terms;
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)_terms.size(); i++) {
    const SOSTerm<T,LABEL>& term=_terms[i];
    SOSTerm<T,LABEL> termRet(term._coef);
    for(sizeType d=0; d<(sizeType)term._id.size(); d++) {
      sizeType newId=ids.at(term._id[d]);
      ASSERT(newId>=0)
      termRet.add(newId,term._order[d]);
    }
    terms.push_back(termRet);
  }
  SOSPolynomial ret;
  ret._terms.insert(ret._terms.end(),terms.begin(),terms.end());
  ret.makeConsistent();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::setAllCoef(T coef) const
{
  SOSPolynomial ret=*this;
  for(sizeType i=0; i<(sizeType)ret._terms.size(); i++)
    ret._terms[i]._coef=coef;
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::removeZero(typename ScalarOfT<T>::Type eps) const
{
  return IsZero<SOSPolynomial<T,LABEL>>::removeZero(*this,eps);
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::removeVariableId(sizeType id) const
{
  SOSPolynomial ret=*this;
  for(sizeType i=0; i<(sizeType)ret._terms.size(); i++) {
    SOSTerm<T,LABEL>& t=ret._terms[i];
    for(sizeType j=0; j<(sizeType)t._id.size(); j++) {
      ASSERT(t._id[j]!=id)
      if(t._id[j]>id)
        t._id[j]--;
    }
  }
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::linearConstraint(sizeType id,const SOSPolynomial& cons) const
{
  SOSPolynomial ret;
  for(sizeType i=0; i<(sizeType)_terms.size(); i++) {
    SOSTerm<T,LABEL> t=_terms[i];
    if(t.hasId(id)) {
      sizeType order=t.order(id);
      SOSTerm<T,LABEL> t2=t.removeId(id);
      SOSPolynomial p=cons;
      for(sizeType i=1; i<order; i++)
        p*=cons;
      ret+=p*t2;
    } else {
      ret+=t;
    }
  }
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::linearTransform(const std::unordered_map<sizeType,typename ScalarOfT<T>::Type>& cons) const
{
  ParallelVector<SOSTerm<T,LABEL>> terms;
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<(sizeType)_terms.size(); i++) {
    const SOSTerm<T,LABEL>& term=_terms[i];
    SOSTerm<T,LABEL> termRet(term._coef);
    for(sizeType v=0; v<(sizeType)term._id.size(); v++) {
      typename std::unordered_map<sizeType,typename ScalarOfT<T>::Type>::const_iterator it=cons.find(term._id[v]);
      if(it!=cons.end()) {
        termRet._coef*=ScalarOfT<T>::convert(std::pow(it->second,term._order[v]));
      } else {
        termRet._id.push_back(term._id[v]);
        termRet._order.push_back(term._order[v]);
      }
    }
    terms.push_back(termRet);
  }
  SOSPolynomial ret;
  ret._terms.insert(ret._terms.end(),terms.begin(),terms.end());
  ret.makeConsistent();
  return ret;
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL> SOSPolynomial<T,LABEL>::linearTransform(const std::unordered_map<sizeType,SOSPolynomial>& cons,sizeType reportInterval) const
{
  SOSPolynomial ret;
  std::vector<std::vector<SOSPolynomial>> cache(nrVar(),std::vector<SOSPolynomial>(1,ScalarOfT<SOSPolynomial>::convert(1)));
  for(sizeType i=0; i<(sizeType)_terms.size(); i++) {
    if(reportInterval>0 && (i%reportInterval)==0) {
      INFOV("linearTransform: %ld/%ld",i,(sizeType)_terms.size())
    }
    const SOSTerm<T,LABEL>& t=_terms[i];
    SOSPolynomial sRet=SOSTerm<T,LABEL>(t._coef);
    for(sizeType d=0; d<(sizeType)t._id.size(); d++) {
      std::vector<SOSPolynomial>& cacheId=cache[t._id[d]];
      if((sizeType)cacheId.size()<=t._order[d]) {
        SOSPolynomial tmpP=cacheId.back();
        const SOSPolynomial& p=cons.find(t._id[d])->second;
        for(sizeType o=(sizeType)cacheId.size()-1; o<(sizeType)t._order[d]; o++) {
          tmpP*=p;
          cacheId.push_back(tmpP);
        }
      }
      sRet*=cacheId[t._order[d]];
    }
    ret+=sRet;
  }
  return ret;
}
template <typename T,char LABEL>
typename SOSPolynomial<T,LABEL>::COLD SOSPolynomial<T,LABEL>::solve(const std::vector<SOSPolynomial>& LHS,const std::vector<SOSPolynomial>& RHS)
{
  sizeType nrEQ=0;
  ParallelVector<Eigen::Triplet<typename ScalarOfT<T>::Type,sizeType>> Lss,Rss;
  for(sizeType i=0; i<(sizeType)LHS.size(); i++)
    Solve<SOSPolynomial>::solve(nrEQ,Lss,Rss,LHS[i],RHS[i]);

  Eigen::SparseMatrix<typename ScalarOfT<T>::Type,0,sizeType> LHSM,RHSM;
  LHSM.resize(nrEQ,nrEQ);
  LHSM.setFromTriplets(Lss.begin(),Lss.end());
  RHSM.resize(nrEQ,1);
  RHSM.setFromTriplets(Rss.begin(),Rss.end());

  COLD RHSD=RHSM.toDense();
  MATD LHSD=LHSM.toDense();
  //std::cout << LHSD << std::endl << RHSD << std::endl;
  //analyze nullspace rows
  sizeType nrRow=0;
  for(sizeType i=0; i<LHSD.rows(); i++)
    if(LHSD.row(i).unaryExpr([&](const T& in) {
    return (T)std::abs(in);
    }).maxCoeff()!=0.0) {
    LHSD.row(nrRow)=LHSD.row(i);
    RHSD[nrRow++]=RHSD[i];
  }
  else {
    ASSERT(RHSD[i]==0.0)
  }
  LHSD=LHSD.block(0,0,nrRow,nrRow).eval();
  RHSD=RHSD.segment(0,nrRow).eval();
  COLD ret=RobustInversion<T,LABEL>::eval(LHSD)*RHSD;
  //std::cout << (LHSD*ret-RHSD) << std::endl;
  return ret;
}
//io
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::read(std::vector<VEC>& v,std::istream& is)
{
  sizeType sz;
  readBinaryData(sz,is);
  v.resize(sz);
  for(sizeType i=0; i<sz; i++)
    read(v[i],is);
  return is.good();
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::read(std::vector<MAT>& m,std::istream& is)
{
  sizeType sz;
  readBinaryData(sz,is);
  m.resize(sz);
  for(sizeType i=0; i<sz; i++)
    read(m[i],is);
  return is.good();
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::write(const std::vector<VEC>& v,std::ostream& os)
{
  sizeType sz=v.size();
  writeBinaryData(sz,os);
  for(sizeType i=0; i<sz; i++)
    write(v[i],os);
  return os.good();
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::write(const std::vector<MAT>& m,std::ostream& os)
{
  sizeType sz=m.size();
  writeBinaryData(sz,os);
  for(sizeType i=0; i<sz; i++)
    write(m[i],os);
  return os.good();
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::read(VEC& v,std::istream& is)
{
  MAT m;
  bool ret=read(m,is);
  v=m;
  return ret;
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::read(MAT& m,std::istream& is)
{
  sizeType rr,cc;
  readBinaryData(rr,is);
  readBinaryData(cc,is);
  m.resize(rr,cc);
  for(sizeType r=0; r<m.rows(); r++)
    for(sizeType c=0; c<m.cols(); c++) {
      scalarD val;
      readBinaryData(val,is);
      m(r,c)=val;
    }
  return is.good();
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::write(const VEC& v,std::ostream& os)
{
  return write((MAT)v,os);
}
template <typename T,char LABEL>
bool SOSPolynomial<T,LABEL>::write(const MAT& m,std::ostream& os)
{
  writeBinaryData((sizeType)m.rows(),os);
  writeBinaryData((sizeType)m.cols(),os);
  for(sizeType r=0; r<m.rows(); r++)
    for(sizeType c=0; c<m.cols(); c++) {
      scalarD val=std::to_double(m(r,c));
      writeBinaryData(val,os);
    }
  return os.good();
}
template <typename T,char LABEL>
SOSPolynomial<T,LABEL>::operator T() const {
  if(_terms.empty())
    return ScalarOfT<T>::convert(0);
  else {
    ASSERT(_terms.size()==1 && _terms[0]._id.empty())
    return _terms[0]._coef;
  }
}
//string
template <typename T,char LABEL>
std::string SOSPolynomial<T,LABEL>::toString(const std::unordered_map<sizeType,std::string>* varNames) const {
  if(_terms.empty())
    return std::string();
  std::string ret=_terms[0].toString(varNames);
  for(sizeType i=1; i<(sizeType)_terms.size(); i++)
    ret+="+"+_terms[i].toString(varNames);
  return ret;
}
template <typename T,char LABEL>
std::string SOSPolynomial<T,LABEL>::formattedString(const std::unordered_map<sizeType,std::string>* varNames) const {
  if(_terms.empty())
    return std::string();
  std::string ret=_terms[0].formattedString(varNames);
  for(sizeType i=1; i<(sizeType)_terms.size(); i++)
    ret+="+"+_terms[i].formattedString(varNames);
  return ret;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::readFormattedString(std::string str) {
  _terms.clear();
  str.erase(std::remove(str.begin(),str.end(),' '),str.end());
  str.erase(std::remove(str.begin(),str.end(),'\n'),str.end());
  int depth=0;
  std::vector<sizeType> cut;
  cut.push_back(-1);
  for(sizeType i=0; i<(sizeType)str.size(); i++) {
    if(str[i]=='(')
      depth++;
    else if(str[i]==')')
      depth--;
    else if(str[i]=='+'&&depth==0)
      cut.push_back(i);
  }
  cut.push_back(str.size());
  for(sizeType i=1; i<(sizeType)cut.size(); i++) {
    sizeType p0=cut[i-1]+1,p1=cut[i];
    SOSTerm<T,LABEL> t;
    t.readFormattedString(str.substr(p0,p1-p0));
    _terms.push_back(t);
  }
  std::sort(_terms.begin(),_terms.end());
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::operator<<(std::string str) {
  str.erase(std::remove(str.begin(),str.end(),' '),str.end());
  str.erase(std::remove(str.begin(),str.end(),'\n'),str.end());
  SOSPolynomialLexer<T,LABEL> parser;
  //parser.testLex(str);
  *this=parser.parse(str);
}
//debug
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::debugIntegrate() const {
  bool bracketNeeded=false;
  COLD L=COLD::Random(nrVar()),U=COLD::Random(nrVar()),LU=concat(L,U);
  T polyIEval=integrate().eval(LU),polyIEvalRef=polyIEval-integrate(L,U);
  std::cout << "Integrate: " << convert(polyIEval,bracketNeeded) << " err: " << convert(polyIEvalRef,bracketNeeded) << std::endl;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::debugGradient() const {
  bool bracketNeeded=false;
  COLD x=COLD::Random(nrVar());
  VEC g=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
  eval(x,&g);
  std::vector<SOSPolynomial> grad=gradient();
  for(sizeType i=0; i<(sizeType)grad.size(); i++) {
    T gradV=grad[i].eval(x),gradVRef=gradV-g[i];
    std::cout << "GradientA: " << convert(gradV,bracketNeeded) << " err: " << convert(gradVRef,bracketNeeded) << std::endl;
  }
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::debugEval() const {
  bool bracketNeeded=false;
  DEFINE_NUMERIC_DELTA
  COLD x=COLD::Random(nrVar());
  COLD dx=COLD::Random(nrVar());
  VEC g=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0)),g2=g;
  MAT h=MAT::Constant(nrVar(),nrVar(),ScalarOfT<T>::convert(0)),h2=h;
  T f=eval(x,&g,&h);
  T f2=eval(x+dx*(typename ScalarOfT<T>::Type)(DELTA),&g2,&h2);
  SMATP hess=hessianSparse();
  {
    T gdxV=ScalarOfT<T>::convert(0),gdxVRef=ScalarOfT<T>::convert(0);
    for(sizeType i=0; i<g.size(); i++) {
      gdxV+=g[i]*ScalarOfT<T>::convert(dx[i]);
      gdxVRef+=g[i]*ScalarOfT<T>::convert(dx[i]);
    }
    gdxVRef-=(f2-f)*ScalarOfT<T>::convert(1/DELTA);
    std::cout << "Gradient: " << convert(gdxV,bracketNeeded) << " err: " << convert(gdxVRef,bracketNeeded) << std::endl;
  }
  {
    VEC hdxV=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
    VEC hdxVRef=VEC::Constant(nrVar(),ScalarOfT<T>::convert(0));
    for(sizeType i=0; i<g.size(); i++)
      for(sizeType j=0; j<g.size(); j++) {
        hdxV[i]+=h(i,j)*ScalarOfT<T>::convert(dx[j]);
        hdxVRef[i]+=h(i,j)*ScalarOfT<T>::convert(dx[j]);
      }
    hdxVRef-=(g2-g)*ScalarOfT<T>::convert(1/DELTA);
    T hdxVNorm=ScalarOfT<T>::convert(0),hdxVRefNorm=ScalarOfT<T>::convert(0);
    for(sizeType i=0; i<hdxV.size(); i++) {
      hdxVNorm+=hdxV[i]*hdxV[i];
      hdxVRefNorm+=hdxVRef[i]*hdxVRef[i];
    }
    std::cout << "Hessian: " << convert(hdxVNorm,bracketNeeded) << " err: " << convert(hdxVRefNorm,bracketNeeded) << std::endl;
  }
  {
    MAT hdx=hessianGradDir(hess,x,dx).toDense();
    MAT hdxVRef=hdx-(h2-h)*ScalarOfT<T>::convert(1/DELTA);
    T hdxVNorm=ScalarOfT<T>::convert(0),hdxVRefNorm=ScalarOfT<T>::convert(0);
    for(sizeType r=0; r<hdx.rows(); ++r)
      for(sizeType c=0; c<hdx.cols(); ++c)
        hdxVNorm+=hdx(r,c)*hdx(r,c);
    for(sizeType r=0; r<hdxVRef.rows(); ++r)
      for(sizeType c=0; c<hdxVRef.cols(); ++c)
        hdxVRefNorm+=hdxVRef(r,c)*hdxVRef(r,c);
    std::cout << "HessianGradDir: " << convert(hdxVNorm,bracketNeeded) << " err: " << convert(hdxVRefNorm,bracketNeeded) << std::endl;
  }
}
//helper
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::add(const SOSTerm<T,LABEL>& other) {
  typename std::vector<SOSTerm<T,LABEL>>::iterator it=std::lower_bound(_terms.begin(),_terms.end(),other);
  if(it==_terms.end()||*it!=other)
    _terms.insert(it,other);
  else *it+=other;
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::consistencyCheck() const {
#ifndef NO_CONSISTENCY_CHECK
  for(sizeType i=0; i<(sizeType)_terms.size()-1; i++) {
    ASSERT(_terms[i]<_terms[i+1])
  }
  for(sizeType i=0; i<(sizeType)_terms.size(); i++)
    _terms[i].consistencyCheck();
#endif
}
template <typename T,char LABEL>
void SOSPolynomial<T,LABEL>::makeConsistent() {
  //sort
  std::sort(_terms.begin(),_terms.end());
  //compact
  sizeType j=0;
  for(sizeType i=1; i<(sizeType)_terms.size(); i++)
    if(_terms[i]==_terms[j])
      _terms[j]._coef+=_terms[i]._coef;
    else _terms[++j]=_terms[i];
  _terms.resize(j+1);
  consistencyCheck();
}
//force instance
PRJ_BEGIN
#define INSTANCE_POLY(TT,POSTFIX)   \
template class SOSTerm<TT,'a'>; \
template class SOSPolynomial<TT,'a'>;
INSTANCE_POLY(double,D)
#ifdef ALL_TYPES
INSTANCE_POLY(__float128,Q)
INSTANCE_POLY(mpfr::mpreal,M)
#endif
PRJ_END
