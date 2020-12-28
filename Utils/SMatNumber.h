#ifndef SMAT_NUMBER_H
#define SMAT_NUMBER_H

#include "SparseUtils.h"

PRJ_BEGIN

template <typename T>
struct SMatNumber : public SerializableBase
{
public:
  typedef Eigen::SparseMatrix<T,0,sizeType> SMat;
  SMatNumber() {
    _SMat.resize(1,1);
    _SMat.coeffRef(0,0)=0;
  }
  SMatNumber(int other) {
    *this=other;
  }
  SMatNumber(sizeType other) {
    *this=other;
  }
#ifdef USE_QUADMATH
  SMatNumber(double other) {
    *this=other;
  }
#endif
  SMatNumber(T other) {
    *this=other;
  }
  SMatNumber(const SMatNumber& other) {
    *this=other;
  }
  SMatNumber(const SMat& sMat) {
    *this=sMat;
  }
  virtual bool read(std::istream& is,IOData* dat) {
    readBinaryData(_SMat,is);
    return is.good();
  }
  virtual bool write(std::ostream& os,IOData* dat) const {
    writeBinaryData(_SMat,os);
    return os.good();
  }
  virtual std::shared_ptr<SerializableBase> copy() const {
    return std::shared_ptr<SerializableBase>(new SMatNumber);
  }
  virtual std::string type() const {
    return typeid(SMatNumber<T>).name();
  }
  SMatNumber operator+(const SMatNumber& other) const {
    if(isScalar() && operator T()==0)
      return other;
    else if(other.isScalar() && other.operator T()==0)
      return *this;
    else return SMatNumber(_SMat+other._SMat);
  }
  SMatNumber& operator+=(const SMatNumber& other) {
    *this=operator+(other);
    return *this;
  }
  SMatNumber operator-(const SMatNumber& other) const {
    return SMatNumber(_SMat-other._SMat);
  }
  SMatNumber& operator-=(const SMatNumber& other) {
    *this=operator-(other);
    return *this;
  }
  SMatNumber operator*(const SMatNumber& other) const {
    if(isScalar())
      return SMatNumber(other._SMat*operator T());
    else if(other.isScalar())
      return SMatNumber(_SMat*other.operator T());
    else return SMatNumber(_SMat*other._SMat);
  }
  SMatNumber& operator*=(const SMatNumber& other) {
    *this=operator*(other);
    return *this;
  }
  SMatNumber operator*(T other) const {
    return SMatNumber(_SMat*other);
  }
  SMatNumber& operator*=(T other) {
    *this=operator*(other);
    return *this;
  }
  SMatNumber operator/(const SMatNumber& other) const {
    return SMatNumber(_SMat*(1/other.operator T()));
  }
  SMatNumber& operator/=(const SMatNumber& other) {
    *this=operator/(other);
    return *this;
  }
  SMatNumber& operator=(int other) {
    _SMat.resize(1,1);
    _SMat.coeffRef(0,0)=other;
    return *this;
  }
  SMatNumber& operator=(sizeType other) {
    _SMat.resize(1,1);
    _SMat.coeffRef(0,0)=other;
    return *this;
  }
#ifdef USE_QUADMATH
  SMatNumber& operator=(double other) {
    _SMat.resize(1,1);
    _SMat.coeffRef(0,0)=other;
    return *this;
  }
#endif
  SMatNumber& operator=(T other) {
    _SMat.resize(1,1);
    _SMat.coeffRef(0,0)=other;
    return *this;
  }
  SMatNumber& operator=(const SMatNumber& other) {
    _SMat=other._SMat;
    return *this;
  }
  SMatNumber& operator=(const SMat& other) {
    _SMat=other;
    return *this;
  }
  bool operator<(const SMatNumber& other) const {
    //rows
    if(_SMat.rows()<other._SMat.rows())
      return true;
    else if(_SMat.rows()>other._SMat.rows())
      return false;
    //cols
    if(_SMat.cols()<other._SMat.cols())
      return true;
    else if(_SMat.cols()>other._SMat.cols())
      return false;
    //number compare
    for(sizeType c=0; c<_SMat.cols(); c++) {
      typename SMat::InnerIterator it(_SMat,c),it2(other._SMat,c);
      while(it && it2) {
        //rows
        if(it.row()<it2.row())
          return true;
        else if(it.row()>it2.row())
          return false;
        //cols
        if(it.col()<it2.col())
          return true;
        else if(it.col()>it2.col())
          return false;
        //value
        if(it.value()<it2.value())
          return true;
        else if(it.value()>it2.value())
          return false;
      }
      if(it2)
        return true;
      else if(it)
        return false;
    }
    return false;
  }
  bool operator>(const SMatNumber& other) const {
    return other.operator<(*this);
  }
#ifdef USE_QUADMATH
  bool operator==(double other) const {
    return operator==(SMatNumber(other));
  }
  bool operator!=(double other) const {
    return operator!=(SMatNumber(other));
  }
#endif
  bool operator==(T other) const {
    return operator==(SMatNumber(other));
  }
  bool operator!=(T other) const {
    return operator!=(SMatNumber(other));
  }
  bool operator==(const SMatNumber& other) const {
    return !(operator<(other)) && !(other.operator<(*this));
  }
  bool operator!=(const SMatNumber& other) const {
    return !operator==(other);
  }
  const SMat& getSMat() const {
    return _SMat;
  }
  operator T() const {
    ASSERT(isScalar())
    return _SMat.coeff(0,0);
  }
  bool isScalar() const {
    return _SMat.size()==1;
  }
private:
  SMat _SMat;
};

PRJ_END

#endif
