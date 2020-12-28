#include "BBoxExact.h"
#include "MPQZIO.h"

USE_PRJ_NAMESPACE

BBoxExact::BBoxExact() {}
BBoxExact::BBoxExact(const PT& a,const PT& b):_minC(a.cwiseMin(b)),_maxC(a.cwiseMax(b)) {}
BBoxExact::BBoxExact(const PT& a,const PT& b,const PT& c):_minC(a.cwiseMin(b).cwiseMin(c)),_maxC(a.cwiseMax(b).cwiseMax(c)) {}
BBoxExact::BBoxExact(const BBox<scalar>& bb):_minC(bb._minC.cast<double>().cast<T>()),_maxC(bb._maxC.cast<double>().cast<T>()) {}
void BBoxExact::setUnion(const PT& other)
{
  for(sizeType d=0; d<3; d++) {
    _minC[d]=std::min<T>(_minC[d],other[d]);
    _maxC[d]=std::max<T>(_maxC[d],other[d]);
  }
}
BBoxExact::PT2 BBoxExact::project(const PT& a) const
{
  PT ctr=(_minC+_maxC)*(T)0.5f;
  T ctrD=a.dot(ctr);
  T delta=0.0f;
  ctr=_maxC-ctr;
  for(sizeType i=0; i<3; i++) {
    T coef=ctr[i]*a[i];
    delta+=coef<0?-coef:coef;
  }
  return PT2(ctrD-delta,ctrD+delta);
}
bool BBoxExact::read(std::istream& is,IOData*)
{
  readBinaryData(_minC,is);
  readBinaryData(_maxC,is);
  return is.good();
}
bool BBoxExact::write(std::ostream& os,IOData*) const
{
  writeBinaryData(_minC,os);
  writeBinaryData(_maxC,os);
  return os.good();
}
std::shared_ptr<SerializableBase> BBoxExact::copy() const
{
  return std::shared_ptr<SerializableBase>(new BBoxExact);
}
std::string BBoxExact::type() const
{
  return typeid(BBoxExact).name();
}
BBoxExact::T BBoxExact::distToSqr(const PT& pt) const {
  PT dist=PT::Zero();
  for(int i=0; i<3; i++) {
    if (pt[i] < _minC[i])
      dist[i] = pt[i] - _minC[i];
    else if (pt[i] > _maxC[i])
      dist[i] = pt[i] - _maxC[i];
  }
  return dist.squaredNorm();
}
