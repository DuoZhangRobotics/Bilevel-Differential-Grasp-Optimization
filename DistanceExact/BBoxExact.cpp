#include "BBoxExact.h"
#include "MPQZIO.h"

USE_PRJ_NAMESPACE

BBoxExact::BBoxExact() {}
BBoxExact::BBoxExact(const BBox<scalar>& bb):_minC(bb._minC.cast<T>()),_maxC(bb._maxC.cast<T>()) {}
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
EIGEN_DEVICE_FUNC BBoxExact::T BBoxExact::distToSqr(const PT& pt) const {
  PT dist=PT::Zero();
  for(int i=0; i<3; i++) {
    if (pt[i] < _minC[i])
      dist[i] = pt[i] - _minC[i];
    else if (pt[i] > _maxC[i])
      dist[i] = pt[i] - _maxC[i];
  }
  return dist.squaredNorm();
}
