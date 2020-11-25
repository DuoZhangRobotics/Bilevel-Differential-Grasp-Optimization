#ifndef BBOX_EXACT_H
#define BBOX_EXACT_H

#include "CommonFile/CollisionDetection.h"
#include "gmpxx.h"

PRJ_BEGIN

struct BBoxExact : public SerializableBase
{
  typedef mpq_class T;
  typedef Eigen::Matrix<T,3,1> PT;
  BBoxExact();
  BBoxExact(const BBox<scalar>& bb);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  T distToSqr(const PT& pt) const;
  ALIGN_16 PT _minC;
  ALIGN_16 PT _maxC;
};

PRJ_END

#endif
