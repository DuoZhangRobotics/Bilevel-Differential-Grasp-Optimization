#ifndef BBOX_EXACT_H
#define BBOX_EXACT_H

#include <CommonFile/CollisionDetection.h>
#include <Utils/Scalar.h>
#include <gmpxx.h>

PRJ_BEGIN

struct BBoxExact : public SerializableBase
{
  typedef mpq_class T;
  typedef Eigen::Matrix<T,3,1> PT;
  typedef Eigen::Matrix<T,2,1> PT2;
  BBoxExact();
  BBoxExact(const BBox<scalar>& bb);
  BBoxExact(const PT& a,const PT& b);
  BBoxExact(const PT& a,const PT& b,const PT& c);
  void setUnion(const PT& other);
  PT2 project(const PT& d) const;
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
