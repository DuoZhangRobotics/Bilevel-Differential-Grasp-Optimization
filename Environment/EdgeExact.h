#ifndef EDGE_EXACT_H
#define EDGE_EXACT_H

#include "BBoxExact.h"
#include "CommonFile/CollisionDetection.h"
#include <Utils/Scalar.h>
#include <gmpxx.h>

PRJ_BEGIN

class ALIGN_16 EdgeExact : public SerializableBase
{
public:
  typedef mpq_class T;
  typedef Eigen::Matrix<T,3,1> PT;
  typedef Eigen::Matrix<T,3,3> MAT3;
  typedef std::vector<PT,Eigen::aligned_allocator<PT>> PTss;
  EdgeExact();
  EdgeExact(const LineSegTpl<scalar>& t);
  EdgeExact(const PT& a,const PT& b);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  std::pair<sizeType,PT> moveFromVertex(const PT& d0,const PT& D) const;
  T dirGrad(sizeType vid,const PT& D) const;
  Vec2i _vid,_tNId;
private:
  void reset();
  ALIGN_16 PT _v[2],_d;
  ALIGN_16 MAT3 _dTdN;
  ALIGN_16 T _lenSqr;
};

PRJ_END

#endif
