#ifndef TRIANGLE_H
#define TRIANGLE_H

#include "BBoxExact.h"
#include "CommonFile/CollisionDetection.h"
#include <Utils/Scalar.h>
#include <gmpxx.h>

PRJ_BEGIN

class ALIGN_16 TriangleExact : public SerializableBase
{
public:
  typedef mpq_class T;
  typedef Eigen::Matrix<T,3,1> PT;
  typedef Eigen::Matrix<T,2,1> PT2;
  typedef Eigen::Matrix<T,2,2> MAT2;
  typedef std::vector<PT,Eigen::aligned_allocator<PT>> PTss;
  TriangleExact();
  TriangleExact(const TriangleTpl<scalar>& t);
  TriangleExact(const PT& a,const PT& b,const PT& c);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  const PT& normal() const;
  const PT& v(sizeType i) const;
  PT2 project(const PT& d) const;
  bool intersect(const BBoxExact& bb) const;
  void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT& b,Vec2i& feat) const;
  void writeVTK(const std::string& path,const PTss* pt=NULL,const PTss* cp=NULL) const;
  static void debugPointDist(const std::string& path,sizeType nrIter=100);
private:
  void reset();
  ALIGN_16 PT _v[3],_n[3],_nO[3],_nOSqr,_nt;
  ALIGN_16 MAT2 _invM;
};

PRJ_END

#endif
