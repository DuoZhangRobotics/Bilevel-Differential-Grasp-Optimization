#ifndef OBJ_MESH_GEOM_CELL_EXACT_H
#define OBJ_MESH_GEOM_CELL_EXACT_H

#include "CommonFile/geom/ObjMeshGeomCell.h"
#include "TriangleExact.h"
#include "BBoxExact.h"

PRJ_BEGIN

struct ALIGN_16 ObjMeshGeomCellExact : public SerializableBase
{
  typedef mpq_class T;
  typedef Eigen::Matrix<T,3,1> PT;
  typedef Eigen::Matrix<T,3,3> MAT3;
  EIGEN_DEVICE_FUNC ObjMeshGeomCellExact();
  EIGEN_DEVICE_FUNC ObjMeshGeomCellExact(const ObjMeshGeomCell& exact);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  DEVICE_ONLY_FUNC const BBoxExact& getBB() const;
  DEVICE_ONLY_FUNC virtual bool closest(const PT& pt,PT& n,PT& normal,MAT3& hessian,Vec2i& feat) const;
  DEVICE_ONLY_FUNC virtual scalar closest(const Vec3& pt,Vec3& n,Vec3& normal,Mat3& hessian,Vec2i& feat) const;
  void writePointDistVTK(const std::string& path,sizeType res=10) const;
  void debugPointDist(sizeType nrIter=100) const;
  sizeType nrV() const;
  sizeType nrI() const;
protected:
  ALIGN_16 std::vector<TriangleExact> _tss;
  ALIGN_16 std::vector<PT,Eigen::aligned_allocator<PT>> _vss;
  ALIGN_16 std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>> _iss;
  ALIGN_16 std::vector<Node<sizeType,BBoxExact>> _bvh;
};

PRJ_END

#endif
