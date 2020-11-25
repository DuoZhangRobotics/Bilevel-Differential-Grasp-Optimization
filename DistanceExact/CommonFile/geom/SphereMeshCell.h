#ifndef SPHERE_MESH_CELL_H
#define SPHERE_MESH_CELL_H

#include "StaticGeom.h"
#include "ObjMeshGeomCell.h"

PRJ_BEGIN

struct ALIGN_16 TwoSphereMeshCell : public ObjMeshGeomCell {
  using Serializable::read;
  using Serializable::write;
  EIGEN_DEVICE_FUNC TwoSphereMeshCell();
  EIGEN_DEVICE_FUNC TwoSphereMeshCell(sizeType dim,const Vec3& ctr1,const Vec3& ctr2,scalar rad1,scalar rad2);
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  scalarD rad1() const;
  scalarD rad2() const;
  Vec3 ctr1() const;
  Vec3 ctr2() const;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBBInner() const;
  ALIGN_16 scalar _rad1,_rad2,_lenY;
};
struct ALIGN_16 ThreeSphereMeshCell : public ObjMeshGeomCell {
  using Serializable::read;
  using Serializable::write;
  EIGEN_DEVICE_FUNC ThreeSphereMeshCell();
  EIGEN_DEVICE_FUNC ThreeSphereMeshCell(sizeType dim,const Vec3& ctr1,const Vec3& ctr2,const Vec3 ctr3,scalar rad1,scalar rad2,scalar rad3);
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  scalarD rad1() const;
  scalarD rad2() const;
  scalarD rad3() const;
  Vec3 ctr1() const;
  Vec3 ctr2() const;
  Vec3 ctr3() const;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBBInner() const;
  ALIGN_16 scalar _rad1,_rad2,_rad3;
  ALIGN_16 Vec2 _ctr2,_ctr3;
};

PRJ_END

#endif

