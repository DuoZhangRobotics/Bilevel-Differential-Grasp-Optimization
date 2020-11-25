#ifndef CAPSULE_GEOM_CELL_H
#define CAPSULE_GEOM_CELL_H

#include "CylinderGeomCell.h"

PRJ_BEGIN

struct ALIGN_16 CapsuleGeomCell : public CylinderGeomCell {
  using Serializable::read;
  using Serializable::write;
  EIGEN_DEVICE_FUNC CapsuleGeomCell();
  EIGEN_DEVICE_FUNC CapsuleGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar y);
  std::shared_ptr<SerializableBase> copy() const;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBBInner() const;
  DEVICE_ONLY_FUNC virtual bool distInner(const Vec3& pt,Vec3& n) const;
  DEVICE_ONLY_FUNC virtual bool closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const;
  DEVICE_ONLY_FUNC virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
};

PRJ_END

#endif
