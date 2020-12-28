#ifndef CAPSULE_GEOM_CELL_H
#define CAPSULE_GEOM_CELL_H

#include "CylinderGeomCell.h"

PRJ_BEGIN

struct ALIGN_16 CapsuleGeomCell : public CylinderGeomCell {
  using Serializable::read;
  using Serializable::write;
  CapsuleGeomCell();
  CapsuleGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar y);
  std::shared_ptr<SerializableBase> copy() const;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual bool distInner(const Vec3& pt,Vec3& n) const;
  virtual bool closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
};

PRJ_END

#endif
