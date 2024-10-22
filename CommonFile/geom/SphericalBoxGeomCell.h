#ifndef SPHERICAL_BOX_CELL_H
#define SPHERICAL_BOX_CELL_H

#include "BoxGeomCell.h"

PRJ_BEGIN

struct ALIGN_16 SphericalBoxGeomCell : public BoxGeomCell {
  using Serializable::read;
  using Serializable::write;
  SphericalBoxGeomCell();
  SphericalBoxGeomCell(const Mat4& T,sizeType dim,const Vec4& ext);
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  scalar getRad() const;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual bool distInner(const Vec3& pt,Vec3& n) const;
  virtual bool closestInner(const Vec3& pt,Vec3& n,Vec3* normal=NULL) const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  ALIGN_16 scalar _rad;
};

PRJ_END

#endif
