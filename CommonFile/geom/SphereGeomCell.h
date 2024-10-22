#ifndef SPHERE_GEOM_CELL_H
#define SPHERE_GEOM_CELL_H

#include "StaticGeom.h"
#include <stack>

PRJ_BEGIN

struct ALIGN_16 SphereGeomCell : public StaticGeomCell {
  friend struct CapsuleGeomCell;
  using Serializable::read;
  using Serializable::write;
  SphereGeomCell();
  SphereGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar depth=0.0f);
  virtual BBox<scalar> getBB(bool ref=false) const;
  virtual bool dist(const Vec3& pt,Vec3& n) const;
  virtual bool closest(const Vec3& pt,Vec3& n,Vec3* normal=NULL) const;
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  scalar getRad() const;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  virtual BBox<scalar> getBBInner() const;
  virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  virtual void generateUVInner(ObjMesh& mesh,scalar scale) const;
  ALIGN_16 scalar _rad;
  ALIGN_16 scalar _depth;
};

PRJ_END

#endif
