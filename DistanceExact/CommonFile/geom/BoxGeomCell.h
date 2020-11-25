#ifndef BOX_GEOM_CELL_H
#define BOX_GEOM_CELL_H

#include "StaticGeom.h"

PRJ_BEGIN

struct ALIGN_16 BoxGeomCell : public StaticGeomCell {
  using Serializable::read;
  using Serializable::write;
  EIGEN_DEVICE_FUNC BoxGeomCell();
  EIGEN_DEVICE_FUNC BoxGeomCell(const Mat4& T,sizeType dim,const Vec3& ext,scalar depth=0.0f);
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  const Vec3& getExt() const;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBBInner() const;
  DEVICE_ONLY_FUNC virtual bool distInner(const Vec3& pt,Vec3& n) const;
  DEVICE_ONLY_FUNC virtual bool closestInner(const Vec3& pt,Vec3& n,Vec3* normal=NULL) const;
  DEVICE_ONLY_FUNC virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  virtual void generateUVInner(ObjMesh& mesh,scalar scale) const;
  ALIGN_16 Vec3 _ext;
  ALIGN_16 scalar _depth;
};

PRJ_END

#endif

