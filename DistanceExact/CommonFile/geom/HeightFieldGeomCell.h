#ifndef HEIGHT_FIELD_GEOM_CELL_H
#define HEIGHT_FIELD_GEOM_CELL_H

#include "ObjMeshGeomCell.h"
#include <stack>

PRJ_BEGIN

struct ALIGN_16 HeightFieldGeomCell : public ObjMeshGeomCell {
  using Serializable::read;
  using Serializable::write;
  EIGEN_DEVICE_FUNC HeightFieldGeomCell();
  EIGEN_DEVICE_FUNC HeightFieldGeomCell(const Mat4& T,const ScalarField& h);
  EIGEN_DEVICE_FUNC HeightFieldGeomCell(sizeType dim,sizeType dimH,scalar h0,scalar hr,scalar sz,scalar cellSz);
  bool read(std::istream& is,IOData* dat);
  bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
protected:
  virtual Vec3i getStride(bool cell) const;
  virtual void getMeshInner(ObjMesh& mesh) const;
  DEVICE_ONLY_FUNC virtual bool distInner(const Vec3& pt,Vec3& n) const;
  DEVICE_ONLY_FUNC virtual bool closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBBInner() const;
  virtual void build(bool buildBVHAsWell=true);
  virtual void debugWrite() const;
  Mat4 genT(sizeType dim,sizeType dimH) const;
  void pointDistQuery(Vec3 pt,Vec3& cp,Vec3& n,scalar& dist,scalar& minDist) const;
  void addRing(std::stack<sizeType>& ss,const Vec3i& id,sizeType ring) const;
  ALIGN_16 BBox<scalar> _bb;
  ALIGN_16 ScalarField _h;
};

PRJ_END

#endif
