#ifndef OBJ_MESH_GEOM_CELL_H
#define OBJ_MESH_GEOM_CELL_H

#include "StaticGeom.h"
#include "../GridBasic.h"

PRJ_BEGIN

template <typename T,typename TI,typename TG>struct Grid;
typedef Grid<scalar,scalar,std::vector<scalar,Eigen::aligned_allocator<scalar> > > ScalarField;

struct ALIGN_16 ObjMeshGeomCell : public StaticGeomCell {
  using Serializable::read;
  using Serializable::write;
  EIGEN_DEVICE_FUNC ObjMeshGeomCell();
  EIGEN_DEVICE_FUNC ObjMeshGeomCell(const Mat4& T,const ObjMesh& mesh,scalar depth,bool buildBVHAsWell=false,bool bottomUp=true);
  bool read(std::istream& is,IOData* dat);
  bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  scalar depth() const;
  scalar& depth();
  void calcMinDist2D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist,Vec2i& feat) const;
  void calcMinDist3D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist,Vec2i& feat) const;
  void buildLevelSet(scalar cellSz,scalar off,scalar eps=0.0f);
  scalar distLevelSet(const Vec3& pos) const;
  const ScalarField& getLevelSet() const;
  virtual void setRes(sizeType res);
protected:
  ObjMeshGeomCell(const std::string& name);
  ObjMeshGeomCell(const Mat4& T,sizeType dim,const std::string& name);
  virtual void getMeshInner(ObjMesh& mesh) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBBInner() const;
  DEVICE_ONLY_FUNC virtual bool distInner(const Vec3& pt,Vec3& n) const;
  DEVICE_ONLY_FUNC virtual bool closestInner(const Vec3& pt,Vec3& n,Vec3* normal=NULL) const;
  DEVICE_ONLY_FUNC virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  ALIGN_16 ScalarField _grid;
  ALIGN_16 scalar _depth;
};

PRJ_END

#endif
