#include "SphericalBoxGeomCell.h"
#include "../MakeMesh.h"

USE_PRJ_NAMESPACE

//SphericalBoxGeomCell
EIGEN_DEVICE_FUNC SphericalBoxGeomCell::SphericalBoxGeomCell()
{
  setType(typeid(SphericalBoxGeomCell).name());
}
EIGEN_DEVICE_FUNC SphericalBoxGeomCell::SphericalBoxGeomCell(const Mat4& T,sizeType dim,const Vec4& ext)
  :BoxGeomCell(T,dim,ext.segment<3>(0))
{
  setType(typeid(SphericalBoxGeomCell).name());
  _rad=ext[3];
  build(false);
}
bool SphericalBoxGeomCell::read(std::istream& is,IOData* dat)
{
  BoxGeomCell::read(is,dat);
  readBinaryData(_rad,is);
  return is.good();
}
bool SphericalBoxGeomCell::write(std::ostream& os,IOData* dat) const
{
  BoxGeomCell::write(os,dat);
  writeBinaryData(_rad,os);
  return os.good();
}
std::shared_ptr<SerializableBase> SphericalBoxGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new SphericalBoxGeomCell(*this));
}
scalar SphericalBoxGeomCell::getRad() const
{
  return _rad;
}
//helper
void SphericalBoxGeomCell::getMeshInner(ObjMesh& mesh) const
{
  if(_dim==2)
    MakeMesh::makeCapsule2D(mesh,_rad,_ext[0],_ext[1],_res);
  else MakeMesh::makeCapsule3D(mesh,_rad,_ext[0],_ext[1],_ext[2],_res);
}
DEVICE_ONLY_FUNC BBox<scalar> SphericalBoxGeomCell::getBBInner() const
{
  return BoxGeomCell::getBBInner().enlarge(_rad,_dim);
}
DEVICE_ONLY_FUNC bool SphericalBoxGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  return closestInner(pt,n,NULL);
}
DEVICE_ONLY_FUNC bool SphericalBoxGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  Vec3 normalTmp;
  bool inside=BoxGeomCell::closestInner(pt,n,&normalTmp);
  if(inside)
    n+=normalTmp*_rad;
  else if(n.squaredNorm() > _rad*_rad)
    n-=normalTmp*_rad;
  else {
    n+=normalTmp*_rad;
    inside=true;
  }
  if(normal)
    *normal=normalTmp;
  return inside;
}
DEVICE_ONLY_FUNC scalar SphericalBoxGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  FUNCTION_NOT_IMPLEMENTED
  return 0;
}
