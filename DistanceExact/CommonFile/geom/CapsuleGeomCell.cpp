#include "CapsuleGeomCell.h"
#include "SphereGeomCell.h"
#include "../MakeMesh.h"

USE_PRJ_NAMESPACE

//CapsuleGeomCell
EIGEN_DEVICE_FUNC CapsuleGeomCell::CapsuleGeomCell():CylinderGeomCell(Mat4::Identity(),3,typeid(CapsuleGeomCell).name()) {}
EIGEN_DEVICE_FUNC CapsuleGeomCell::CapsuleGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar y)
  :CylinderGeomCell(T,dim,typeid(CapsuleGeomCell).name())
{
  _rad=rad;
  _y=y;
  build(false);
}
std::shared_ptr<SerializableBase> CapsuleGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new CapsuleGeomCell(*this));
}
//helper
void CapsuleGeomCell::getMeshInner(ObjMesh& mesh) const
{
  if(_dim == 3)
    MakeMesh::makeCapsule3D(mesh,_rad,_y,_res);
  else MakeMesh::makeCapsule2D(mesh,_rad,_y,_res);
}
DEVICE_ONLY_FUNC BBox<scalar> CapsuleGeomCell::getBBInner() const
{
  Vec3 cor(_rad,_y+_rad,_rad);
  if(_dim == 2)cor[2]=0.0f;
  return BBox<scalar>(-cor,cor);
}
DEVICE_ONLY_FUNC bool CapsuleGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  n=pt;
  if(!(n[1] >= _y || n[1] <= -_y))
    n[1]=0.0f;
  else if(n[1] >= _y)
    n-=Vec3::Unit(1)*_y;
  else n+=Vec3::Unit(1)*_y;

  scalar norm=std::max<scalar>(ScalarUtil<scalar>::scalar_eps(),n.norm());
  if(norm > _rad)
    return false;
  n*=((_rad-norm)/norm);
  return true;
}
DEVICE_ONLY_FUNC bool CapsuleGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  if(pt[1] < -_y) {
    n=pt+Vec3::UnitY()*_y;
  } else if(pt[1] > _y) {
    n=pt-Vec3::UnitY()*_y;
  } else {
    n=pt;
    n[1]=0;
  }
  scalar len=n.norm();
  Vec3 nn=n/std::max<scalar>(len,ScalarUtil<scalar>::scalar_eps());
  if(normal)
    *normal=nn;
  n=nn*(_rad-len);
  return len < _rad;
}
DEVICE_ONLY_FUNC scalar CapsuleGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  SphereGeomCell c(Mat4::Identity(),_dim,_rad);
  scalar s0=CylinderGeomCell::rayQueryInner(x0,dir);
  scalar s1=c.rayQueryInner(x0-Vec3::Unit(1)*_y,dir);
  scalar s2=c.rayQueryInner(x0+Vec3::Unit(1)*_y,dir);
  return std::min(s0,std::min(s1,s2));
}
