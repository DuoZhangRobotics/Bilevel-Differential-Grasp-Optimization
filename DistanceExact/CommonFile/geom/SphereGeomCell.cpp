#include "SphereGeomCell.h"
#include "../MakeMesh.h"

USE_PRJ_NAMESPACE

//SphereGeomCell
EIGEN_DEVICE_FUNC SphereGeomCell::SphereGeomCell():StaticGeomCell(typeid(SphereGeomCell).name()) {}
EIGEN_DEVICE_FUNC SphereGeomCell::SphereGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar depth)
  :StaticGeomCell(T,dim,typeid(SphereGeomCell).name()),_rad(rad),_depth(depth)
{
  build(false);
}
DEVICE_ONLY_FUNC BBox<scalar> SphereGeomCell::getBB(bool ref) const
{
  if(ref)
    return getBBInner();
  else {
    BBox<scalar> bb=getBBInner();
    Vec3 ctr=_T.block<3,1>(0,3);
    return BBox<scalar>(bb._minC+ctr,bb._maxC+ctr);
  }
}
DEVICE_ONLY_FUNC bool SphereGeomCell::dist(const Vec3& pt,Vec3& n) const
{
  Vec3 dir=pt-_T.block<3,1>(0,3);
  if((dir.squaredNorm() <= _rad*_rad)) {
    closest(pt,n);
    return true;
  } else return false;
}
DEVICE_ONLY_FUNC bool SphereGeomCell::closest(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  return Sphere<scalar>(_T.block<3,1>(0,3),_rad).closest(pt,n,normal);
}
bool SphereGeomCell::read(std::istream& is,IOData* dat)
{
  StaticGeomCell::read(is,dat);
  readBinaryData(_rad,is);
  readBinaryData(_depth,is);
  return is.good();
}
bool SphereGeomCell::write(std::ostream& os,IOData* dat) const
{
  StaticGeomCell::write(os,dat);
  writeBinaryData(_rad,os);
  writeBinaryData(_depth,os);
  return os.good();
}
std::shared_ptr<SerializableBase> SphereGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new SphereGeomCell(*this));
}
scalar SphereGeomCell::getRad() const
{
  return _rad;
}
//helper
void SphereGeomCell::getMeshInner(ObjMesh& mesh) const
{
  if(_dim == 2)MakeMesh::makeSphere2D(mesh,_rad,_res);
  else MakeMesh::makeSphere3D(mesh,_rad,_res);
}
DEVICE_ONLY_FUNC BBox<scalar> SphereGeomCell::getBBInner() const
{
  Vec3 rad=Vec3::Zero();
  rad.block(0,0,_dim,1).setConstant(_rad+_depth);
  return BBox<scalar>(-rad,rad);
}
DEVICE_ONLY_FUNC scalar SphereGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  //solve ||x0+s*dir|| == _rad
  scalar a=dir.squaredNorm();
  scalar b=dir.dot(x0)*2;
  scalar c=x0.squaredNorm()-_rad*_rad;
  scalar delta=b*b-4*a*c;
  if(delta <= 0)
    return 1;
  else {
    scalar s0=(-b-std::sqrt(delta))/(2*a);
    scalar s1=(-b+std::sqrt(delta))/(2*a);
    if(s0 > 0 && s0 < 1)
      return s0;
    else if(s1 > 0 && s1 < 1)
      return s1;
    else return 1;
  }
}
void SphereGeomCell::generateUVInner(ObjMesh& mesh,scalar scale) const
{
  ASSERT(_dim == 3)
  std::vector<Vec2,Eigen::aligned_allocator<Vec2> >& tss=mesh._tex._uv;
  tss.resize(mesh.getV().size());
  //generate texture coordinates
  for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++) {
    const Vec3& v=mesh.getV()[i];
    tss[i][1]=atan2(v[2],v.segment<2>(0).norm())*scale;
    tss[i][0]=atan2(v[1],v[0])*scale;
  }
}
