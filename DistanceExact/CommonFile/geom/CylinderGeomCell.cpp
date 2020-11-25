#include "CylinderGeomCell.h"
#include "../MakeMesh.h"

USE_PRJ_NAMESPACE

//CylinderGeomCell
EIGEN_DEVICE_FUNC CylinderGeomCell::CylinderGeomCell():StaticGeomCell(typeid(CylinderGeomCell).name()) {}
EIGEN_DEVICE_FUNC CylinderGeomCell::CylinderGeomCell(const Mat4& T,sizeType dim,scalar rad,scalar y)
  :StaticGeomCell(T,dim,typeid(CylinderGeomCell).name()),_rad(rad),_y(y)
{
  build(false);
}
bool CylinderGeomCell::read(std::istream& is,IOData* dat)
{
  StaticGeomCell::read(is,dat);
  readBinaryData(_rad,is);
  readBinaryData(_y,is);
  return is.good();
}
bool CylinderGeomCell::write(std::ostream& os,IOData* dat) const
{
  StaticGeomCell::write(os,dat);
  writeBinaryData(_rad,os);
  writeBinaryData(_y,os);
  return os.good();
}
std::shared_ptr<SerializableBase> CylinderGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new CylinderGeomCell(*this));
}
scalar CylinderGeomCell::getRad() const
{
  return _rad;
}
scalar CylinderGeomCell::getY() const
{
  return _y;
}
//helper
CylinderGeomCell::CylinderGeomCell(const Mat4& T,sizeType dim,const std::string& type):StaticGeomCell(T,dim,type) {}
void CylinderGeomCell::getMeshInner(ObjMesh& mesh) const
{
  MakeMesh::makeCylinder3D(mesh,_rad,_y,_res,1);
}
DEVICE_ONLY_FUNC BBox<scalar> CylinderGeomCell::getBBInner() const
{
  Vec3 cor(_rad,_y,_rad);
  if(_dim == 2)cor[2]=0.0f;
  return BBox<scalar>(-cor,cor);
}
DEVICE_ONLY_FUNC bool CylinderGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  scalar len=Vec3(pt[0],0.0f,pt[2]).norm();
  //boundary
  scalar dist=_rad-len;
  if(dist < 0.0f)return false;
  n=Vec3(pt[0],0.0f,pt[2])*dist/std::max<scalar>(len,1E-6f);
  //bottom
  scalar dist2=_y+pt[1];
  if(dist2 < 0.0f)return false;
  if(dist2 < dist) {
    dist=dist2;
    n=-Vec3::Unit(1)*dist2;
  }
  //top
  scalar dist3=_y-pt[1];
  if(dist3 < 0.0f)return false;
  if(dist3 < dist) {
    dist=dist3;
    n=Vec3::Unit(1)*dist3;
  }
  return true;
}
DEVICE_ONLY_FUNC bool CylinderGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  Vec3 ptN(pt[0],0.0f,pt[2]);
  scalar len=ptN.norm();
  if(len < _rad && std::abs(pt[1]) < _y) {
    scalar distLen=_rad-len;
    scalar distY1=_y-pt[1];
    scalar distY2=pt[1]+_y;
    if(distLen < distY1 && distLen < distY2) {
      ptN/=std::max<scalarD>(1E-5f,len);
      if(normal)*normal=ptN;
      ptN*=_rad;
      ptN[1]=pt[1];
    } else if(distY1 < distY2 && distY1 < distLen) {
      ptN[1]=_y;
      if(normal)*normal=Vec3::Unit(1);
    } else {
      ptN[1]=-_y;
      if(normal)*normal=-Vec3::Unit(1);
    }
    n=ptN-pt;
    return true;
  } else if(std::abs(pt[1]) < _y) {
    ptN/=std::max<scalarD>(1E-5f,len);
    if(normal)*normal=ptN;
    ptN*=_rad;
    ptN[1]=pt[1];
  } else if(len < _rad) {
    if(pt[1] >= _y) {
      ptN[1]=_y;
      if(normal)*normal=Vec3::Unit(1);
    } else {
      ptN[1]=-_y;
      if(normal)*normal=-Vec3::Unit(1);
    }
  } else {
    ptN*=_rad/std::max<scalarD>(1E-5f,len);
    if(pt[1] >= _y)
      ptN[1]=_y;
    else ptN[1]=-_y;
    if(normal) {
      *normal=pt-ptN;
      *normal/=std::max<scalarD>(1E-5f,normal->norm());
    }
  }
  n=ptN-pt;
  return false;
}
DEVICE_ONLY_FUNC scalar CylinderGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  //solve ||x0+s*dir|| == _rad
  scalar a=Vec2(dir[0],dir[2]).squaredNorm();
  scalar b=Vec2(dir[0],dir[2]).dot(Vec2(x0[0],x0[2]))*2;
  scalar c=Vec2(x0[0],x0[2]).squaredNorm()-_rad*_rad;

  scalar s0=0,s1=1;
  scalar delta=b*b-4*a*c;
  if(a < ScalarUtil<scalar>::scalar_eps() || delta <= 0) {
    if(c >= 0)
      return 1;
  } else {
    s0=std::max<scalar>((-b-std::sqrt(delta))/(2*a),0);
    s1=std::min<scalar>((-b+std::sqrt(delta))/(2*a),1);
    if(s0 >= s1)
      return 1;
  }

  scalar y0=x0[1]+dir[1]*s0;
  scalar y1=x0[1]+dir[1]*s1;
  if(y0 > _y) {
    if(y1 >= _y)
      return 1;
    return (_y-x0[1])/dir[1];
  } else if(y0 < -_y) {
    if(y1 <= -_y)
      return 1;
    return (-_y-x0[1])/dir[1];
  } else return s0;
}
