#include "SphereMeshCell.h"
#include "../MakeMesh.h"
#include "../CameraModel.h"
using std::max;
using std::min;

USE_PRJ_NAMESPACE

//TwoSphereMeshCell
EIGEN_DEVICE_FUNC TwoSphereMeshCell::TwoSphereMeshCell()
{
  setType(typeid(TwoSphereMeshCell).name());
}
EIGEN_DEVICE_FUNC TwoSphereMeshCell::TwoSphereMeshCell(sizeType dim,const Vec3& ctr1,const Vec3& ctr2,scalar rad1,scalar rad2)
  :ObjMeshGeomCell(Mat4::Identity(),dim,typeid(TwoSphereMeshCell).name()),_rad1(rad1),_rad2(rad2),_lenY((ctr1-ctr2).norm())
{
  setType(typeid(TwoSphereMeshCell).name());
  _dim=dim;
  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=ctr1;
  Vec3 ctr21=(ctr2-ctr1).normalized();
  if(dim==2) {
    T(2,3)=0;
    ctr21[2]=0;
  }
  T.block<3,3>(0,0)=ScalarUtil<scalar>::ScalarQuat::FromTwoVectors(Vec3::UnitY(),ctr21).toRotationMatrix();
  setT(T);

  ObjMesh mesh;
  getMeshInner(mesh);
  _depth=mesh.getBB().getExtent().maxCoeff();
  _vss=mesh.getV();
  _iss=mesh.getI();
  build(true);
}
bool TwoSphereMeshCell::read(std::istream& is,IOData* dat)
{
  ObjMeshGeomCell::read(is,dat);
  readBinaryData(_rad1,is);
  readBinaryData(_rad2,is);
  readBinaryData(_lenY,is);
  return is.good();
}
bool TwoSphereMeshCell::write(std::ostream& os,IOData* dat) const
{
  ObjMeshGeomCell::write(os,dat);
  writeBinaryData(_rad1,os);
  writeBinaryData(_rad2,os);
  writeBinaryData(_lenY,os);
  return os.good();
}
std::shared_ptr<SerializableBase> TwoSphereMeshCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new TwoSphereMeshCell(*this));
}
scalarD TwoSphereMeshCell::rad1() const
{
  return _rad1;
}
scalarD TwoSphereMeshCell::rad2() const
{
  return _rad2;
}
Vec3 TwoSphereMeshCell::ctr1() const
{
  return Vec3(0,0,0);
}
Vec3 TwoSphereMeshCell::ctr2() const
{
  return Vec3(0,_lenY,0);
}
//helper
void TwoSphereMeshCell::getMeshInner(ObjMesh& mesh) const
{
  if(_dim == 3)
    MakeMesh::makeTwoSphereMesh3D(mesh,_rad1,_rad2,_lenY,_res);
  else MakeMesh::makeTwoSphereMesh2D(mesh,_rad1,_rad2,_lenY,_res);
}
DEVICE_ONLY_FUNC BBox<scalar> TwoSphereMeshCell::getBBInner() const
{
  Vec3 R=Vec3::Constant(max(_rad1,_rad2));
  return BBox<scalar>(-R,R+Vec3(0,_lenY,0));
}

//ThreeSphereMeshCell
EIGEN_DEVICE_FUNC ThreeSphereMeshCell::ThreeSphereMeshCell()
{
  setType(typeid(ThreeSphereMeshCell).name());
}
EIGEN_DEVICE_FUNC ThreeSphereMeshCell::ThreeSphereMeshCell(sizeType dim,const Vec3& ctr1,const Vec3& ctr2,const Vec3 ctr3,scalar rad1,scalar rad2,scalar rad3)
  :ObjMeshGeomCell(Mat4::Identity(),dim,typeid(ThreeSphereMeshCell).name()),_rad1(rad1),_rad2(rad2),_rad3(rad3)
{
  setType(typeid(ThreeSphereMeshCell).name());
  _dim=dim;
  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=ctr1;
  T.block<3,3>(0,0)=ScalarUtil<scalar>::ScalarQuat::FromTwoVectors(Vec3::UnitZ(),(ctr2-ctr1).cross(ctr3-ctr1).normalized()).toRotationMatrix();
  setT(T);
  _ctr2=transformHomo<scalar>(Mat4(T.inverse()),ctr2).segment<2>(0);
  _ctr3=transformHomo<scalar>(Mat4(T.inverse()),ctr3).segment<2>(0);

  ObjMesh mesh;
  getMeshInner(mesh);
  _depth=mesh.getBB().getExtent().maxCoeff();
  _vss=mesh.getV();
  _iss=mesh.getI();
  build(true);
}
bool ThreeSphereMeshCell::read(std::istream& is,IOData* dat)
{
  ObjMeshGeomCell::read(is,dat);
  readBinaryData(_rad1,is);
  readBinaryData(_rad2,is);
  readBinaryData(_rad3,is);
  readBinaryData(_ctr2,is);
  readBinaryData(_ctr3,is);
  return is.good();
}
bool ThreeSphereMeshCell::write(std::ostream& os,IOData* dat) const
{
  ObjMeshGeomCell::write(os,dat);
  writeBinaryData(_rad1,os);
  writeBinaryData(_rad2,os);
  writeBinaryData(_rad3,os);
  writeBinaryData(_ctr2,os);
  writeBinaryData(_ctr3,os);
  return os.good();
}
std::shared_ptr<SerializableBase> ThreeSphereMeshCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new ThreeSphereMeshCell(*this));
}
scalarD ThreeSphereMeshCell::rad1() const
{
  return _rad1;
}
scalarD ThreeSphereMeshCell::rad2() const
{
  return _rad2;
}
scalarD ThreeSphereMeshCell::rad3() const
{
  return _rad3;
}
Vec3 ThreeSphereMeshCell::ctr1() const
{
  return Vec3(0,0,0);
}
Vec3 ThreeSphereMeshCell::ctr2() const
{
  return Vec3(_ctr2[0],_ctr2[1],0);
}
Vec3 ThreeSphereMeshCell::ctr3() const
{
  return Vec3(_ctr3[0],_ctr3[1],0);
}
//helper
void ThreeSphereMeshCell::getMeshInner(ObjMesh& mesh) const
{
  if(_dim == 3)
    MakeMesh::makeThreeSphereMesh3D(mesh,_rad1,_rad2,_rad3,Vec2::Zero(),_ctr2,_ctr3,_res);
  else MakeMesh::makeThreeSphereMesh3D(mesh,_rad1,_rad2,_rad3,Vec2::Zero(),_ctr2,_ctr3,_res);
}
DEVICE_ONLY_FUNC BBox<scalar> ThreeSphereMeshCell::getBBInner() const
{
  BBox<scalar> bb;
  bb.setPoints(Vec3(0,0,0),Vec3(_ctr2[0],_ctr2[1],0),Vec3(_ctr3[0],_ctr3[1],0));
  bb.enlarged(std::max(_rad1,std::max(_rad2,_rad3)));
  return bb;
}
