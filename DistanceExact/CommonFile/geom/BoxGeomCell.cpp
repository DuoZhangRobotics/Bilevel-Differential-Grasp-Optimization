#include "BoxGeomCell.h"
#include "../MakeMesh.h"

USE_PRJ_NAMESPACE

//BoxGeomCell
EIGEN_DEVICE_FUNC BoxGeomCell::BoxGeomCell():StaticGeomCell(typeid(BoxGeomCell).name()) {}
EIGEN_DEVICE_FUNC BoxGeomCell::BoxGeomCell(const Mat4& T,sizeType dim,const Vec3& ext,scalar depth)
  :StaticGeomCell(T,dim,typeid(BoxGeomCell).name()),_ext(ext),_depth(depth)
{
  build(false);
}
bool BoxGeomCell::read(std::istream& is,IOData* dat)
{
  StaticGeomCell::read(is,dat);
  readBinaryData(_ext,is);
  readBinaryData(_depth,is);
  return is.good();
}
bool BoxGeomCell::write(std::ostream& os,IOData* dat) const
{
  StaticGeomCell::write(os,dat);
  writeBinaryData(_ext,os);
  writeBinaryData(_depth,os);
  return os.good();
}
std::shared_ptr<SerializableBase> BoxGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new BoxGeomCell(*this));
}
const Vec3& BoxGeomCell::getExt() const
{
  return _ext;
}
//helper
void BoxGeomCell::getMeshInner(ObjMesh& mesh) const
{
  if(_dim == 2)MakeMesh::makeBox2D(mesh,_ext);
  else MakeMesh::makeBox3D(mesh,_ext);
}
DEVICE_ONLY_FUNC BBox<scalar> BoxGeomCell::getBBInner() const
{
  Vec3 depth=Vec3::Zero();
  depth.block(0,0,_dim,1).setConstant(_depth);
  return BBox<scalar>(-_ext-depth,_ext+depth);
}
DEVICE_ONLY_FUNC bool BoxGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  BBox<scalar> box(-_ext,_ext);
  if(box.contain(pt,_dim)) {
    closestInner(pt,n);
    return true;
  } else return false;
}
DEVICE_ONLY_FUNC bool BoxGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  if(_dim==2) {
    Vec2 n2,normal2;
    return OBBTpl<scalar,2>(_T.block<3,3>(0,0),_T.block<3,1>(0,3),_ext).closestInner(pt.segment<2>(0),n2,normal?&normal2:NULL);
    n=Vec3(n2[0],n2[1],0);
    if(normal)
      *normal=Vec3(normal2[0],normal2[1],0);
  } else {
    return OBBTpl<scalar,3>(_T.block<3,3>(0,0),_T.block<3,1>(0,3),_ext).closestInner(pt,n,normal);
  }
}
DEVICE_ONLY_FUNC scalar BoxGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  scalar s,t;
  BBox<scalar> bb=getBBInner();
  if(bb.intersect(x0,x0+dir,s,t,_dim))
    return s;
  else return 1;
}
void BoxGeomCell::generateUVInner(ObjMesh& mesh,scalar scale) const
{
  ASSERT(_dim == 3)
  mesh.smooth();
  mesh=mesh.cutOpen(M_PI/4);
  std::vector<Vec2,Eigen::aligned_allocator<Vec2> >& tss=mesh._tex._uv;
  tss.resize(mesh.getV().size());
  //generate texture coordinates
  BBox<scalar> bb=mesh.getBB();
  Vec3 ext=bb.getExtent();
  for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++) {
    Vec2& t=tss[i];
    Vec3 n=mesh.getN()[i];
    Vec3 v=mesh.getV()[i]-bb._minC;
    v.array()/=ext.array();
    for(sizeType d=0; d<3; d++)
      if(n[d] < -0.99f || n[d] > 0.99f) {
        t[0]=v[(d+1)%3]*scale;
        t[1]=v[(d+2)%3]*scale;
      }
  }
}
