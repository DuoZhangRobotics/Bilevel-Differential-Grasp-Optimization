#include "CompositeGeomCell.h"
#include "../IO.h"

USE_PRJ_NAMESPACE

//CompositeGeomCell
EIGEN_DEVICE_FUNC CompositeGeomCell::CompositeGeomCell():StaticGeomCell(typeid(CompositeGeomCell).name()) {}
EIGEN_DEVICE_FUNC CompositeGeomCell::CompositeGeomCell(const CompositeGeomCell& other):StaticGeomCell(typeid(CompositeGeomCell).name())
{
  _children.clear();
  for(sizeType i=0; i<(sizeType)other._children.size(); i++)
    _children.push_back(std::dynamic_pointer_cast<StaticGeomCell>(other._children[i]->copy()));
  build(false);
}
EIGEN_DEVICE_FUNC CompositeGeomCell::CompositeGeomCell(const Mat4& T,std::vector<std::shared_ptr<StaticGeomCell> > children)
  :StaticGeomCell(T,children[0]->dim(),typeid(CompositeGeomCell).name()),_children(children)
{
  build(false);
}
bool CompositeGeomCell::read(std::istream& is,IOData* dat)
{
  StaticGeomCell::read(is,dat);
  readBinaryData(_children,is,dat);
  return is.good();
}
bool CompositeGeomCell::write(std::ostream& os,IOData* dat) const
{
  StaticGeomCell::write(os,dat);
  writeBinaryData(_children,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> CompositeGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new CompositeGeomCell(*this));
}
void CompositeGeomCell::setRes(sizeType res)
{
  StaticGeomCell::setRes(res);
  for(sizeType i=0; i<(sizeType)_children.size(); i++)
    _children[i]->setRes(res);
}
std::shared_ptr<StaticGeomCell> CompositeGeomCell::getChild(sizeType i) const
{
  return _children[i];
}
sizeType CompositeGeomCell::nrChildren() const
{
  return (sizeType)_children.size();
}
void CompositeGeomCell::getMeshInner(ObjMesh& mesh) const
{
  ObjMesh m;
  mesh=ObjMesh();
  for(sizeType i=0; i<(sizeType)_children.size(); i++) {
    _children[i]->getMesh(m);
    mesh.addMesh(m,"c"+std::to_string(i));
  }
}
DEVICE_ONLY_FUNC BBox<scalar> CompositeGeomCell::getBBInner() const
{
  BBox<scalar> ret;
  for(sizeType i=0; i<(sizeType)_children.size(); i++)
    ret.setUnion(_children[i]->getBB(false));
  return ret;
}
DEVICE_ONLY_FUNC bool CompositeGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  bool ret=false;
  Vec3 nTmp;
  n.setZero();
  n.segment(0,_dim).setConstant(ScalarUtil<scalar>::scalar_max());
  for(sizeType i=0; i<(sizeType)_children.size(); i++)
    if(_children[i]->dist(pt,nTmp)) {
      if(!ret || nTmp.norm() < n.norm())
        n=nTmp;
      ret=true;
    }
  return ret;
}
DEVICE_ONLY_FUNC bool CompositeGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  bool ret=false,retI;
  Vec3 nTmp,normalTmp;
  n.setZero();
  n.segment(0,_dim).setConstant(ScalarUtil<scalar>::scalar_max());
  for(sizeType i=0; i<(sizeType)_children.size(); i++) {
    retI=_children[i]->closest(pt,nTmp,normal ? &normalTmp : NULL);
    if(!ret) {
      if(retI || nTmp.norm() < n.norm()) {
        n=nTmp;
        if(normal)
          *normal=normalTmp;
      }
    } else if(retI && nTmp.norm() < n.norm()) {
      n=nTmp;
      if(normal)
        *normal=normalTmp;
    }
    ret=ret || retI;
  }
  return ret;
}
DEVICE_ONLY_FUNC scalar CompositeGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  scalar s=ScalarUtil<scalar>::scalar_max(),sc;
  for(sizeType i=0; i<(sizeType)_children.size(); i++) {
    sc=_children[i]->rayQuery(x0,dir);
    if(sc < s)
      s=sc;
  }
  return s;
}
