#include "ObjMeshGeomCell.h"
#include "../ImplicitFunc.h"
#include "../CameraModel.h"
#include "BVHBuilder.h"

USE_PRJ_NAMESPACE

//ObjMeshGeomCell
struct LineCallback {
  LineCallback(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss,
               const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss,
               const LineSeg& l,sizeType dim):_vss(vss),_iss(iss),_l(l),_dim(dim),_s(1) {}
  bool validNode(const Node<sizeType>& node) {
    return node._bb.intersect(_l._x,_l._y,_dim);
  }
  void updateDist(const Node<sizeType>& node) {
    scalar s;
    const Vec3i& iss=_iss[node._cell];
    if(_dim == 3) {
      Triangle t(_vss[iss[0]],_vss[iss[1]],_vss[iss[2]]);
      if(t.intersect(_l,s) && s < _s)
        _s=s;
    } else {
      LineSeg l(_vss[iss[0]],_vss[iss[1]]);
      if(l.intersect(_l,s) && s < _s)
        _s=s;
    }
  }
  const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& _vss;
  const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& _iss;
  LineSeg _l;
  sizeType _dim;
  scalar _s;
};
EIGEN_DEVICE_FUNC ObjMeshGeomCell::ObjMeshGeomCell()
  :StaticGeomCell(typeid(ObjMeshGeomCell).name()),_depth(0.0f) {}
EIGEN_DEVICE_FUNC ObjMeshGeomCell::ObjMeshGeomCell(const Mat4& trans,const ObjMesh& mesh,scalar depth,bool buildBVHAsWell,bool bottomUp)
  :StaticGeomCell(trans,mesh.getDim(),typeid(ObjMeshGeomCell).name()),_depth(depth)
{
  _vss=mesh.getV();
  _iss=mesh.getI();
  build(buildBVHAsWell,bottomUp&&mesh.getDim()==3);
  //writeBVHByLevel<sizeType,BBox<scalar>>(_bvh,-1);
  if(_depth == 0.0f)
    _depth=mesh.getBB().getExtent().norm();
}
bool ObjMeshGeomCell::read(std::istream& is,IOData* dat)
{
  StaticGeomCell::read(is,dat);
  _grid.read(is);
  readBinaryData(_depth,is);
  return is.good();
}
bool ObjMeshGeomCell::write(std::ostream& os,IOData* dat) const
{
  StaticGeomCell::write(os,dat);
  _grid.write(os);
  writeBinaryData(_depth,os);
  return os.good();
}
std::shared_ptr<SerializableBase> ObjMeshGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new ObjMeshGeomCell(*this));
}
scalar ObjMeshGeomCell::depth() const
{
  return _depth;
}
scalar& ObjMeshGeomCell::depth()
{
  return _depth;
}
void ObjMeshGeomCell::calcMinDist2D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist,Vec2i& feat) const
{
  Vec3 cpTmp,b;
  scalar distNew;
  LineSegTpl<scalar> l(_vss[I[0]],_vss[I[1]]);
  static const scalarF eps=1E-5f;

  l.calcPointDist(pt,distNew,cpTmp,b);
  distNew=std::sqrt(distNew);
  *minDist=std::min(*minDist,distNew);

  Vec2i featCurr(-1,-1);
  //detect vertex feature
  for(sizeType d=0; d<2; d++)
    if(b[d]>1-eps)
      featCurr=Vec2i(I[d],-1);
  //update distance
  if(featCurr[0]>=0 && feat==featCurr) {
    scalar align=std::abs((pt-cp).dot(n));
    scalar alignCurr=std::abs((pt-cpTmp).dot(l.normal()));
    //this is a better element to tell inside from outside
    if(alignCurr>align) {
      //INFOV("Using feature rule: (%d,%d)!",feat[0],feat[1])
      cp=cpTmp;
      dist=distNew;
      n=l.normal();
    }
  } else if(distNew < dist) {
    cp=cpTmp;
    dist=distNew;
    n=l.normal();
  }
}
void ObjMeshGeomCell::calcMinDist3D(const Vec3i& I,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist,Vec2i& feat) const
{
  Vec3 cpTmp,b;
  scalar distNew;
  Triangle t(_vss[I[0]],_vss[I[1]],_vss[I[2]]);
  static const scalarF eps=1E-5f;

  t.calcPointDist(pt,distNew,cpTmp,b);
  distNew=std::sqrt(distNew);
  *minDist=std::min(*minDist,distNew);

  Vec2i featCurr(-1,-1);
  //detect vertex feature
  for(sizeType d=0; d<3; d++)
    if(b[d]>1-eps)
      featCurr=Vec2i(I[d],-1);
  //detect edge feature
  if(featCurr[0]==-1) {
    for(sizeType d=0; d<3; d++)
      if(b[d]<eps) {
        featCurr=Vec2i(I[(d+1)%3],I[(d+2)%3]);
        if(featCurr[0]>featCurr[1])
          std::swap(featCurr[0],featCurr[1]);
      }
  }
  //update distance
  if(featCurr[0]>=0 && feat==featCurr) {
    scalar align=std::abs((pt-cp).dot(n));
    scalar alignCurr=std::abs((pt-cpTmp).dot(t.normal()));
    //this is a better element to tell inside from outside
    if(alignCurr>align) {
      //INFOV("Using feature rule: (%d,%d)!",feat[0],feat[1])
      cp=cpTmp;
      dist=distNew;
      n=t.normal();
    }
  } else if(distNew<dist) {
    feat=featCurr;
    cp=cpTmp;
    dist=distNew;
    n=t.normal();
  }
}
void ObjMeshGeomCell::buildLevelSet(scalar cellSz,scalar off,scalar eps)
{
  ImplicitFuncOffset offset;
  offset._off-=off;
  offset._inner.reset(new ImplicitFuncMeshRef(*this,eps));
  _grid=ImplicitFuncReinit(cellSz,offset,true)._ls;

  _grid.add(off);
  //GridOp<scalar,scalar>::write3DScalarGridVTK("./levelset.vtk",_grid);
}
scalar ObjMeshGeomCell::distLevelSet(const Vec3& pos) const
{
  return _grid.sampleSafe(transformHomo<scalar>(_invT,pos));
}
const ScalarField& ObjMeshGeomCell::getLevelSet() const
{
  return _grid;
}
void ObjMeshGeomCell::setRes(sizeType res)
{
  FUNCTION_NOT_IMPLEMENTED
}
//helper
ObjMeshGeomCell::ObjMeshGeomCell(const std::string& name):StaticGeomCell(name) {}
ObjMeshGeomCell::ObjMeshGeomCell(const Mat4& T,sizeType dim,const std::string& name):StaticGeomCell(T,dim,name) {}
void ObjMeshGeomCell::getMeshInner(ObjMesh& mesh) const
{
  mesh.getV()=_vss;
  mesh.getI()=_iss;
  mesh.setDim((int)_dim);
  mesh.smooth();
}
class CallbackMesh
{
public:
  CallbackMesh(const ObjMeshGeomCell& cell,sizeType dim):_cell(cell),_iss(cell.iss()),_dim(dim),_feat(-1,-1) {}
  void updateDist(const Node<sizeType>& node,const Vec3& pt,Vec3& cp,Vec3& n,scalar& dist,scalar* minDist) {
    if(_dim == 2)
      _cell.calcMinDist2D(_iss[node._cell],pt,cp,n,dist,minDist,_feat);
    else
      _cell.calcMinDist3D(_iss[node._cell],pt,cp,n,dist,minDist,_feat);
  }
  scalar depth() const {
    return _cell.depth();
  }
  const ObjMeshGeomCell& _cell;
  const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& _iss;
  const sizeType _dim;
  Vec2i _feat;
};
DEVICE_ONLY_FUNC BBox<scalar> ObjMeshGeomCell::getBBInner() const
{
  if(!_grid.data().empty())
    return _grid.getBB();
  return _bvh.back()._bb;
}
DEVICE_ONLY_FUNC bool ObjMeshGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  if(!_grid.data().empty()) {
    scalar phi=_grid.sampleSafe(pt);
    if(phi > 0)
      return false;
    else {
      n=-_grid.sampleSafeGrad(pt).normalized()*phi;
      return true;
    }
  }

  Vec3 cp;
  CallbackMesh cb(*this,_dim);
  scalar dist=ScalarUtil<scalar>::scalar_max(),minDist=dist;
  BVHQuery<sizeType>(_bvh,_dim,-1).pointDistQuery(pt,cb,cp,n,dist,&minDist);
  cp-=pt;
  if(cp.dot(n) < 1E-6f)
    return false;
  n=cp;
  return dist < ScalarUtil<scalar>::scalar_max();
}
DEVICE_ONLY_FUNC bool ObjMeshGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  if(!_grid.data().empty()) {
    scalar phi=_grid.sampleSafe(pt);
    n=-_grid.sampleSafeGrad(pt).normalized()*phi;
    return phi < 0;
  }

  Vec3 cp,nor;
  scalar dist=ScalarUtil<scalar>::scalar_max(),minDist=dist;
  cp.block(0,0,_dim,1).setConstant(dist);
  CallbackMesh cb(*this,_dim);
  BVHQuery<sizeType>(_bvh,_dim,-1).pointDistQuery(pt,cb,cp,nor,dist,&minDist);
  if(normal)*normal=nor;
  n=cp-pt;

  if(dist < _depth) {
    n*=minDist/std::max<scalar>(n.norm(),ScalarUtil<scalar>::scalar_eps());
    return n.dot(nor) > 0.0f;
  } else {
    n.block(0,0,_dim,1).setConstant(_depth);
    if(minDist < _depth)
      n*=minDist/std::max<scalar>(n.norm(),ScalarUtil<scalar>::scalar_eps());
    return false;
  }
}
DEVICE_ONLY_FUNC scalar ObjMeshGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  LineCallback cb(_vss,_iss,LineSeg(x0,x0+dir),_dim);
  BVHQuery<sizeType>(_bvh,_dim,-1).pointQuery(cb);
  return cb._s;
}
