#include "../CollisionDetection.h"
#include "../MakeMesh.h"
#include "../CameraModel.h"
#include "../Heap.h"
#include "BVHBuilder.h"
#include "StaticGeom.h"
#include "StaticGeomCell.h"
#include <experimental/filesystem>
#include <unordered_map>

USE_PRJ_NAMESPACE

//Geom
EIGEN_DEVICE_FUNC StaticGeomCell::StaticGeomCell(const std::string& type):Serializable(type) {}
EIGEN_DEVICE_FUNC StaticGeomCell::StaticGeomCell(const Mat4& T,sizeType dim,const std::string& type)
  :Serializable(type),_T(T),_invT(T.inverse()),_dim(dim),_res(16) {}
void StaticGeomCell::subdivideMesh(sizeType nrSubd)
{
  if(nrSubd <= 0)
    return;
  ObjMesh mesh;
  getMesh(mesh,true);
  mesh.smooth();
  mesh.subdivide(nrSubd);
  _tss.clear();
  _vss=mesh.getV();
  _iss=mesh.getI();
  _bvh.clear();
}
void StaticGeomCell::generateUV(scalar scale)
{
  ObjMesh mesh;
  getMesh(mesh,true);
  generateUVInner(mesh,scale);
  _tss=mesh._tex._uv;
  _vss=mesh.getV();
  _iss=mesh.getI();
}
void StaticGeomCell::getMesh(ObjMesh& mesh,bool ref,bool render) const
{
  if(render) {
    getMeshInner(mesh);
  } else {
    if(_tss.size() == _vss.size()) {
      mesh._tex._uv=_tss;
      mesh._tex._fuv=_iss;
    }
    mesh.getV()=_vss;
    mesh.getI()=_iss;
  }
  mesh.setDim((int)_dim);
  if(!ref) {
    mesh.getT()=_T.block<3,3>(0,0);
    mesh.getPos()=_T.block<3,1>(0,3);
    mesh.applyTrans(Vec3::Zero());
  }
}
DEVICE_ONLY_FUNC BBox<scalar> StaticGeomCell::getBB(bool ref) const
{
  Vec3 pt;
  BBox<scalar> tmp=getBBInner(),ret;
  if(ref)
    return tmp;
  for(sizeType x=0; x<2; x++)
    for(sizeType y=0; y<2; y++)
      for(sizeType z=0; z<2; z++) {
        pt[0]=(x==0) ? tmp._minC[0] : tmp._maxC[0];
        pt[1]=(y==0) ? tmp._minC[1] : tmp._maxC[1];
        pt[2]=(z==0) ? tmp._minC[2] : tmp._maxC[2];
        ret.setUnion(transformHomo<scalar>(_T,pt));
      }
  return ret;
}
DEVICE_ONLY_FUNC bool StaticGeomCell::dist(const Vec3& pt,Vec3& n) const
{
  Vec3 pt0=transformHomo<scalar>(_invT,pt);
  if(distInner(pt0,n)) {
    n=_T.block<3,3>(0,0)*n;
    return true;
  }
  return false;
}
DEVICE_ONLY_FUNC bool StaticGeomCell::closest(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  Vec3 pt0=transformHomo<scalar>(_invT,pt);
  bool inside=closestInner(pt0,n,normal);
  n=_T.block<3,3>(0,0)*n;
  if(normal)
    *normal=_T.block<3,3>(0,0)**normal;
  return inside;
}
DEVICE_ONLY_FUNC scalar StaticGeomCell::rayQuery(Vec3 x0,Vec3 dir) const
{
  x0=transformHomo<scalar>(_invT,x0);
  dir=(_invT.block<3,3>(0,0)*dir).eval();
  return rayQueryInner(x0,dir);
}
bool StaticGeomCell::read(std::istream& is,IOData* dat)
{
  readBinaryData(_tss,is);
  readBinaryData(_vss,is);
  readBinaryData(_iss,is);
  readBinaryData(_bvh,is);

  readBinaryData(_T,is);
  readBinaryData(_invT,is);
  readBinaryData(_dim,is);
  readBinaryData(_index,is);
  return is.good();
}
bool StaticGeomCell::write(std::ostream& os,IOData* dat) const
{
  writeBinaryData(_tss,os);
  writeBinaryData(_vss,os);
  writeBinaryData(_iss,os);
  writeBinaryData(_bvh,os);

  writeBinaryData(_T,os);
  writeBinaryData(_invT,os);
  writeBinaryData(_dim,os);
  writeBinaryData(_index,os);
  return os.good();
}
const std::vector<Node<sizeType>>& StaticGeomCell::bvh() const
{
  return _bvh;
}
const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& StaticGeomCell::vss() const
{
  return _vss;
}
const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& StaticGeomCell::iss() const
{
  return _iss;
}
void StaticGeomCell::debugDistQuery(bool closestTest,scalar nScale)
{
#define NR_TEST 1000
  //test dist query
  Vec3 n,normal;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > pss,lss;
  std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
  std::vector<scalar> css;
  BBox<scalar> bb=getBB().enlargeEps(1);
  for(sizeType i=0; i<NR_TEST; i++) {
    Vec3 pt=Vec3::Zero();
    for(sizeType d=0; d<_dim; d++)
      pt[d]=RandEngine::randR(bb._minC[d],bb._maxC[d]);
    pss.push_back(Vec3i::Constant((sizeType)vss.size()-1));
    vss.push_back(pt);
    css.push_back(0);
    if(closestTest) {
      bool inside=closest(pt,n,&normal);
      lss.push_back(Vec3i((sizeType)vss.size()-1,(sizeType)vss.size(),0));
      lss.push_back(Vec3i((sizeType)vss.size()-1,(sizeType)vss.size()+1,0));

      vss.push_back(pt+n);
      css.push_back(1);
      vss.push_back(pt+normal*nScale);
      css.push_back(inside ? -2 : 2);
    } else if(dist(pt,n)) {
      lss.push_back(Vec3i((sizeType)vss.size()-1,(sizeType)vss.size(),0));
      vss.push_back(pt+n);
      css.push_back(1);
    }
  }
  VTKWriter<scalar> os("DistTest","./testDist/result"+std::string(closestTest ? "Closest" : "Dist")+".vtk",true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCustomPointData("Color",css.begin(),css.end());
  os.appendCells(pss.begin(),pss.end(),VTKWriter<scalar>::POINT);
  os.appendCells(lss.begin(),lss.end(),VTKWriter<scalar>::LINE);
#undef NR_TEST
}
void StaticGeomCell::debugDistQuery(scalar nScale)
{
  std::experimental::filesystem::v1::create_directory("./testDist/");
  //write geometry
  ObjMesh m;
  getMesh(m);
  m.writeVTK("./testDist/geom.vtk",true);
  //test
  debugDistQuery(true,nScale);
  debugDistQuery(false,nScale);
}
const Mat4& StaticGeomCell::getInvT() const
{
  return _invT;
}
const Mat4& StaticGeomCell::getT() const
{
  return _T;
}
void StaticGeomCell::setT(const Mat4& T)
{
  _T=T;
  _invT=T.inverse();
}
sizeType StaticGeomCell::getRes() const
{
  return _res;
}
void StaticGeomCell::setRes(sizeType res)
{
  _res=res;
  build(false);
}
sizeType StaticGeomCell::dim() const
{
  return _dim;
}
void StaticGeomCell::getMeshInner(ObjMesh& mesh) const
{
  FUNCTION_NOT_IMPLEMENTED
}
DEVICE_ONLY_FUNC BBox<scalar> StaticGeomCell::getBBInner() const
{
  FUNCTION_NOT_IMPLEMENTED
  return BBox<scalar>();
}
DEVICE_ONLY_FUNC bool StaticGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  return closestInner(pt,n);
}
DEVICE_ONLY_FUNC bool StaticGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  FUNCTION_NOT_IMPLEMENTED
  return false;
}
DEVICE_ONLY_FUNC scalar StaticGeomCell::rayQueryInner(const Vec3& x0,const Vec3& dir) const
{
  FUNCTION_NOT_IMPLEMENTED
  return 1;
}
struct StaticGeomCellEdgeHash {
  size_t operator()(const Vec2i& key) const {
    std::hash<sizeType> h;
    return h(key[0])+h(key[1]);
  }
};
void StaticGeomCell::buildBVHBottomUp()
{
  std::unordered_map<Vec2i,Vec2i,StaticGeomCellEdgeHash> edgeMap;
  for(sizeType i=0; i<(sizeType)_iss.size(); i++) {
    ASSERT_MSG(_iss[i][2]!=-1,"You cannot use BottomUp BVH building in 3D Meshes!")
    for(sizeType d=0; d<3; d++) {
      //edge index
      Vec2i e(_iss[i][d],_iss[i][(d+1)%3]);
      if(e[0]>e[1])
        std::swap(e[0],e[1]);
      //insert edge
      std::unordered_map<Vec2i,Vec2i,StaticGeomCellEdgeHash>::iterator it=edgeMap.find(e);
      if(it==edgeMap.end())
        edgeMap[e]=Vec2i(i,-1);
      else {
        ASSERT_MSG(it->second[1]==-1,"Non-manifold mesh detected!")
        it->second[1]=i;
      }
    }
  }
  //initialize hash
  std::vector<scalar> cost;
  std::vector<sizeType> heap;
  std::vector<sizeType> heapOffsets;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> ess;
  for(std::unordered_map<Vec2i,Vec2i,StaticGeomCellEdgeHash>::const_iterator beg=edgeMap.begin(),end=edgeMap.end(); beg!=end; beg++)
  {
    heapOffsets.push_back(-1);
    ess.push_back(beg->second);
    BBox<scalar> bb=_bvh[beg->second[0]]._bb;
    if(beg->second[1]>=0)
      bb.setUnion(_bvh[beg->second[1]]._bb);
    scalar c=SurfaceArea<3>::area(bb);
    cost.push_back(c);
  }
  for(sizeType i=0; i<(sizeType)ess.size(); i++)
    pushHeapDef(cost,heapOffsets,heap,i);
  //merge BVH
  sizeType err;
  while(!heap.empty()) {
    sizeType i=popHeapDef(cost,heapOffsets,heap,err);
    sizeType t0=ess[i][0],t1=ess[i][1];
    //boundary edge
    if(t1==-1)
      continue;
    //find parent
    while(_bvh[t0]._parent>=0)
      t0=_bvh[t0]._parent;
    while(_bvh[t1]._parent>=0)
      t1=_bvh[t1]._parent;
    //check already merged
    if(t0==t1)
      continue;
    //merge
    BBox<scalar> bb=_bvh[t0]._bb;
    bb.setUnion(_bvh[t1]._bb);
    scalar c=SurfaceArea<3>::area(bb);
    if(c>cost[i]) {
      cost[i]=c;
      pushHeapDef(cost,heapOffsets,heap,i);
    } else {
      Node<sizeType> n;
      n._l=t0;
      n._r=t1;
      n._parent=-1;
      n._cell=-1;
      n._bb=bb;
      n._nrCell=_bvh[n._l]._nrCell+_bvh[n._r]._nrCell;
      _bvh[t0]._parent=(sizeType)_bvh.size();
      _bvh[t1]._parent=(sizeType)_bvh.size();
      _bvh.push_back(n);
    }
  }
  ASSERT_MSG(_bvh.size()==_iss.size()*2-1,"Multi-component mesh detected!")
}
void StaticGeomCell::build(bool buildBVHAsWell,bool bottomUp)
{
  ObjMesh mesh;
  getMeshInner(mesh);
  _vss=mesh.getV();
  _iss=mesh.getI();
  if(buildBVHAsWell) {
    //bvh triangle
    _bvh.resize(_iss.size());
    for(sizeType i=0; i<(sizeType)_bvh.size(); i++) {
      Node<sizeType>& n=_bvh[i];
      n._nrCell=1;
      n._cell=i;
      n._bb.reset();
      for(sizeType v=0; v<_dim; v++)
        n._bb.setUnion(_vss[_iss[i][v]]);
      n._bb.enlarged(0.01f*n._bb.getExtent().maxCoeff());
    }
    if(bottomUp)
      buildBVHBottomUp();
    else buildBVH<sizeType>(_bvh,_dim,-1);
  } else {
    _bvh.clear();
  }
}
void StaticGeomCell::generateUVInner(ObjMesh& mesh,scalar scale) const
{
  FUNCTION_NOT_IMPLEMENTED
}

//helper
template <typename T>
typename ScalarUtil<T>::ScalarMat3 makeRandomRotation(sizeType dim)
{
  return dim == 3 ?
         Eigen::AngleAxis<scalar>(RandEngine::randR(0,M_PI*2),Vec3::Random().normalized()).toRotationMatrix() :
         Eigen::AngleAxis<scalar>(RandEngine::randR(0,M_PI*2),Vec3::UnitZ()).toRotationMatrix();
}
//Static Geometry
PRJ_BEGIN
struct GeomCallback {
  GeomCallback(const Vec3& pt,Vec3& n):_pt(pt),_dist(ScalarUtil<scalar>::scalar_max()),_n(n) {}
  virtual ~GeomCallback() {}
  virtual bool validNode(const Node<std::shared_ptr<StaticGeomCell> >& node) {
    return node._bb.contain(_pt);
  }
  virtual void updateDist(const Node<std::shared_ptr<StaticGeomCell> >& node) {
    Vec3 n;
    scalar dist;
    if(node._cell->dist(_pt,n) && (dist=n.norm()) < _dist) {
      _cell=node._cell;
      _dist=dist;
      _n=n;
    }
  }
  const Vec3& _pt;
  std::shared_ptr<StaticGeomCell> _cell;
  scalar _dist;
  Vec3& _n;
};
struct GeomClosestCallback : public GeomCallback {
  GeomClosestCallback(const Vec3& pt,Vec3& n,Vec3* normal):GeomCallback(pt,n),_normal(normal) {}
  virtual bool validNode(const Node<std::shared_ptr<StaticGeomCell> >& node) {
    return node._bb.distTo(_pt) < _dist;
  }
  virtual void updateDist(const Node<std::shared_ptr<StaticGeomCell> >& node) {
    Vec3 n,normal;
    scalar dist;
    node._cell->closest(_pt,n,_normal?&normal:NULL);
    if((dist=n.norm()) < _dist) {
      _cell=node._cell;
      _dist=dist;
      _n=n;
      if(_normal)
        *_normal=normal;
    }
  }
  Vec3* _normal;
};
struct LineCallback {
  LineCallback(const Vec3& x0,Vec3& dir,sizeType dim):_x0(x0),_dir(dir),_dim(dim) {}
  bool validNode(const Node<std::shared_ptr<StaticGeomCell> >& node) {
    return true;//node._bb.intersect(_x0,_x0+_dir,_dim);
  }
  void updateDist(const Node<std::shared_ptr<StaticGeomCell> >& node) {
    scalar s=node._cell->rayQuery(_x0,_dir);
    if(s < 1) {
      _dir*=s;
      _cell=node._cell;
    }
  }
  Vec3 _x0;
  Vec3& _dir;
  sizeType _dim;
  std::shared_ptr<StaticGeomCell> _cell;
};
struct VertexCallback {
  VertexCallback(StaticGeomCallback& cb):_cb(cb) {}
  bool validNode(const Node<sizeType>& node) {
    if(_gA->dim() == 2) {
      OBB2D obb(_BTA.block<3,3>(0,0),_BTA.block<3,1>(0,3),node._bb);
      return obb.intersect(_bbA);
    } else {
      OBB3D obb(_BTA.block<3,3>(0,0),_BTA.block<3,1>(0,3),node._bb);
      return obb.intersect(_bbA);
    }
    return true;
  }
  void updateDist(const Node<sizeType>& node) {
    Vec3 n,vRef,vG;
    sizeType dim=_gA->dim();
    const Vec3i& I=_gB->_iss[node._cell];
    for(sizeType v=0; v<dim; v++)
      if(!_mask[I[v]]) {
        _mask[I[v]]=true;
        vG=transformHomo<scalar>(_gB->getT(),_gB->_vss[I[v]]);
        if(_gA->dist(vG,n))
          _cb.onCollideVertex(vG,n,_gB);
      }
  }
  void onCell(const Node<std::shared_ptr<StaticGeomCell> >& nA,
              const Node<std::shared_ptr<StaticGeomCell> >& nB) {
    _gA=nA._cell;
    _gB=nB._cell;
    _mask.assign(_gB->_vss.size(),false);
    _bbA=_gA->getBBInner();
    _BTA=_gA->getInvT()*_gB->getT();

    BVHQuery<sizeType,BBox<scalar> > queryNarrow(nB._cell->_bvh,nB._cell->dim(),-1);
    queryNarrow.pointQuery(*this);
  }
  StaticGeomCallback& _cb;
  std::shared_ptr<StaticGeomCell> _gA,_gB;
  std::vector<bool> _mask;
  BBox<scalar> _bbA;
  Mat4 _BTA;
};
class DebugStaticGeomCallback : public StaticGeomCallback
{
public:
  DebugStaticGeomCallback(const std::string& path):_os("DebugCallback",path,true) {}
  void onCollideVertex(const Vec3& x,const Vec3& n,std::shared_ptr<StaticGeomCell> c) {
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
    vss.push_back(x);
    vss.push_back(x+n);

    _os.setRelativeIndex();
    _os.appendPoints(vss.begin(),vss.end());
    _os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                    VTKWriter<scalar>::IteratorIndex<Vec3i>(1,2,0),
                    VTKWriter<scalar>::LINE,true);
  }
  VTKWriter<scalar> _os;
};
PRJ_END
EIGEN_DEVICE_FUNC StaticGeom::StaticGeom():Serializable(typeid(StaticGeom).name())
{
  _bvh.reset(new std::vector<Node<std::shared_ptr<StaticGeomCell> > >);
}
EIGEN_DEVICE_FUNC StaticGeom::StaticGeom(sizeType dim):Serializable(typeid(StaticGeom).name()),_dim(dim)
{
  _bvh.reset(new std::vector<Node<std::shared_ptr<StaticGeomCell> > >);
}
const std::vector<Node<std::shared_ptr<StaticGeomCell> > >& StaticGeom::getBVH() const
{
  return *_bvh;
}
std::vector<Node<std::shared_ptr<StaticGeomCell> > >& StaticGeom::getBVH()
{
  return *_bvh;
}
const StaticGeomCell& StaticGeom::getG(sizeType i) const
{
  return *(_css[i]);
}
std::shared_ptr<StaticGeomCell> StaticGeom::getGPtr(sizeType i) const
{
  return _css[i];
}
StaticGeomCell& StaticGeom::getG(sizeType i)
{
  return *(_css[i]);
}
sizeType StaticGeom::nrG() const
{
  return (sizeType)_css.size();
}
sizeType StaticGeom::depth() const
{
  if(!_bvh)
    return 0;
  BVHQuery<std::shared_ptr<StaticGeomCell>,BBox<scalar> > handler(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  return handler.depth();
}
void StaticGeom::clear()
{
  _css.clear();
  _bvh->clear();
}
void StaticGeom::assemble()
{
  _bvh->clear();
  for(sizeType i=0; i<(sizeType)_css.size(); i++) {
    _css[i]->_index=i;
    Node<std::shared_ptr<StaticGeomCell> > n;
    n._l=n._r=n._parent=-1;
    n._nrCell=1;
    n._cell=_css[i];
    n._bb=n._cell->getBB();
    _bvh->push_back(n);
  }
  buildBVH(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
}
void StaticGeom::parityCheck()
{
  if(!_bvh)
    return;
  BVHQuery<std::shared_ptr<StaticGeomCell>,BBox<scalar> > handler(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  handler.parityCheck();
}
bool StaticGeom::update(scalar expand,bool dynamic)
{
  if(!_bvh)
    return false;
  BVHQuery<std::shared_ptr<StaticGeomCell>,BBox<scalar> > handler(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  //compact
  handler.compact();
  //leaf
  std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar> > >& bvh=*_bvh;
  for(sizeType i=0; i<(sizeType)bvh.size(); i++)
    if(bvh[i]._cell)
      bvh[i]._bb=bvh[i]._cell->getBB();
  //non-leaf
  return handler.updateBVH(expand,dynamic);
}
bool StaticGeom::dist(const Vec3& pt,Vec3& n,std::shared_ptr<StaticGeomCell>& cell) const
{
  BVHQuery<std::shared_ptr<StaticGeomCell> > query(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  GeomCallback g(pt,n);
  query.pointQuery(g);
  cell=g._cell;
  return g._dist < ScalarUtil<scalar>::scalar_max();
}
void StaticGeom::closest(const Vec3& pt,Vec3& n,Vec3* normal,std::shared_ptr<StaticGeomCell>& cell) const
{
  BVHQuery<std::shared_ptr<StaticGeomCell> > query(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  GeomClosestCallback g(pt,n,normal);
  query.pointQuery(g);
  cell=g._cell;
}
bool StaticGeom::rayQuery(const Vec3& pt0,Vec3& dir,std::shared_ptr<StaticGeomCell>& cell,Vec3& r) const
{
  BVHQuery<std::shared_ptr<StaticGeomCell> > query(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  LineCallback g(pt0,dir,_dim);
  query.pointQuery(g);
  if(g._cell) {
    cell=g._cell;
    r=transformHomo<scalar>(cell->getInvT(),g._x0+g._dir);
    return true;
  } else return false;
}
void StaticGeom::collideVertex(const StaticGeom& other,StaticGeomCallback& cb) const
{
  forceBuild();
  other.forceBuild();
  BVHQuery<std::shared_ptr<StaticGeomCell>,BBox<scalar> > query(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  BVHQuery<std::shared_ptr<StaticGeomCell>,BBox<scalar> > queryOther(*(other._bvh),_dim,std::shared_ptr<StaticGeomCell>());
  VertexCallback vcb(cb);
  query.interBodyQuery(queryOther,vcb);
}
void StaticGeom::forceBuild() const
{
  for(sizeType i=0; i<(sizeType)_css.size(); i++) {
    StaticGeomCell& c=*(_css[i]);
    sizeType nrI=(sizeType)c._iss.size();
    if((sizeType)c._bvh.size() != (nrI==0?0:(nrI*2-1)))
      c.build(true);
  }
}
void StaticGeom::setRes(sizeType res)
{
  for(sizeType i=0; i<(sizeType)_css.size(); i++)
    _css[i]->setRes(res);
}
//remove geometry
void StaticGeom::removeGeomCell(std::shared_ptr<StaticGeomCell> c)
{
  std::vector<std::shared_ptr<StaticGeomCell> >::iterator it=std::find(_css.begin(),_css.end(),c);
  if(it == _css.end())
    return;
  _css.erase(it);
  BVHQuery<std::shared_ptr<StaticGeomCell>,BBox<scalar> > handler(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  for(sizeType i=0; i<(sizeType)_bvh->size(); i++)
    if(_bvh->at(i)._cell == c) {
      handler.removeLeaf(i);
      return;
    }
}
//add geometry
void StaticGeom::addGeomCell(std::shared_ptr<StaticGeomCell> c)
{
  _css.push_back(c);
  if(!_bvh)
    return;

  BVHQuery<std::shared_ptr<StaticGeomCell>,BBox<scalar> > handler(*_bvh,_dim,std::shared_ptr<StaticGeomCell>());
  handler.insertLeaf(c,handler.findLeaf(c->getBB()));
}
void StaticGeom::addGeomBox(const Mat4& trans,const BBox<scalar>& bb,scalar depth)
{
  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=(bb._maxC+bb._minC)/2.0f;
  addGeomBox(trans*T,bb.getExtent()*0.5f,depth);
}
void StaticGeom::addGeomBox(const Mat4& trans,const Vec3& ext,scalar depth)
{
  depth=std::max<scalar>(depth,0.0f);
  addGeomCell(std::shared_ptr<StaticGeomCell>(new BoxGeomCell(trans,_dim,ext,depth)));
}
void StaticGeom::addGeomBox(const OBBTpl<scalar,2>& obb,scalar depth)
{
  ASSERT(_dim == 2)
  Mat4 m=Mat4::Identity();
  m.block<2,2>(0,0)=obb._rot;
  m.block<2,1>(0,3)=obb._trans;
  addGeomBox(m,Vec3(obb._ext[0],obb._ext[1],0.0f),depth);
}
void StaticGeom::addGeomBox(const OBBTpl<scalar,3>& obb,scalar depth)
{
  ASSERT(_dim == 3)
  Mat4 m=Mat4::Identity();
  m.block<3,3>(0,0)=obb._rot;
  m.block<3,1>(0,3)=obb._trans;
  addGeomBox(m,obb._ext,depth);
}
void StaticGeom::addGeomSphericalBox(const Mat4& trans,const Vec4& ext)
{
  ASSERT(_dim == 3)
  addGeomCell(std::shared_ptr<StaticGeomCell>(new SphericalBoxGeomCell(trans,_dim,ext)));
}
void StaticGeom::addGeomCylinder(const Mat4& trans,scalar rad,scalar y,bool capsule)
{
  if(!capsule)
    addGeomCell(std::shared_ptr<StaticGeomCell>(new CylinderGeomCell(trans,_dim,rad,y)));
  else addGeomCell(std::shared_ptr<StaticGeomCell>(new CapsuleGeomCell(trans,_dim,rad,y)));
}
void StaticGeom::addGeomPlane(const Mat4& trans,const Vec4& plane,scalar ext)
{
  scalar alpha=-plane[3]/plane.block<3,1>(0,0).squaredNorm();
  Vec3 p0=plane.block<3,1>(0,0)*alpha;
  Quatf q;
  q.setFromTwoVectors(Vec3::Unit(1).cast<scalarF>(),plane.block<3,1>(0,0).normalized().cast<scalarF>());

  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=p0-plane.block<3,1>(0,0).normalized()*ext;
  T.block<3,3>(0,0)=q.cast<scalar>().toRotationMatrix();
  if(_dim == 3)
    addGeomCell(std::shared_ptr<StaticGeomCell>(new CylinderGeomCell(trans*T,_dim,ext,ext)));
  else addGeomCell(std::shared_ptr<StaticGeomCell>(new BoxGeomCell(trans*T,_dim,Vec3::Constant(ext))));
}
void StaticGeom::addGeomPlane(const Mat4& trans,const PlaneTpl<scalar>& plane,scalar ext)
{
  Vec4 p;
  p.block<3,1>(0,0)=plane._n;
  p[3]=-plane._n.dot(plane._x0);
  addGeomPlane(trans,p,ext);
}
void StaticGeom::addGeomSphere(const Vec3& ctr,scalar rad,scalar depth)
{
  Mat4 T=Mat4::Identity();
  T.block<3,1>(0,3)=ctr;
  depth=std::max<scalar>(depth,0.0f);
  addGeomCell(std::shared_ptr<StaticGeomCell>(new SphereGeomCell(T,_dim,rad,depth)));
}
void StaticGeom::addGeomMesh(const Mat4& trans,const ObjMesh& mesh,scalar depth)
{
  addGeomCell(std::shared_ptr<StaticGeomCell>(new ObjMeshGeomCell(trans,mesh,depth,true)));
}
void StaticGeom::addGeomMesh(const Mat4& trans,const std::string& path,scalar depth)
{
  ObjMesh mesh;
  std::ifstream is(path);
  mesh.read(is,false,false);
  mesh.smooth();
  mesh.makeUniform();
  addGeomMesh(trans,mesh,depth);
}
void StaticGeom::addGeomHeightField(sizeType dimH,scalar h0,scalar hr,scalar sz,scalar cellSz)
{
  addGeomCell(std::shared_ptr<StaticGeomCell>(new HeightFieldGeomCell(_dim,dimH,h0,hr,sz,cellSz)));
}
void StaticGeom::addGeomSolidBox(const Mat4& trans,const Vec3& ext,scalar thick)
{
  ObjMesh mesh;
  if(_dim == 2)MakeMesh::makeBox2D(mesh,ext,thick);
  else MakeMesh::makeBox3D(mesh,ext,thick);
  addGeomMesh(trans,mesh);
}
void StaticGeom::addGeomSolidSphere(const Mat4& trans,scalar rad,scalar thick)
{
  ObjMesh mesh;
  if(_dim == 2)MakeMesh::makeSphere2D(mesh,rad,32,thick);
  else MakeMesh::makeSphere3D(mesh,rad,32,thick);
  addGeomMesh(trans,mesh);
}
//IO
void StaticGeom::writeVTK(const std::string& path) const
{
  ObjMesh mesh;
  VTKWriter<scalar> os("Geom",path,true);
  std::experimental::filesystem::v1::path components=std::experimental::filesystem::v1::path(path).parent_path()/"geomComponents/";
  std::experimental::filesystem::v1::create_directory(components);
  for(sizeType i=0; i<(sizeType)_css.size(); i++) {
    _css[i]->getMesh(mesh);
    std::ostringstream oss;
    oss << components.string() << "/comp" << i << ".obj";
    mesh.write(oss.str());
    std::ofstream povSS(std::experimental::filesystem::v1::path(oss.str()).replace_extension(".pov"));
    mesh.writePov(povSS,false);
    std::ofstream spovSS(std::experimental::filesystem::v1::path(oss.str()).replace_extension(".spov"));
    mesh.writePov(spovSS,true);
    mesh.writeVTK(os,false,false);
  }
}
void StaticGeom::writeBVH() const
{
  typedef std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar> > > TV;
  writeBVHByLevel<std::shared_ptr<StaticGeomCell>,BBox<scalar>,TV>(*_bvh,std::shared_ptr<StaticGeomCell>());
}
bool StaticGeom::read(std::istream& is,IOData* dat)
{
  _css.clear();
  _bvh->clear();
  registerType(dat);
  readBinaryData(_dim,is);
  readBinaryData(_css,is,dat);
  readBinaryData(*_bvh,is,dat);
  return is.good();
}
bool StaticGeom::write(std::ostream& os,IOData* dat) const
{
  registerType(dat);
  writeBinaryData(_dim,os);
  writeBinaryData(_css,os,dat);
  writeBinaryData(*_bvh,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> StaticGeom::copy() const
{
  return std::shared_ptr<SerializableBase>(new StaticGeom);
}
bool StaticGeom::write(const std::shared_ptr<StaticGeomCell>& cell,std::ostream& os)
{
  std::shared_ptr<IOData> dat=getIOData();
  registerType(dat.get());
  writeBinaryData(cell,os,dat.get());
  return os.good();
}
bool StaticGeom::read(std::shared_ptr<StaticGeomCell>& cell,std::istream& is)
{
  std::shared_ptr<IOData> dat=getIOData();
  registerType(dat.get());
  readBinaryData(cell,is,dat.get());
  return is.good();
}
void StaticGeom::registerType(IOData* dat)
{
  NAMESPACE::registerType<BoxGeomCell>(dat);
  NAMESPACE::registerType<SphereGeomCell>(dat);
  NAMESPACE::registerType<CylinderGeomCell>(dat);
  NAMESPACE::registerType<SphericalBoxGeomCell>(dat);
  NAMESPACE::registerType<CapsuleGeomCell>(dat);
  NAMESPACE::registerType<ObjMeshGeomCell>(dat);
  NAMESPACE::registerType<HeightFieldGeomCell>(dat);
  NAMESPACE::registerType<CompositeGeomCell>(dat);
  NAMESPACE::registerType<TwoSphereMeshCell>(dat);
  NAMESPACE::registerType<ThreeSphereMeshCell>(dat);
}
void StaticGeom::debugRayQuery(const std::string& path,std::shared_ptr<StaticGeomCell> cell,sizeType nr)
{
  StaticGeom geom(cell->dim());
  for(sizeType i=0; i<nr; i++) {
    std::shared_ptr<StaticGeomCell> newCell=
      std::dynamic_pointer_cast<StaticGeomCell>(cell->copy());

    Mat4 T=Mat4::Identity();
    T.block<3,3>(0,0)=makeRandomRotation<scalar>(cell->dim());
    T.block(0,3,cell->dim(),1).setRandom();
    newCell->setT(T);
    geom.addGeomCell(newCell);
  }
  geom.assemble();
  {
    std::ostringstream oss;
    oss << "./data/geom" << path << ".vtk";
    geom.writeVTK(oss.str());
  }

  std::ostringstream oss;
  oss << "./data/ray" << path << ".vtk";
  VTKWriter<scalar> os("rayQuery",oss.str(),true);
  std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
  std::vector<scalar> css;

  Vec3 x,dir,r;
  BBox<scalar> bb=geom._bvh->back()._bb;
  scalar norm=bb.getExtent().norm();
  for(sizeType j=0; j<100; j++) {
    x.setZero();
    x.block(0,0,cell->dim(),1).setRandom();
    x=x.normalized()*norm*2;
    dir=Vec3::Zero()-x;

    std::shared_ptr<StaticGeomCell> ICell;
    geom.rayQuery(x,dir,ICell,r);

    vss.push_back(x);
    vss.push_back(x+dir);

    css.push_back(ICell ? 1.0f : 0.0f);
    css.push_back(ICell ? 1.0f : 0.0f);
  }
  os.appendPoints(vss.begin(),vss.end());
  os.appendCustomPointData("intersect",css.begin(),css.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>(vss.size()/2,2,0),
                 VTKWriter<scalar>::LINE);
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>(vss.size(),0,0),
                 VTKWriter<scalar>::POINT);
}
void StaticGeom::debugVertexQuery(const std::string& path,std::shared_ptr<StaticGeomCell> cell,sizeType nr)
{
  StaticGeom geomA(cell->dim());
  StaticGeom geomB(cell->dim());
  for(sizeType i=0; i<nr; i++) {
    std::shared_ptr<StaticGeomCell> newCell=
      std::dynamic_pointer_cast<StaticGeomCell>(cell->copy());

    Mat4 T=Mat4::Identity();
    T.block<3,3>(0,0)=makeRandomRotation<scalar>(cell->dim());
    T.block(0,3,cell->dim(),1).setRandom();
    newCell->setT(T);
    if(i >= nr/2)
      geomA.addGeomCell(newCell);
    else geomB.addGeomCell(newCell);
  }
  geomA.assemble();
  geomB.assemble();
  {
    std::ostringstream oss;
    oss << "./data/geomA" << path << ".vtk";
    geomA.writeVTK(oss.str());
  }
  {
    std::ostringstream oss;
    oss << "./data/geomB" << path << ".vtk";
    geomB.writeVTK(oss.str());
  }
  std::ostringstream oss;
  oss << "./data/vert" << path << ".vtk";
  DebugStaticGeomCallback dsgcb(oss.str());
  geomA.collideVertex(geomB,dsgcb);
}
void StaticGeom::debugVertexDistQuery(const std::string& path,std::shared_ptr<StaticGeomCell> cell,sizeType nr)
{
  std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > pss,lss;
  for(sizeType i=0; i<nr; i++) {
    Vec3 pt=Vec3::Zero(),n;
    pt.segment(0,cell->dim()).setRandom();
    pt*=cell->getBB().getExtent().maxCoeff();
    vss.push_back(pt);
    pss.push_back(Vec3i((sizeType)vss.size()-1,0,0));

    bool has=cell->dist(pt,n);
    if(has) {
      vss.push_back(pt+n);
      pss.push_back(Vec3i((sizeType)vss.size()-1,0,0));
      lss.push_back(Vec3i((sizeType)vss.size()-2,(sizeType)vss.size()-1,0));
    }
  }
  VTKWriter<scalar> os("VertDist","./data/vertDist"+path+".vtk",true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(pss.begin(),pss.end(),VTKWriter<scalar>::POINT);
  os.appendCells(lss.begin(),lss.end(),VTKWriter<scalar>::LINE);
}
void StaticGeom::debugRayVertexQuery(bool box,bool sphere,bool cylinder,bool sphericalBox,bool capsule,bool objMesh,bool height,bool twoSphere,bool threeSphere)
{
  ObjMesh mesh;
  if(std::experimental::filesystem::v1::exists("./data"))
    std::experimental::filesystem::v1::remove_all("./data");
  std::experimental::filesystem::v1::create_directory("./data");
  //box
  if(box) {
    std::shared_ptr<BoxGeomCell> box2(new BoxGeomCell(Mat4::Identity(),2,Vec3(0.1f,0.1f,0.0f)));
    debugRayQuery("box2",box2);
    debugVertexQuery("box2",box2);
    std::shared_ptr<BoxGeomCell> box3(new BoxGeomCell(Mat4::Identity(),3,Vec3(0.1f,0.1f,0.1f)));
    debugRayQuery("box3",box3);
    debugVertexQuery("box3",box3);
  }
  //sphere
  if(sphere) {
    std::shared_ptr<SphereGeomCell> sphere2(new SphereGeomCell(Mat4::Identity(),2,0.1f));
    debugRayQuery("sphere2",sphere2);
    debugVertexQuery("sphere2",sphere2);
    std::shared_ptr<SphereGeomCell> sphere3(new SphereGeomCell(Mat4::Identity(),3,0.1f));
    debugRayQuery("sphere3",sphere3);
    debugVertexQuery("sphere3",sphere3);
  }
  //cylinder
  if(cylinder) {
    std::shared_ptr<CylinderGeomCell> cylinder3(new CylinderGeomCell(Mat4::Identity(),3,0.1f,0.2f));
    debugRayQuery("cylinder3",cylinder3);
    debugVertexQuery("cylinder3",cylinder3);
  }
  //sphericalBox
  if(sphericalBox) {
    std::shared_ptr<SphericalBoxGeomCell> sbox2(new SphericalBoxGeomCell(Mat4::Identity(),2,Vec4(0.02f,0.02f,0.0f,0.1f)));
    //debugRayQuery("sphericalBox2",sbox2);
    debugVertexQuery("sphericalBox2",sbox2);
    std::shared_ptr<SphericalBoxGeomCell> sbox3(new SphericalBoxGeomCell(Mat4::Identity(),3,Vec4(0.02f,0.02f,0.02f,0.1f)));
    //debugRayQuery("sphericalBox3",sbox3);
    debugVertexQuery("sphericalBox3",sbox3);
  }
  //capsule
  if(capsule) {
    std::shared_ptr<CapsuleGeomCell> capsule2(new CapsuleGeomCell(Mat4::Identity(),2,0.1f,0.2f));
    debugRayQuery("capsule2",capsule2);
    debugVertexQuery("capsule2",capsule2);
    std::shared_ptr<CapsuleGeomCell> capsule3(new CapsuleGeomCell(Mat4::Identity(),3,0.1f,0.2f));
    debugRayQuery("capsule3",capsule3);
    debugVertexQuery("capsule3",capsule3);
  }
  //mesh
  if(objMesh && std::experimental::filesystem::v1::exists("./data/bunny.obj")) {
    std::ifstream mis("./data/bunny.obj");
    mesh.read(mis,false,false);
    mesh.getScale()=5.0f;
    mesh.applyTrans(Vec3::Zero());
    std::shared_ptr<ObjMeshGeomCell> mesh3(new ObjMeshGeomCell(Mat4::Identity(),mesh,1000.0f,false));
    debugRayQuery("mesh",mesh3,4);
    debugVertexQuery("mesh",mesh3,16);
  }
  //height
  if(height) {
    std::shared_ptr<HeightFieldGeomCell> height;
    //2D
    height.reset(new HeightFieldGeomCell(2,0,2,0.5f,30,1));
    height->getMesh(mesh);
    mesh.writeVTK("./data/height20.vtk",true);
    debugVertexDistQuery("height20",height,5000);

    height.reset(new HeightFieldGeomCell(2,1,2,0.5f,30,1));
    height->getMesh(mesh);
    mesh.writeVTK("./data/height21.vtk",true);
    debugVertexDistQuery("height21",height,5000);

    //3D
    height.reset(new HeightFieldGeomCell(3,0,2,0.5f,30,1));
    height->getMesh(mesh);
    mesh.writeVTK("./data/height30.vtk",true);
    debugVertexDistQuery("height30",height,5000);

    height.reset(new HeightFieldGeomCell(3,1,2,0.5f,30,1));
    height->getMesh(mesh);
    mesh.writeVTK("./data/height31.vtk",true);
    debugVertexDistQuery("height31",height,5000);

    height.reset(new HeightFieldGeomCell(3,2,2,0.5f,30,1));
    height->getMesh(mesh);
    mesh.writeVTK("./data/height32.vtk",true);
    debugVertexDistQuery("height32",height,5000);
  }
  //twoSphere
  if(twoSphere) {
    std::shared_ptr<TwoSphereMeshCell> twoSphere2;
    twoSphere2.reset(new TwoSphereMeshCell(2,Vec3::Random(),Vec3::Random(),RandEngine::randR01(),RandEngine::randR01()));
    debugRayQuery("twoSpehre2",twoSphere2);
    std::shared_ptr<TwoSphereMeshCell> twoSphere3;
    twoSphere3.reset(new TwoSphereMeshCell(3,Vec3::Random(),Vec3::Random(),RandEngine::randR01(),RandEngine::randR01()));
    debugRayQuery("twoSpehre3",twoSphere3);
  }
  //threeSphere
  if(threeSphere) {
    std::shared_ptr<ThreeSphereMeshCell> threeSphere2;
    threeSphere2.reset(new ThreeSphereMeshCell(2,Vec3::Random(),Vec3::Random(),Vec3::Random(),RandEngine::randR01(),RandEngine::randR01(),RandEngine::randR01()));
    debugRayQuery("threeSpehre2",threeSphere2);
    std::shared_ptr<ThreeSphereMeshCell> threeSphere3;
    threeSphere3.reset(new ThreeSphereMeshCell(3,Vec3::Random(),Vec3::Random(),Vec3::Random(),RandEngine::randR01(),RandEngine::randR01(),RandEngine::randR01()));
    debugRayQuery("threeSpehre3",threeSphere3);
  }
}
