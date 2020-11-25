#include "HeightFieldGeomCell.h"
#include "BVHBuilder.h"
#include <experimental/filesystem>

USE_PRJ_NAMESPACE

//HeightFieldGeomCell
//#define HEIGHT_FIELD_AS_MESH
#define GI(X,Y) Vec3i(id[0]+X,id[1]+Y,0).dot(stride)
EIGEN_DEVICE_FUNC HeightFieldGeomCell::HeightFieldGeomCell()
  :ObjMeshGeomCell(typeid(HeightFieldGeomCell).name()) {}
EIGEN_DEVICE_FUNC HeightFieldGeomCell::HeightFieldGeomCell(const Mat4& T,const ScalarField& h)
  :ObjMeshGeomCell(typeid(HeightFieldGeomCell).name()),_h(h)
{
  _T=T;
  _invT=_T.inverse();
  _dim=h.getDim()+1;
  ASSERT(!h.isCenter())

  build(true);
  //debugWrite();
}
EIGEN_DEVICE_FUNC HeightFieldGeomCell::HeightFieldGeomCell(sizeType dim,sizeType dimH,scalar h0,scalar hr,scalar sz,scalar cellSz)
  :ObjMeshGeomCell(typeid(HeightFieldGeomCell).name())
{
  _T=genT(dim-1,dimH);
  _invT=_T.inverse();
  _dim=dim;

  Vec3i id(0,0,0),nrCell(0,0,0);
  BBox<scalar> bb(Vec3::Zero(),Vec3::Zero());
  for(sizeType d=0; d<dim-1; d++) {
    bb._minC[d]=-sz;
    bb._maxC[d]=sz;
    nrCell[d]=(sizeType)(bb.getExtent()[d]/cellSz);
  }
  _h.reset(nrCell,bb,h0,false);
  for(id[0]=0; id[0]<_h.getNrPoint()[0]; id[0]++)
    for(id[1]=0; id[1]<_h.getNrPoint()[1]; id[1]++)
      _h.get(id)+=(RandEngine::randR01()*2-1)*hr;
  ASSERT(!_h.isCenter())

  build(true);
  //debugWrite();
}
bool HeightFieldGeomCell::read(std::istream& is,IOData* dat)
{
  ObjMeshGeomCell::read(is,dat);
  readBinaryData(_bb,is);
  _h.read(is);
  return is.good();
}
bool HeightFieldGeomCell::write(std::ostream& os,IOData* dat) const
{
  ObjMeshGeomCell::write(os,dat);
  writeBinaryData(_bb,os);
  _h.write(os);
  return os.good();
}
std::shared_ptr<SerializableBase> HeightFieldGeomCell::copy() const
{
  return std::shared_ptr<SerializableBase>(new HeightFieldGeomCell(*this));
}
//helper
Vec3i HeightFieldGeomCell::getStride(bool cell) const
{
  if(_dim == 2)
    return Vec3i(1,0,0);
  else return Vec3i(cell ? _h.getNrCell()[1] : _h.getNrPoint()[1],1,0);
}
void HeightFieldGeomCell::getMeshInner(ObjMesh& mesh) const
{
  Vec3i id(0,0,0),stride=getStride(false);
  Vec3 unit=Vec3::Unit(_dim-1);
  mesh.setDim((int)_dim);
  //get vertices
  mesh.getV().clear();
  for(id[0]=0; id[0]<_h.getNrPoint()[0]; id[0]++)
    for(id[1]=0; id[1]<_h.getNrPoint()[1]; id[1]++)
      mesh.getV().push_back(_h.getPt(id)+unit*_h.get(id));
  //get indices
  mesh.getI().clear();
  if(_dim == 2) {
    for(id[0]=0; id[0]<_h.getNrPoint()[0]-1; id[0]++)
      mesh.getI().push_back(Vec3i(GI(1,0),GI(0,0),0));
  } else {
    for(id[0]=0; id[0]<_h.getNrPoint()[0]-1; id[0]++)
      for(id[1]=0; id[1]<_h.getNrPoint()[1]-1; id[1]++) {
        mesh.getI().push_back(Vec3i(GI(0,0),GI(1,0),GI(1,1)));
        mesh.getI().push_back(Vec3i(GI(0,0),GI(1,1),GI(0,1)));
      }
  }
}
DEVICE_ONLY_FUNC bool HeightFieldGeomCell::distInner(const Vec3& pt,Vec3& n) const
{
  if(!_bb.contain(pt,_dim))
    return false;
  Vec3 cp;
  scalar dist=ScalarUtil<scalar>::scalar_max(),minDist=dist;
#ifdef HEIGHT_FIELD_AS_MESH
  BVHQuery<sizeType>(_bvh,_dim,-1).pointDistQuery(pt,*this,cp,n,dist,&minDist);
#else
  pointDistQuery(pt,cp,n,dist,minDist);
#endif
  cp-=pt;
  if(cp.dot(n) < 1E-6f)
    return false;
  n=cp;
  return dist < ScalarUtil<scalar>::scalar_max();
}
DEVICE_ONLY_FUNC bool HeightFieldGeomCell::closestInner(const Vec3& pt,Vec3& n,Vec3* normal) const
{
  Vec3 cp,nor;
  scalar dist=ScalarUtil<scalar>::scalar_max(),minDist=dist;
  cp.block(0,0,_dim,1).setConstant(dist);

#ifdef HEIGHT_FIELD_AS_MESH
  BVHQuery<sizeType>(_bvh,_dim,-1).pointDistQuery(pt,*this,cp,nor,dist,&minDist);
#else
  pointDistQuery(pt,cp,nor,dist,minDist);
#endif
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
DEVICE_ONLY_FUNC BBox<scalar> HeightFieldGeomCell::getBBInner() const
{
  BBox<scalar> bb=_h.getBB();
  _h.minMax(bb._minC[_dim-1],bb._maxC[_dim-1]);
  return bb;
}
void HeightFieldGeomCell::build(bool buildBVHAsWell)
{
  StaticGeomCell::build(buildBVHAsWell);
  _bb=getBBInner();
  _depth=_bb.getExtent().maxCoeff();
  _bb._minC[_dim-1]-=_depth;
}
void HeightFieldGeomCell::debugWrite() const
{
  std::string path="./heightFieldVTK";
  if(std::experimental::filesystem::v1::exists(path))
    std::experimental::filesystem::v1::remove_all(path);
  std::experimental::filesystem::v1::create_directory(path);
  //debug mesh
  ObjMesh mesh;
  getMesh(mesh,true);
  mesh.writeVTK(path+"/mesh.vtk",true);
  //debug BVH
  writeBVHByLevel<sizeType,BBox<scalar> >(_bvh,-1,path);
  //debug ring
#define NR_TEST 1000
#define NR_RING 5
  for(sizeType i=0; i<NR_TEST; i++) {
    //find point for test
    Vec3 pos;
    while(true) {
      pos=Vec3::Random()*getBB().getExtent().maxCoeff();
      pos.segment(_dim,3-_dim).setZero();
      if(getBBInner().contain(pos,_dim))
        break;
    }
    //write bounding box
    std::stack<sizeType> ss;
    std::vector<scalar> css;
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
    sizeType nrCell=_h.getNrCell().segment(0,_dim-1).prod();
    for(sizeType j=0; j<NR_RING; j++) {
      addRing(ss,floorV(_h.getIndexFrac(pos)),j);
      while(!ss.empty()) {
        sizeType bvhId=ss.top();
        for(sizeType b=0; b<_dim-1; b++) {
          ASSERT(bvhId+b >= 0 && bvhId+b < nrCell*2)
          const BBox<scalar>& bb=_bvh[bvhId+b]._bb;
          vss.push_back(bb._minC);
          vss.push_back(bb._maxC);
          css.push_back((scalar)j);
        }
        ss.pop();
      }
    }
    VTKWriter<scalar> os("Ring",path+"/ring"+std::to_string(i)+".vtk",true);
    os.appendVoxels(vss.begin(),vss.end(),true);
    os.appendCustomData("ring",css.begin(),css.end());
  }
#undef NR_RING
#undef NR_TEST
}
Mat4 HeightFieldGeomCell::genT(sizeType dim,sizeType dimH) const
{
  Mat4 ret=Mat4::Identity();
  typedef ScalarUtil<scalar>::ScalarQuat QUAT;
  ret.block<3,3>(0,0)=QUAT::FromTwoVectors(Vec3::Unit(dim),Vec3::Unit(dimH)).toRotationMatrix();
  return ret;
}
void HeightFieldGeomCell::pointDistQuery(Vec3 pt,Vec3& cp,Vec3& n,scalar& dist,scalar& minDist) const
{
  sizeType ring=0;
  std::stack<sizeType> ss;
  pt.segment(_dim,3-_dim).setZero();
  Vec3i id=floorV(_h.getIndexFrac(pt));
  addRing(ss,id,ring++);
  while(true) {
    //check current ring
    bool more=false;
    while(!ss.empty()) {
      sizeType bvhId=ss.top();
      ss.pop();

      const BBox<scalar>& bb=_bvh[bvhId]._bb;
      if(bb.distTo(pt,_dim-1) > std::min<scalar>(depth(),dist))
        continue;

      Vec2i feat;
      more=true;
      if(_dim == 3) {
        calcMinDist3D(_iss[bvhId+0],pt,cp,n,dist,&minDist,feat);
        calcMinDist3D(_iss[bvhId+1],pt,cp,n,dist,&minDist,feat);
      } else {
        calcMinDist2D(_iss[bvhId+0],pt,cp,n,dist,&minDist,feat);
      }
    }
    if(!more)
      break;
    //add new ring
    addRing(ss,id,ring++);
  }
}
void HeightFieldGeomCell::addRing(std::stack<sizeType>& ss,const Vec3i& id,sizeType ring) const
{
  sizeType dim=_dim-1;
  Vec3i idr=id,nrC=_h.getNrCell(),stride=getStride(true);
  for(sizeType i=0; i<dim; i++)
    idr[i]=std::max<sizeType>(std::min<sizeType>(idr[i],nrC[i]-1),0);

  if(ring == 0)
    ss.push(stride.dot(idr)*dim);
  else if(_dim == 2) {
    //left
    idr[0]=id[0]-ring;
    if(idr[0] >= 0)
      ss.push(stride.dot(idr)*dim);
    //right
    idr[0]=id[0]+ring;
    if(idr[0] < nrC[0])
      ss.push(stride.dot(idr)*dim);
  } else {
    //left
    idr[0]=id[0]-ring;
    if(idr[0] >= 0)
      for(idr[1] =std::max<sizeType>(id[1]-ring,0);
          idr[1]<=std::min<sizeType>(id[1]+ring,nrC[1]-1); idr[1]++)
        ss.push(stride.dot(idr)*dim);
    //right
    idr[0]=id[0]+ring;
    if(idr[0] < nrC[0])
      for(idr[1] =std::max<sizeType>(id[1]-ring,0);
          idr[1]<=std::min<sizeType>(id[1]+ring,nrC[1]-1); idr[1]++)
        ss.push(stride.dot(idr)*dim);
    //bottom
    idr[1]=id[1]-ring;
    if(idr[1] >= 0)
      for(idr[0] =std::max<sizeType>(id[0]-ring+1,0);
          idr[0]<=std::min<sizeType>(id[0]+ring-1,nrC[0]-1); idr[0]++)
        ss.push(stride.dot(idr)*dim);
    //top
    idr[1]=id[1]+ring;
    if(idr[1] < nrC[1])
      for(idr[0] =std::max<sizeType>(id[0]-ring+1,0);
          idr[0]<=std::min<sizeType>(id[0]+ring-1,nrC[0]-1); idr[0]++)
        ss.push(stride.dot(idr)*dim);
  }
}
#undef GI
