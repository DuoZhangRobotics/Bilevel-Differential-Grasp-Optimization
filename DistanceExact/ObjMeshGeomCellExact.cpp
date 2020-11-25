#include "ObjMeshGeomCellExact.h"
#include "Utils/DebugGradient.h"
#include "CommonFile/IO.h"
#include "MPQZIO.h"
#include <stack>

USE_PRJ_NAMESPACE

ObjMeshGeomCellExact::ObjMeshGeomCellExact() {}
ObjMeshGeomCellExact::ObjMeshGeomCellExact(const ObjMeshGeomCell& exact)
{
  //vss
  _vss.resize(exact.vss().size());
  for(sizeType i=0; i<(sizeType)_vss.size(); i++)
    _vss[i]=exact.vss()[i].cast<T>();
  //iss
  _iss=exact.iss();
  //tss
  _tss.resize(_iss.size());
  for(sizeType i=0; i<(sizeType)_tss.size(); i++)
    _tss[i]=TriangleExact(_vss[_iss[i][0]],_vss[_iss[i][1]],_vss[_iss[i][2]]);
  //bvh
  _bvh.resize(exact.bvh().size());
  for(sizeType i=0; i<(sizeType)_bvh.size(); i++) {
    Node<sizeType,BBoxExact>& n=_bvh[i];
    const Node<sizeType>& nRef=exact.bvh()[i];
    n._l=nRef._l;
    n._r=nRef._r;
    n._parent=nRef._parent;
    n._nrCell=nRef._nrCell;
    n._cell=nRef._cell;
    n._bb=BBoxExact(nRef._bb);
  }
}
bool ObjMeshGeomCellExact::read(std::istream& is,IOData* dat)
{
  readBinaryData(_vss,is,dat);
  readBinaryData(_iss,is,dat);
  readBinaryData(_tss,is,dat);
  readBinaryData(_bvh,is,dat);
  return is.good();
}
bool ObjMeshGeomCellExact::write(std::ostream& os,IOData* dat) const
{
  writeBinaryData(_vss,os,dat);
  writeBinaryData(_iss,os,dat);
  writeBinaryData(_tss,os,dat);
  writeBinaryData(_bvh,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> ObjMeshGeomCellExact::copy() const
{
  return std::shared_ptr<SerializableBase>(new ObjMeshGeomCellExact);
}
std::string ObjMeshGeomCellExact::type() const
{
  return typeid(ObjMeshGeomCellExact).name();
}
const BBoxExact& ObjMeshGeomCellExact::getBB() const
{
  return _bvh.back()._bb;
}
bool ObjMeshGeomCellExact::closest(const PT& pt,PT& n,PT& normal,MAT3& hessian,Vec2i& feat) const
{
  Vec2i featTmp;
  PT cp,cpTmp,bTmp;
  feat=Vec2i(-1,-1);
  T distSqr=ScalarUtil<double>::scalar_max(),distSqrTmp;
  cp.setConstant(distSqr);
  //main loop
  std::stack<std::pair<T,sizeType>> ss;
  ss.push(std::make_pair(_bvh.back()._bb.distToSqr(pt),(sizeType)_bvh.size()-1));
  while(!ss.empty()) {
    T distSqrToBB=ss.top().first;
    const Node<sizeType,BBoxExact>& node=_bvh[ss.top().second];
    ss.pop();
    if(distSqrToBB>distSqr)
      continue;
    else if(node._cell>=0) {
      const TriangleExact& te=_tss[node._cell];
      te.calcPointDist(pt,distSqrTmp,cpTmp,bTmp,featTmp);
      //get feature id
      for(int d=0; d<2; d++)
        if(featTmp[d]>=0)
          featTmp[d]=_iss[node._cell][featTmp[d]];
      //make feature id Unique
      if(featTmp[1]>=0 && featTmp[0]>featTmp[1])
        std::swap(featTmp[0],featTmp[1]);
      //update distance
      if(featTmp[0]>=0 && feat==featTmp) {
        //update due to same feature
        ASSERT(distSqr==distSqrTmp)
        T align=abs((cp-pt).dot(normal));
        T alignCurr=abs((cpTmp-pt).dot(te.normal()));
        if(alignCurr>align) {
          distSqr=distSqrTmp;
          normal=te.normal();
          cp=cpTmp;
          n=cp-pt;
        }
      } else if(distSqrTmp<distSqr) {
        //update due to new distance
        feat=featTmp;
        distSqr=distSqrTmp;
        normal=te.normal();
        cp=cpTmp;
        n=cp-pt;
      }
    } else {
      const Node<sizeType,BBoxExact>& nl=_bvh[node._l];
      const Node<sizeType,BBoxExact>& nr=_bvh[node._r];
      T distSqrToBBL=nl._bb.distToSqr(pt);
      T distSqrToBBR=nr._bb.distToSqr(pt);
      if(distSqrToBBL<distSqrToBBR) {
        ss.push(std::make_pair(distSqrToBBR,node._r));
        ss.push(std::make_pair(distSqrToBBL,node._l));
      } else {
        ss.push(std::make_pair(distSqrToBBL,node._l));
        ss.push(std::make_pair(distSqrToBBR,node._r));
      }
    }
  }
  //adjust hessian
  if(distSqr==0) {
    //surface
    hessian.setZero();
  } else if(feat[0]==-1) {
    //surface
    hessian.setZero();
  } else if(feat[1]==-1) {
    //vertex
    ASSERT(feat[0]>=0)
    hessian=n*n.transpose();
    hessian/=hessian.trace();
    hessian-=MAT3::Identity();
  } else {
    //edge
    ASSERT(feat[0]>=0 && feat[1]>=0)
    PT e=_vss[feat[0]]-_vss[feat[1]];
    T eDotE=e.dot(e);
    PT d=n-n.dot(e)*e/eDotE;
    T dDotD=d.dot(d);
    PT nd=e.dot(d)*e/eDotE/dDotD;
    hessian=e*e.transpose()/eDotE-MAT3::Identity();
    hessian+=d*d.transpose()/dDotD;
    hessian-=d*nd.transpose();
  }
  return n.dot(normal)>0;
}
scalar ObjMeshGeomCellExact::closest(const Vec3& pt,Vec3& n,Vec3& normal,Mat3& hessian,Vec2i& feat) const
{
  bool inside;
  MAT3 hessianR;
  PT ptR=pt.cast<T>(),nR,normalR;
  inside=closest(ptR,nR,normalR,hessianR,feat);
  //cast
  n=castRational<Vec3,PT>(nR);
  normal=castRational<Vec3,PT>(normalR);
  hessian=castRational<Mat3,MAT3>(hessianR);
  //post process
  scalar nLen=n.norm();
  hessian/=std::max(std::numeric_limits<scalar>::epsilon(),nLen);
  if(nR[0]==0 && nR[1]==0 && nR[2]==0) {
    nLen=0;
    normal/=std::max(std::numeric_limits<scalar>::epsilon(),normal.norm());
  } else {
    normal=n/std::max(std::numeric_limits<scalar>::epsilon(),nLen);
    if(inside)
      nLen*=-1;
    else {
      normal*=-1;
      hessian*=-1;
    }
  }
  return nLen;
}
void ObjMeshGeomCellExact::writePointDistVTK(const std::string& path,sizeType res) const
{
  VTKWriter<scalar> os("PointDist",path,true);
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> iss;
  std::vector<scalar> css;
  Vec3 p,n,normal,alpha;
  Mat3 hessian;
  Vec2i feat;
  BBoxExact bb=getBB();
  BBox<scalar> bbs(castRational<Vec3>(bb._minC),castRational<Vec3>(bb._maxC));
  bbs.enlargedEps(1);
  for(alpha[0]=0; alpha[0]<=res; alpha[0]++)
    for(alpha[1]=0; alpha[1]<=res; alpha[1]++)
      for(alpha[2]=0; alpha[2]<=res; alpha[2]++) {
        for(sizeType d=0; d<3; d++)
          p[d]=bbs._minC[d]*(1-alpha[d]/res)+bbs._maxC[d]*alpha[d]/res;
        bool in=closest(p,n,normal,hessian,feat)<0;
        vss.push_back(p);
        vss.push_back(n+p);
        iss.push_back(Vec2i(vss.size()-2,vss.size()-1));
        css.push_back(in?1:0);
      }
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<scalar>::LINE);
  os.appendCustomData("InsideOutside",css.begin(),css.end());
}
void ObjMeshGeomCellExact::debugPointDist(sizeType nrIter) const
{
  DEFINE_NUMERIC_DELTA_T(scalar)
  Vec2i feat,feat2;
  Vec3 n,normal,n2,normal2;
  Mat3 hessian,hessian2;
  BBox<scalar> bb;
  bb._minC=castRational<Vec3,PT>(getBB()._minC);
  bb._maxC=castRational<Vec3,PT>(getBB()._maxC);
  for(sizeType i=0; i<nrIter; i++) {
    Vec3 pt=Vec3::Random()*0.5+Vec3::Constant(0.5);
    pt.array()=bb._minC.array()*(1-pt.array())+bb._maxC.array()*pt.array();
    scalar dist=closest(pt,n,normal,hessian,feat);
    std::string type=dist<0?"I":"O";
    type+=feat[0]<0?"T":feat[1]<0?"V":"E";
    for(sizeType d=0; d<3; d++) {
      std::string name="N"+type+std::to_string(d);
      scalar dist2=closest(pt+Vec3::Unit(d)*DELTA,n2,normal2,hessian2,feat2);
      DEBUG_GRADIENT(name.c_str(),normal[d],normal[d]-(dist2-dist)/DELTA)
    }
    if(!hessian.isZero())
      for(sizeType d=0; d<3; d++) {
        std::string name="H"+type+std::to_string(d);
        closest(pt+Vec3::Unit(d)*DELTA,n2,normal2,hessian2,feat2);
        DEBUG_GRADIENT(name.c_str(),(hessian*Vec3::Unit(d)).norm(),(hessian*Vec3::Unit(d)-(normal2-normal)/DELTA).norm())
      }
  }
}
sizeType ObjMeshGeomCellExact::nrV() const
{
  return (sizeType)_vss.size();
}
sizeType ObjMeshGeomCellExact::nrI() const
{
  return (sizeType)_iss.size();
}
