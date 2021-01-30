#include "ConvexHullExact.h"
#include "ObjMeshGeomCellExact.h"
#include <CommonFile/CameraModel.h>
#include <Utils/DebugGradient.h>
#include <Utils/Utils.h>
#include <Utils/Hash.h>
#include <CommonFile/IO.h>
#include <iomanip>
#include <stack>

USE_PRJ_NAMESPACE

ConvexHullExact::ConvexHullExact()
{
  buildBD();
}
ConvexHullExact::ConvexHullExact(const ObjMeshGeomCell& exact)
{
  //vss
  BBoxExact bb;
  _vss.resize(exact.vss().size());
  for(sizeType i=0; i<(sizeType)_vss.size(); i++) {
    _vss[i]=(exact.getT().block<3,3>(0,0)*exact.vss()[i]+exact.getT().block<3,1>(0,3)).cast<double>().cast<T>();
    bb.setUnion(_vss[i]);
  }
  if(!_vss.empty()) {

    _bvh.resize(1);
    _bvh[0]._bb=bb;
  }
  //iss
  _iss=exact.iss();
  //tss
  _tss.resize(_iss.size());

  for(sizeType i=0; i<(sizeType)_tss.size(); i++) {
    _tss[i]=TriangleExact(_vss[_iss[i][0]],_vss[_iss[i][1]],_vss[_iss[i][2]]);
    _tss[i]._vid=_iss[i];
  }

  //ess
  std::unordered_map<Vec2i,EdgeExact,Hash> ess;
  for(sizeType i=0; i<(sizeType)_tss.size(); i++)
    for(sizeType d=0; d<3; d++) {
      Vec2i id(_iss[i][d],_iss[i][(d+1)%3]);
      if(id[0]>id[1])
        std::swap(id[0],id[1]);
      if(ess.find(id)==ess.end()) {
        EdgeExact E(_vss[_iss[i][d]],_vss[_iss[i][(d+1)%3]]);
        E._vid=Vec2i(_iss[i][d],_iss[i][(d+1)%3]);
        E._tNId=Vec2i(i,-1);
        ess[id]=E;
      } else
        ess.find(id)->second._tNId[1]=i;
    }

  _eNss.resize(_vss.size());
  for(const std::pair<Vec2i,EdgeExact>& E:ess) {
    for(sizeType d=0; d<2; d++) {
      sizeType tid=E.second._tNId[d];
      TriangleExact& T=_tss[tid];
      bool found=false;
      for(sizeType d2=0; d2<3; d2++)
        if(T._vid[d2]!=E.second._vid[0] && T._vid[d2]!=E.second._vid[1]) {
          T._eNId[d2]=(sizeType)_ess.size();
          found=true;
          break;
        }
      ASSERT_MSGV(found,"Cannot find Edge (%d,%d) in Triangle %d",E.first[0],E.first[1],tid)
    }
    _eNss[E.second._vid[0]].push_back((sizeType)_ess.size());
    _eNss[E.second._vid[1]].push_back((sizeType)_ess.size());
    _ess.push_back(E.second);

  }
  parityCheck();
  buildBD();

}
ConvexHullExact::~ConvexHullExact()
{
  removeBD();
}
bool ConvexHullExact::read(std::istream& is,IOData* dat)
{
  ObjMeshGeomCellExact::read(is,dat);
  readBinaryData(_ess,is,dat);
  readBinaryData(_eNss,is,dat);
  removeBD();
  buildBD();
  return is.good();
}
bool ConvexHullExact::write(std::ostream& os,IOData* dat) const
{
  ObjMeshGeomCellExact::write(os,dat);
  writeBinaryData(_ess,os,dat);
  writeBinaryData(_eNss,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> ConvexHullExact::copy() const
{
  return std::shared_ptr<SerializableBase>(new ConvexHullExact);
}
std::string ConvexHullExact::type() const
{
  return typeid(ConvexHullExact).name();
}
void ConvexHullExact::parityCheck() const
{
  for(sizeType tid=0; tid<(sizeType)_tss.size(); tid++)
    for(sizeType d=0; d<3; d++) {
      sizeType eid=_tss[tid]._eNId[d];
      ASSERT_MSGV(eid>=0,"Triangle %d's Edge %d is -1",tid,d)
      const EdgeExact& e=_ess[eid];
      ASSERT_MSGV(e._tNId[0]==tid || e._tNId[1]==tid,"Edge %d does not contain Traingle %d",eid,tid)
    }
  for(sizeType eid=0; eid<(sizeType)_ess.size(); eid++)
    for(sizeType d=0; d<2; d++) {
      sizeType tid=_ess[eid]._tNId[d];
      ASSERT_MSGV(tid>=0,"Edge %d's Triangle %d is -1",eid,d)
      const TriangleExact& t=_tss[tid];
      ASSERT_MSGV(t._eNId[0]==eid || t._eNId[1]==eid || t._eNId[2]==eid,"Triangle %d does not contain Edge %d",tid,eid)
      sizeType vid=_ess[eid]._vid[d];
      ASSERT_MSGV(std::find(_eNss[vid].begin(),_eNss[vid].end(),eid)!=_eNss[vid].end(),"Cannot find Edge %d as neighbor of Vertex %d",eid,vid)
    }
  for(sizeType vid=0; vid<(sizeType)_eNss.size(); vid++)
    for(sizeType eid:_eNss[vid]) {
      const EdgeExact& e=_ess[eid];
      ASSERT_MSGV(e._vid[0]==vid || e._vid[1]==vid,"Edge %d does not contain Vertex %d",eid,vid)
    }
}
bool ConvexHullExact::closest(const PT& pt,PT& n,PT& normal,MAT3& hessian,Vec2i& feat,bool cache,std::vector<PT,Eigen::aligned_allocator<PT>>* history) const
{
  //initialize
  sizeType bestVid=-1;
  T distSqr,bestDistSqr;
  if(cache) {
    if(feat[0]>=0)
      bestVid=feat[0];
    else bestVid=_tss[feat[1]]._vid[0];
  } else {
    for(sizeType vid=0; vid<(sizeType)_vss.size(); vid++) {
      distSqr=(_vss[vid]-pt).squaredNorm();
      if(bestVid==-1 || distSqr<bestDistSqr) {
        bestDistSqr=distSqr;
        bestVid=vid;
      }
    }
  }
  //local optimization
  PT p=_vss[bestVid];
  if(history)
    history->push_back(p);
  sizeType currentVid=bestVid;
  sizeType currentEid=-1;
  sizeType currentTid=-1;
  while(true)
    if(currentTid>=0) {
      ASSERT(currentVid==-1 && currentEid>=0)
      const EdgeExact& E=_ess[currentEid];
      const TriangleExact& T=_tss[currentTid];
      std::pair<sizeType,PT> pT=T.moveFromEdge(E._vid,p,pt-p);
      p=pT.second;
      if(history)
        history->push_back(p);
      if(pT.first==-1) {
        normal=T.normal();
        feat=Vec2i(-1,currentTid);
        break;
      } else if(T._vid[pT.first]!=E._vid[0] && T._vid[pT.first]!=E._vid[1]) {
        currentVid=T._vid[pT.first];
        currentEid=-1;
        currentTid=-1;
      } else {
        currentVid=-1;
        currentEid=T._eNId[pT.first];
        currentTid=-1;
      }
    } else if(currentVid>=0) {
      ASSERT(currentEid==-1 && currentTid==-1)
      //check edge
      T bestDirGrad=0,dirGrad;
      for(sizeType eNId:_eNss[currentVid]) {
        dirGrad=_ess[eNId].dirGrad(currentVid,pt-p);
        if(dirGrad>bestDirGrad) {
          currentEid=eNId;
          bestDirGrad=dirGrad;
        }
      }
      if(bestDirGrad>0)
        currentVid=-1;
      //no edge, return
      if(currentVid>=0) {
        sizeType eNId=_eNss[currentVid][0];
        normal=_tss[_ess[eNId]._tNId[0]].normal();
        feat=Vec2i(currentVid,-1);
        break;
      }
    } else if(currentEid>=0) {
      ASSERT(currentVid==-1 && currentTid==-1)
      const EdgeExact& E=_ess[currentEid];
      std::pair<sizeType,PT> pE=E.moveFromVertex(p,pt-p);
      p=pE.second;
      if(history)
        history->push_back(p);
      if(pE.first>=0) {
        currentVid=pE.first;
        currentEid=-1;
      } else {
        //check triangle
        for(sizeType d=0; d<2; d++) {
          const TriangleExact& T=_tss[E._tNId[d]];
          if(T.dirGrad(E._vid,pt-p)>0) {
            currentTid=E._tNId[d];
            break;
          }
        }
        //no triangle, return
        if(currentTid==-1) {
          normal=_tss[E._tNId[0]].normal();
          feat=E._vid;
          break;
        }
      }
    } else {
      ASSERT_MSG(false,"Invalid configuration with vid=eid=tid=-1")
    }
  n=p-pt;
  distSqr=n.squaredNorm();
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
void ConvexHullExact::writeHistoryDistVTK(const std::string& path,sizeType res) const
{
  recreate(path);
  Vec3 p,n,normal,alpha;
  Mat3 hessian;
  Vec2i feat;
  sizeType id=0;
  BBoxExact bb=getBB();
  BBox<scalar> bbs(castRational<Vec3>(bb._minC),castRational<Vec3>(bb._maxC));
  bbs.enlargedEps(1);
  for(alpha[0]=0; alpha[0]<=res; alpha[0]++)
    for(alpha[1]=0; alpha[1]<=res; alpha[1]++)
      for(alpha[2]=0; alpha[2]<=res; alpha[2]++) {
        for(sizeType d=0; d<3; d++)
          p[d]=bbs._minC[d]*(1-alpha[d]/res)+bbs._maxC[d]*alpha[d]/res;
        std::vector<Vec3,Eigen::aligned_allocator<Vec3>> history;
        feat=Vec2i(RandEngine::randI(0,(sizeType)_vss.size()-1),-1);
        bool in=ObjMeshGeomCellExact::closest<scalar>(p,n,normal,hessian,feat,true,&history)<0;
        if(!in) {
          {
            VTKWriter<scalar> os("PointDist",path+"/frm"+std::to_string(id)+".vtk",true);
            std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
            std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> iss;
            std::vector<scalar> css;
            vss.push_back(p);
            vss.push_back(n+p);
            iss.push_back(Vec2i(vss.size()-2,vss.size()-1));
            css.push_back(in?1:0);
            os.appendPoints(vss.begin(),vss.end());
            os.appendCells(iss.begin(),iss.end(),VTKWriter<scalar>::LINE);
            os.appendCustomData("InsideOutside",css.begin(),css.end());
          }
          {
            VTKWriter<scalar> os("PointDist",path+"/history"+std::to_string(id)+".vtk",true);
            std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> iss;
            os.appendPoints(history.begin(),history.end());
            os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,0,1),
                           VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)history.size()-1,0,1),
                           VTKWriter<scalar>::LINE);
          }
          id++;
        }
      }
  getMesh().writeVTK(path+"/mesh.vtk",true);
}
void ConvexHullExact::compare(const ObjMeshGeomCellExact& mExact,sizeType res) const
{
  Vec3 p,alpha;
  Vec3d n3,normal3;
  PT n,n2,normal,normal2;
  MAT3 hessian,hessian2;
  Vec2i feat,feat2;
  BBoxExact bb=getBB();
  BBox<scalar> bbs(castRational<Vec3>(bb._minC),castRational<Vec3>(bb._maxC));
  bbs.enlargedEps(1);
  for(alpha[0]=0; alpha[0]<=res; alpha[0]++)
    for(alpha[1]=0; alpha[1]<=res; alpha[1]++)
      for(alpha[2]=0; alpha[2]<=res; alpha[2]++) {
        for(sizeType d=0; d<3; d++)
          p[d]=bbs._minC[d]*(1-alpha[d]/res)+bbs._maxC[d]*alpha[d]/res;
        PT P=p.cast<scalarD>().cast<T>();
        bool in=closest(P,n,normal,hessian,feat);
        if(!in) {
          bool in2=mExact.closest(P,n2,normal2,hessian2,feat2);
          ASSERT(in==in2)
          if(feat[0]==-1) {
            PT b=_tss[feat[1]].bary(P+n);
            PT nRef=_tss[feat[1]].interp(b)-P;
            ASSERT(nRef.cross(normal).isZero())
            ASSERT((nRef-n).isZero())
            Vec3d bd=castRational<Vec3d,PT>(b);
            std::cout << "BarycentricCoordinate-ConvexHullExact(" << feat2[0] << "," << feat2[1] << "): " << bd.transpose() << std::endl;
          }
          if(feat2[0]==-1) {
            PT b=_tss[feat2[1]].bary(P+n2);
            PT nRef=_tss[feat2[1]].interp(b)-P;
            ASSERT(nRef.cross(normal2).isZero())
            ASSERT((nRef-n2).isZero())
            Vec3d bd2=castRational<Vec3d,PT>(b);
            std::cout << "BarycentricCoordinate-ObjMeshGeomCellExact(" << feat2[0] << "," << feat2[1] << "): " << bd2.transpose() << std::endl;
          }
          {
            closestGJK<scalarD>(p.cast<scalarD>(),n3,normal3);
            Vec3 nd=castRational<Vec3,PT>(n);
            std::cout << "n(" << nd[0] << "," << nd[1] << "," << nd[2] << ") ";
            std::cout << "nGJK(" << n3[0] << "," << n3[1] << "," << n3[2] << ")" << std::endl;
          }
        }
      }
}
void ConvexHullExact::buildBD()
{
  typedef double* doublePtr;
  _bd.numpoints=(int)_vss.size();
  if(_bd.numpoints>0) {
    _bd.coord=new doublePtr[_bd.numpoints];
    for(sizeType i=0; i<_bd.numpoints; i++) {
      _bd.coord[i]=new double[3];
      for(sizeType d=0; d<3; d++)
        castRational(_vss[i][d],_bd.coord[i][d]);
    }
  } else {
    _bd.coord=NULL;
  }
  _bd.s[0]=_bd.s[1]=_bd.s[2]=0;
}
void ConvexHullExact::removeBD()
{
  if(_bd.numpoints>0) {
    for(sizeType i=0; i<_bd.numpoints; i++)
      delete [] _bd.coord[i];
    delete [] _bd.coord;
  }
  _bd.numpoints=0;
  _bd.coord=NULL;
}
