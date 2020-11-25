#include "MakeMesh.h"
#include "Interp.h"
#include <map>

USE_PRJ_NAMESPACE

std::map<sizeType,ObjMesh> sphere3DCache;
//MakeMesh
void MakeMesh::makeTet3D(ObjMesh& m)
{
  m.getV().clear();

  m.getV().push_back(Vec3(0.0f,0.0f,0.0f));
  m.getV().push_back(Vec3(1.0f,0.0f,0.0f));
  m.getV().push_back(Vec3(0.0f,1.0f,0.0f));
  m.getV().push_back(Vec3(0.0f,0.0f,1.0f));

  m.getI().push_back(Vec3i(0,2,1));
  m.getI().push_back(Vec3i(0,1,3));
  m.getI().push_back(Vec3i(0,3,2));
  m.getI().push_back(Vec3i(1,2,3));
  m.getDim()=3;
}
void MakeMesh::makeBox3D(ObjMesh& m,const Vec3& ext)
{
  m.getV().clear();
#define ADDV(S1,S2,S3) m.getV().push_back(Vec3(S1*ext[0],S2*ext[1],S3*ext[2]));
  ADDV(-1.0f,-1.0f,-1.0f)
  ADDV(+1.0f,-1.0f,-1.0f)
  ADDV(+1.0f,+1.0f,-1.0f)
  ADDV(-1.0f,+1.0f,-1.0f)
  ADDV(-1.0f,-1.0f,+1.0f)
  ADDV(+1.0f,-1.0f,+1.0f)
  ADDV(+1.0f,+1.0f,+1.0f)
  ADDV(-1.0f,+1.0f,+1.0f)
  m.getI().clear();
#define ADDI(A,B,C) m.getI().push_back(Vec3i(A,B,C));
  //-Z
  ADDI(0,2,1)
  ADDI(0,3,2)
  //+Z
  ADDI(4,5,6)
  ADDI(4,6,7)
  //-X
  ADDI(0,4,7)
  ADDI(0,7,3)
  //+X
  ADDI(1,2,6)
  ADDI(1,6,5)
  //-Y
  ADDI(0,1,5)
  ADDI(0,5,4)
  //+Y
  ADDI(2,3,7)
  ADDI(2,7,6)
  m.applyTrans();
  m.makeUniform();
  m.getDim()=3;
}
void MakeMesh::makeBox3D(ObjMesh& m,const Vec3& ext,scalar thick)
{
  ObjMesh tmp;
  makeBox3D(tmp,ext);
  tmp.insideOut();
  m.addMesh(tmp,"inner");
  makeBox3D(tmp,ext+Vec3(thick,thick,thick));
  m.addMesh(tmp,"outer");
  m.getDim()=3;
}
void MakeMesh::makeDiscreteBox3D(ObjMesh& m,const Vec3& ext)
{
  m.getV().clear();
  m.getI().clear();
  sizeType id;
#define ADDV(S1,S2,S3) m.getV().push_back(Vec3(S1*ext[0],S2*ext[1],S3*ext[2]));
  //top
  id=(sizeType)m.getV().size();
  ADDV(-1.0f,-1.0f,+1.0f)
  ADDV(+1.0f,-1.0f,+1.0f)
  ADDV(+1.0f,+1.0f,+1.0f)
  ADDV(-1.0f,+1.0f,+1.0f)
  ADDI(id+0,id+1,id+2)
  ADDI(id+0,id+2,id+3)
  //bottom
  id=(sizeType)m.getV().size();
  ADDV(-1.0f,-1.0f,-1.0f)
  ADDV(-1.0f,+1.0f,-1.0f)
  ADDV(+1.0f,+1.0f,-1.0f)
  ADDV(+1.0f,-1.0f,-1.0f)
  ADDI(id+0,id+1,id+2)
  ADDI(id+0,id+2,id+3)
  //right
  id=(sizeType)m.getV().size();
  ADDV(+1.0f,-1.0f,-1.0f)
  ADDV(+1.0f,+1.0f,-1.0f)
  ADDV(+1.0f,+1.0f,+1.0f)
  ADDV(+1.0f,-1.0f,+1.0f)
  ADDI(id+0,id+1,id+2)
  ADDI(id+0,id+2,id+3)
  //left
  id=(sizeType)m.getV().size();
  ADDV(-1.0f,-1.0f,-1.0f)
  ADDV(-1.0f,-1.0f,+1.0f)
  ADDV(-1.0f,+1.0f,+1.0f)
  ADDV(-1.0f,+1.0f,-1.0f)
  ADDI(id+0,id+1,id+2)
  ADDI(id+0,id+2,id+3)
  //front
  id=(sizeType)m.getV().size();
  ADDV(-1.0f,+1.0f,-1.0f)
  ADDV(-1.0f,+1.0f,+1.0f)
  ADDV(+1.0f,+1.0f,+1.0f)
  ADDV(+1.0f,+1.0f,-1.0f)
  ADDI(id+0,id+1,id+2)
  ADDI(id+0,id+2,id+3)
  //back
  id=(sizeType)m.getV().size();
  ADDV(-1.0f,-1.0f,-1.0f)
  ADDV(+1.0f,-1.0f,-1.0f)
  ADDV(+1.0f,-1.0f,+1.0f)
  ADDV(-1.0f,-1.0f,+1.0f)
  ADDI(id+0,id+1,id+2)
  ADDI(id+0,id+2,id+3)
  m.applyTrans();
  m.getDim()=3;
}
void MakeMesh::makeBox2D(ObjMesh& m,const Vec3& ext)
{
  m.getV().clear();
  ADDV(-1.0f,-1.0f,0.0f)
  ADDV(+1.0f,-1.0f,0.0f)
  ADDV(+1.0f,+1.0f,0.0f)
  ADDV(-1.0f,+1.0f,0.0f)
  m.getI().clear();
  ADDI(0,1,0)
  ADDI(1,2,0)
  ADDI(2,3,0)
  ADDI(3,0,0)
  m.getDim()=2;
}
void MakeMesh::makeBox2D(ObjMesh& m,const Vec3& ext,scalar thick)
{
  ObjMesh tmp;
  makeBox2D(tmp,ext);
  tmp.insideOut();
  m.addMesh(tmp,"inner");
  makeBox2D(tmp,ext+Vec3(thick,thick,0.0f));
  m.addMesh(tmp,"outer");
  m.getDim()=2;
}
void MakeMesh::makeSphere3D(ObjMesh& m,const scalar rad,const sizeType slice)
{
  //put cache
  if(sphere3DCache.find(slice) == sphere3DCache.end()) {
    makeBox3D(m,Vec3::Ones());
    sizeType nr=std::convert<sizeType>()(ceil(std::log(scalar(slice))));
    m.subdivide((int)nr);
    for(sizeType i=0; i<(sizeType)m.getV().size(); i++)
      m.getV()[i]=m.getV()[i].normalized();
    m.applyTrans();
    m.getDim()=3;
    sphere3DCache[slice]=m;
  }
  //build
  m=sphere3DCache[slice];
  m.getScale()=rad;
  m.applyTrans();
}
void MakeMesh::makeSphere2D(ObjMesh& m,const scalar rad,const sizeType slice)
{
  m.getV().clear();
  for(sizeType i=0; i<slice; i++) {
    scalar ang=2.0f*M_PI*((scalar(i)+0.5f)/scalar(slice));
    m.getV().push_back(Vec3(cos(ang)*rad,sin(ang)*rad,0.0f));
  }
  m.getI().clear();
  for(sizeType i=0; i<slice; i++)
    m.getI().push_back(Vec3i((i+1)%slice,i,0));
  m.getDim()=2;
}
void MakeMesh::makeSphere3D(ObjMesh& m,const scalar rad,const sizeType slice,scalar thick)
{
  ObjMesh tmp;
  makeSphere3D(tmp,rad,slice);
  tmp.insideOut();
  m.addMesh(tmp,"inner");
  makeSphere3D(tmp,rad+thick,slice);
  m.addMesh(tmp,"outer");
  m.getDim()=3;
}
void MakeMesh::makeSphere2D(ObjMesh& m,const scalar rad,const sizeType slice,scalar thick)
{
  ObjMesh tmp;
  makeSphere2D(tmp,rad,slice);
  tmp.insideOut();
  m.addMesh(tmp,"inner");
  makeSphere2D(tmp,rad+thick,slice);
  m.addMesh(tmp,"outer");
  m.getDim()=2;
}
//capsule3D
void separateDim(ObjMesh& m,sizeType d,scalar offDim,sizeType N)
{
  ObjMesh::EdgeMap edge;
  std::map<int,std::vector<int> > vss;
  std::vector<std::pair<int,int> > ess;
  typedef std::map<std::pair<int,int>,ObjMesh::Edge,ObjMesh::EdgeMap::LSS> EMap;
  if(N<1)
    return;
  //build edges
  m.buildEdge(edge);
  sizeType nrV=(sizeType)m.getV().size();
  for(EMap::const_iterator beg=edge._ess.begin(),end=edge._ess.end(); beg!=end; beg++) {
    const std::pair<int,int>& E=beg->first;
    if(std::abs(m.getV(E.first)[d])<ScalarUtil<scalar>::scalar_eps() && std::abs(m.getV(E.second)[d])<ScalarUtil<scalar>::scalar_eps()) {
      ess.push_back(E);
      if(vss.find(E.first) == vss.end()) {
        vss[E.first].push_back(E.first);
        for(sizeType i=0; i<N; i++) {
          vss[E.first].push_back((sizeType)m.getV().size());
          Vec3 off=Vec3::Unit(d)*(offDim*2*(scalarD)(i+1)/(scalarD)N-offDim);
          m.getV().push_back(m.getV()[E.first]+off);
        }
      }
      if(vss.find(E.second) == vss.end()) {
        vss[E.second].push_back(E.second);
        for(sizeType i=0; i<N; i++) {
          vss[E.second].push_back((sizeType)m.getV().size());
          Vec3 off=Vec3::Unit(d)*(offDim*2*(scalarD)(i+1)/(scalarD)N-offDim);
          m.getV().push_back(m.getV()[E.second]+off);
        }
      }
    }
  }
  //split vertices
  for(sizeType t=0; t<(sizeType)m.getI().size(); t++)
    if(m.getTC(t)[d]>0) {
      Vec3i& I=m.getI()[t];
      for(sizeType i=0; i<3; i++)
        if(vss.find(I[i]) != vss.end())
          I[i]=vss.find(I[i])->second.back();
    }
  //move vertices
  for(sizeType i=0; i<nrV; i++)
    if(m.getV(i)[d]<=ScalarUtil<scalar>::scalar_eps())
      m.getV(i)[d]-=offDim;
    else m.getV(i)[d]+=offDim;
  //insert more elements
  for(sizeType i=0; i<(sizeType)ess.size(); i++) {
    std::vector<int> I0=vss[ess[i].first];
    std::vector<int> I1=vss[ess[i].second];
    for(sizeType j=0; j<(sizeType)I0.size()-1; j++) {
      m.getI().push_back(Vec3i(I0[j+0],I0[j+1],I1[j+1]));
      m.getI().push_back(Vec3i(I0[j+0],I1[j+1],I1[j+0]));
    }
  }
  m.smooth();
}
void MakeMesh::makeCapsule3D(ObjMesh& m,const scalar rad,const Vec3& offD,const sizeType slice)
{
  makeSphere3D(m,rad,slice);
  //avgEdge
  scalarD avgEdge=0;
  ObjMesh::EdgeMap edge;
  typedef std::map<std::pair<int,int>,ObjMesh::Edge,ObjMesh::EdgeMap::LSS> EMap;
  m.buildEdge(edge);
  for(EMap::const_iterator beg=edge._ess.begin(),end=edge._ess.end(); beg!=end; beg++) {
    const std::pair<int,int>& E=beg->first;
    avgEdge+=(m.getV(E.first)-m.getV(E.second)).norm();
  }
  avgEdge/=(scalarD)edge._ess.size();
  //move
  for(sizeType d=0; d<3; d++) {
    sizeType N=(offD[d]<=0) ? 0 : (sizeType)ceil(offD[d]/avgEdge);
    separateDim(m,d,offD[d],N);
  }
  m.makeUniform();
  m.applyTrans();
  m.getDim()=3;
  if(m.getVolume()<0) {
    m.insideOut();
    m.smooth();
  }
  //m.writeVTK("mesh.vtk",true,true);
  //exit(EXIT_FAILURE);
}
void MakeMesh::makeCapsule3D(ObjMesh& m,const scalar rad,const scalar x,const scalar y,const sizeType slice)
{
  makeCapsule3D(m,rad,Vec3(x,y,0),slice);
}
void MakeMesh::makeCapsule3D(ObjMesh& m,const scalar rad,const scalar y,const sizeType slice)
{
  makeCapsule3D(m,rad,Vec3(0,y,0),slice);
}
void MakeMesh::makeCapsule3D(ObjMesh& m,const scalar rad,const scalar x,const scalar y,const scalar z,const sizeType slice)
{
  makeCapsule3D(m,rad,Vec3(x,y,z),slice);
}
//capsule2D
void MakeMesh::makeCapsule2D(ObjMesh& m,const scalar rad,const Vec2& offD,const sizeType slice)
{
  makeSphere2D(m,rad,slice);
  for(sizeType i=0; i<(sizeType)m.getV().size(); i++) {
    Vec3& v=m.getV()[i];
    for(sizeType d=0; d<2; d++)
      if(v[d] < 0)
        v-=Vec3::Unit(d)*offD[d];
      else v+=Vec3::Unit(d)*offD[d];
  }
  m.getDim()=2;
}
void MakeMesh::makeCapsule2D(ObjMesh& m,const scalar rad,const scalar y,const sizeType slice)
{
  makeCapsule2D(m,rad,Vec2(0,y),slice);
}
void MakeMesh::makeCapsule2D(ObjMesh& m,const scalar rad,const scalar x,const scalar y,const sizeType slice)
{
  makeCapsule2D(m,rad,Vec2(x,y),slice);
}
//other
void MakeMesh::makeCylinder3D(ObjMesh& m,const scalar rad,const scalar y,const sizeType slice,const sizeType sliceY,bool cap)
{
  m.getV().clear();
  m.getI().clear();
  const scalar deltaY=y*2.0f/scalar(sliceY);
  scalar Y;
  for(Y=-y; Y<=y+ScalarUtil<scalar>::scalar_eps(); Y+=deltaY) {
    sizeType IBeg=(sizeType)m.getV().size();
    for(sizeType i=0; i<slice; i++) {
      if(IBeg > 0) {
        sizeType I=IBeg+i;
        sizeType NI=IBeg+(i+1)%slice;
        m.getI().push_back(Vec3i(I,NI-slice,I-slice));
        m.getI().push_back(Vec3i(I,NI,NI-slice));
      }
      scalar ang=2.0f*M_PI*(scalar(i)/scalar(slice));
      m.getV().push_back(Vec3((scalar)cos(ang)*rad,Y,(scalar)sin(ang)*rad));
    }
  }
  if(cap) {
    sizeType base=m.getV().size()-slice;
    m.getV().push_back(Vec3(0.0f,-y,0.0f));
    m.getV().push_back(Vec3(0.0f,Y-deltaY,0.0f));
    for(sizeType i=0; i<slice; i++) {
      m.getI().push_back(Vec3i(i,(i+1)%slice,m.getV().size()-2));
      m.getI().push_back(Vec3i((i+1)%slice+base,i+base,m.getV().size()-1));
    }
  }
  m.applyTrans();
  m.getDim()=3;
}
void MakeMesh::makeTorus3D(ObjMesh& m,const scalar rad1,const scalar rad2,const sizeType slice1,const sizeType slice2)
{
  m.getV().clear();
  for(sizeType i=0; i<slice1; i++) {
    scalar angI=2.0f*M_PI*(scalar(i)/scalar(slice1));
    Vec3 ctr=Vec3(cos(angI),sin(angI),0.0f)*rad1;
    Vec3 axis1=ctr*rad2/rad1;
    Vec3 axis2=Vec3::Unit(2)*rad2;
    for(sizeType j=0; j<slice2; j++) {
      scalar angJ=2.0f*M_PI*(scalar(j)/scalar(slice2));
      Vec3 pt=ctr+axis1*cos(angJ)+axis2*sin(angJ);
      m.getV().push_back(pt);
    }
  }
  m.getI().clear();
  for(sizeType i=0; i<slice1; i++)
    for(sizeType j=0; j<slice2; j++) {
      m.getI().push_back(Vec3i(GI(i  ,j  ,slice1,slice2),
                               GI(i+1,j+1,slice1,slice2),
                               GI(i  ,j+1,slice1,slice2)));
      m.getI().push_back(Vec3i(GI(i  ,j  ,slice1,slice2),
                               GI(i+1,j  ,slice1,slice2),
                               GI(i+1,j+1,slice1,slice2)));
    }
  m.applyTrans();
  m.getDim()=3;
}
void MakeMesh::makeRing3D(ObjMesh& m,const scalar rad1,const scalar rad2,const scalar rad3,const sizeType slice)
{
  m.getV().clear();
  for(sizeType i=0; i<slice; i++) {
    scalar angI=2.0f*M_PI*(scalar(i)/scalar(slice));
    Vec3 ctr=Vec3(cos(angI),sin(angI),0.0f)*rad1;
    Vec3 axis1=ctr*rad2/rad1;
    Vec3 axis2=Vec3::Unit(2)*rad3;

    m.getV().push_back(ctr-axis1-axis2);
    m.getV().push_back(ctr+axis1-axis2);
    m.getV().push_back(ctr+axis1+axis2);
    m.getV().push_back(ctr-axis1+axis2);
  }
  m.getI().clear();
  for(sizeType i=0; i<slice; i++)
    for(sizeType j=0; j<4; j++) {
      m.getI().push_back(Vec3i(GI(i  ,j  ,slice,4),
                               GI(i+1,j+1,slice,4),
                               GI(i  ,j+1,slice,4)));
      m.getI().push_back(Vec3i(GI(i  ,j  ,slice,4),
                               GI(i+1,j  ,slice,4),
                               GI(i+1,j+1,slice,4)));
    }
  m.applyTrans();
  m.getDim()=3;
}
void MakeMesh::makeGridWithHole(ObjMesh& m,Vec4i slice)
{
  m.getV().clear();
  m.getI().clear();
  std::map<sizeType,sizeType> vMap;
  slice.x()=(slice.x()+slice[2]+slice[3]-1)/(slice[2]+slice[3])*(slice[2]+slice[3])+slice[3];
  slice.y()=(slice.y()+slice[2]+slice[3]-1)/(slice[2]+slice[3])*(slice[2]+slice[3])+slice[3];
  for(sizeType x=0; x<=slice.x(); x++)
    for(sizeType y=0; y<=slice.y(); y++)
      if((x%(slice[2]+slice[3]))<=slice[3] || (y%(slice[2]+slice[3]))<=slice[3]) {
        vMap[GI(x  ,y  ,slice.x()+1,slice.y()+1)]=(sizeType)m.getV().size();
        m.getV().push_back(Vec3(scalar(x)/scalar(slice.x()),scalar(y)/scalar(slice.y()),0.0f));
      }
  for(sizeType x=0; x<=slice.x(); x++)
    for(sizeType y=0; y<=slice.y(); y++)
      if(x<slice.x() && y<slice.y()) {
        m.getI().push_back(
          Vec3i(GI(x  ,y  ,slice.x()+1,slice.y()+1,vMap),
                GI(x+1,y  ,slice.x()+1,slice.y()+1,vMap),
                GI(x+1,y+1,slice.x()+1,slice.y()+1,vMap)));
        if(m.getI().back()[0]==-1 || m.getI().back()[1]==-1 || m.getI().back()[2]==-1)
          m.getI().pop_back();

        m.getI().push_back(
          Vec3i(GI(x  ,y  ,slice.x()+1,slice.y()+1,vMap),
                GI(x+1,y+1,slice.x()+1,slice.y()+1,vMap),
                GI(x  ,y+1,slice.x()+1,slice.y()+1,vMap)));
        if(m.getI().back()[0]==-1 || m.getI().back()[1]==-1 || m.getI().back()[2]==-1)
          m.getI().pop_back();
      }
  m.smooth();
}
void MakeMesh::makeGrid(ObjMesh& m,const Vec2i& slice)
{
  m.getV().clear();
  m.getI().clear();
  for(sizeType x=0; x<=slice.x(); x++)
    for(sizeType y=0; y<=slice.y(); y++) {
      m.getV().push_back(Vec3(scalar(x)/scalar(slice.x()),scalar(y)/scalar(slice.y()),0.0f));
      if(x<slice.x() && y<slice.y()) {
        m.getI().push_back(
          Vec3i(GI(x  ,y  ,slice.x()+1,slice.y()+1),
                GI(x+1,y  ,slice.x()+1,slice.y()+1),
                GI(x+1,y+1,slice.x()+1,slice.y()+1)));
        m.getI().push_back(
          Vec3i(GI(x  ,y  ,slice.x()+1,slice.y()+1),
                GI(x+1,y+1,slice.x()+1,slice.y()+1),
                GI(x  ,y+1,slice.x()+1,slice.y()+1)));
      }
    }
  m.smooth();
}
//two sphere mesh
bool MakeMesh::makeTwoSphereMesh3D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar lenY,const sizeType slice)
{
  return makeTwoSphereMesh(3,m,rad1,rad2,lenY,slice);
}
bool MakeMesh::makeTwoSphereMesh2D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar lenY,const sizeType slice)
{
  bool canDo=makeTwoSphereMesh(2,m,rad1,rad2,lenY,slice);
  if(canDo && m.getVolume()<0)
    m.insideOut();
  return canDo;
}
bool MakeMesh::makeTwoSphereMesh(int dim,ObjMesh& m,const scalar rad1,const scalar rad2,scalar lenY,const sizeType slice)
{
  m.getDim()=dim;
  m.getV().clear();
  m.getI().clear();
  Vec3 x=Vec3::Zero(),y=Vec3::Zero();
  ASSERT(rad1>0 && rad2>0 && lenY>0)
  scalar theta=0;
  if((theta=getTheta(rad1,rad2,lenY))==0)
    return false;
  std::vector<sizeType> id0=addSpherePatch(dim,m,-Vec3::UnitY()*rad1,theta,slice,x,y,Vec3::Zero());
  std::vector<sizeType> id1=addSpherePatch(dim,m, Vec3::UnitY()*rad2,M_PI-theta,slice,x,y,Vec3::UnitY()*lenY);
  hook(dim,m,id0,id1,true);
  m.smooth();
  if(dim==3)
    m.makeUniform();
  return true;
}
//three sphere mesh
void addLine(ObjMesh& m,const Vec2& a,const Vec2& b)
{
  sizeType off=(sizeType)m.getV().size();
  m.getV().push_back(Vec3(a[0],a[1],0));
  m.getV().push_back(Vec3(b[0],b[1],0));
  m.getI().push_back(Vec3i(off,off+1,-1));
}
bool MakeMesh::makeThreeSphereMesh3D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar rad3,const Vec2& ctr1,const Vec2& ctr2,const Vec2& ctr3,const sizeType slice)
{
  bool canDo=makeThreeSphereMesh(3,m,rad1,rad2,rad3,ctr1,ctr2,ctr3,slice);
  if(canDo && m.getVolume()<0)
    m.insideOut();
  return canDo;
}
bool MakeMesh::makeThreeSphereMesh2D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar rad3,const Vec2& ctr1,const Vec2& ctr2,const Vec2& ctr3,const sizeType slice)
{
  return makeThreeSphereMesh(2,m,rad1,rad2,rad3,ctr1,ctr2,ctr3,slice);
}
bool MakeMesh::makeThreeSphereMesh(int dim,ObjMesh& m,scalar rad1,scalar rad2,scalar rad3,Vec2 ctr1,Vec2 ctr2,Vec2 ctr3,const sizeType slice)
{
  {
    Vec2 c21=ctr2-ctr1,c31=ctr3-ctr1;
    if(Vec3(c21[0],c21[1],0).cross(Vec3(c31[0],c31[1],0))[2] < 0) {
      std::swap(ctr2,ctr3);
      std::swap(rad2,rad3);
    }
  }
  //compute normal
  Vec3 nn1,nn2;
  {
    Mat2 DV;
    DV.row(0)=ctr1-ctr2;
    DV.row(1)=ctr1-ctr3;
    nn1.segment<2>(0)=nn2.segment<2>(0)=DV.inverse()*Vec2(rad1-rad2,rad1-rad3);
    scalar zz=1-nn1.segment<2>(0).squaredNorm();
    if(zz<=0)
      return false;
    nn2.segment<2>(0)*=-1;
    nn1[2]=nn2[2]=std::sqrt(zz);
  }
  //create sphere patches
  m.getDim()=dim;
  m.getV().clear();
  m.getI().clear();
//#define DEBUG_THREE_SPHERE
#ifdef DEBUG_THREE_SPHERE
  addLine(m,ctr1,ctr2);
  addLine(m,ctr2,ctr3);
  addLine(m,ctr1,ctr3);
#endif
  scalar d12=(ctr1-ctr2).norm(),d13=(ctr1-ctr3).norm(),d23=(ctr2-ctr3).norm();
  std::pair<std::vector<sizeType>,std::vector<sizeType> > id1,id2,id3;
  {
    //ctr1
    scalar theta21=getTheta(rad2,rad1,d12);
    scalar theta31=getTheta(rad3,rad1,d13);
    if(theta21==0 || theta31==0)
      return false;
    Vec2 n21=rot2D(-theta21)*(ctr2-ctr1).normalized();
    Vec2 n31=rot2D( theta31)*(ctr3-ctr1).normalized();
    id1=addSpherePatch(dim,m,n31*rad1,n21*rad1,-nn1*rad1,nn2*rad1,slice,ctr1);
#ifdef DEBUG_THREE_SPHERE
    addLine(m,ctr1,ctr1+n21*rad1);
    addLine(m,ctr1,ctr1+n31*rad1);
#endif
  }
  {
    //ctr2
    scalar theta12=getTheta(rad1,rad2,d12);
    scalar theta32=getTheta(rad3,rad2,d23);
    if(theta12==0 || theta32==0)
      return false;
    Vec2 n12=rot2D( theta12)*(ctr1-ctr2).normalized();
    Vec2 n32=rot2D(-theta32)*(ctr3-ctr2).normalized();
    id2=addSpherePatch(dim,m,n12*rad2,n32*rad2,-nn1*rad2,nn2*rad2,slice,ctr2);
#ifdef DEBUG_THREE_SPHERE
    addLine(m,ctr2,ctr2+n12*rad2);
    addLine(m,ctr2,ctr2+n32*rad2);
#endif
  }
  {
    //ctr3
    scalar theta13=getTheta(rad1,rad3,d13);
    scalar theta23=getTheta(rad2,rad3,d23);
    if(theta13==0 || theta23==0)
      return false;
    Vec2 n13=rot2D(-theta13)*(ctr1-ctr3).normalized();
    Vec2 n23=rot2D( theta23)*(ctr2-ctr3).normalized();
    id3=addSpherePatch(dim,m,n23*rad3,n13*rad3,-nn1*rad3,nn2*rad3,slice,ctr3);
#ifdef DEBUG_THREE_SPHERE
    addLine(m,ctr3,ctr3+n13*rad3);
    addLine(m,ctr3,ctr3+n23*rad3);
#endif
  }
  if(dim==3) {
    m.getI().push_back(Vec3i(id1.first.front(),id2.first.front(),id3.first.front()));
    m.getI().push_back(Vec3i(id1.first.back(),id2.first.back(),id3.first.back()));
  }
  hook(dim,m,id1.second,id2.first,dim==2);
  hook(dim,m,id1.first,id3.second,dim==2);
  hook(dim,m,id2.second,id3.first,dim==2);
  m.smooth();
  if(dim==3)
    m.makeUniform();
  return true;
}
//helper
Mat2 MakeMesh::rot2D(scalar theta)
{
  Mat2 ret;
  ret << cos(theta),-sin(theta),sin(theta),cos(theta);
  return ret;
}
scalar MakeMesh::getTheta(scalar rad1,scalar rad2,scalar lenY)
{
  scalar theta=M_PI/2;
  if(std::abs(rad1-rad2)>=lenY) {
    INFO("Cannot create SphereMesh!")
    return 0;
  } else if(rad1>rad2)
    theta+=asin((rad1-rad2)/lenY);
  else theta-=asin((rad2-rad1)/lenY);
  return theta;
}
void MakeMesh::hook(int dim,ObjMesh& m,std::vector<sizeType>& id0,std::vector<sizeType>& id1,bool close)
{
  ASSERT(id0.size()==id1.size())
  for(sizeType i=0,nr=(sizeType)id0.size(); i<(close?nr:nr-1); i++)
    if(dim == 2) {
      m.getI().push_back(Vec3i(id0[i],id1[i],-1));
    } else {
      m.getI().push_back(Vec3i(id0[i],id0[(i+1)%nr],id1[(i+1)%nr]));
      m.getI().push_back(Vec3i(id0[i],id1[(i+1)%nr],id1[i]));
    }
}
std::vector<sizeType> MakeMesh::addSpherePatch(int dim,ObjMesh& mm,const Vec3& n,scalar theta,sizeType slice,Vec3& x,Vec3& y,const Vec3& off)
{
#define GI(X,Y) (Y)*slice+((X)%slice)+1
  ObjMesh m;
  m.getDim()=mm.getDim();
  scalar rad=n.norm();
  scalar res=rad*sin(theta)*M_PI*2/(scalar)slice;
  sizeType sliceY=(sizeType)std::ceil(rad*theta/res);
  if(x.isZero()) {
    sizeType id;
    n.cwiseAbs().minCoeff(&id);
    x=Vec3::Unit(id).cross(n).normalized();
    y=n.cross(x).normalized();
  }
  std::vector<sizeType> ret;
  if(dim==2) {
    for(sizeType iy=0; iy<=sliceY*2; iy++) {
      scalar theta0=theta*(1-(scalar)iy/(scalar)sliceY);
      m.getV().push_back(n.normalized()*cos(theta0)+y*sin(theta0));
      m.getV().back()=m.getV().back()*n.norm()+off;
      if(iy<sliceY*2)
        m.getI().push_back(Vec3i(iy,iy+1,-1));
    }
    //hook
    sizeType vOff=(sizeType)mm.getV().size();
    ret.push_back(vOff);
    ret.push_back(vOff+sliceY*2);
  } else {
    //vertices
    m.getV().push_back(n+off);
    for(sizeType iy=0; iy<sliceY; iy++) {
      scalar theta0=theta*(1-(scalar)iy/(scalar)sliceY);
      for(sizeType ix=0; ix<slice; ix++) {
        scalar phi=M_PI*2*(scalar)ix/(scalar)slice;
        m.getV().push_back(n.normalized()*cos(theta0)+x*sin(theta0)*sin(phi)+y*sin(theta0)*cos(phi));
        m.getV().back()=m.getV().back()*n.norm()+off;
      }
    }
    //indices
    for(sizeType y=0; y<sliceY; y++)
      for(sizeType x=0; x<slice; x++)
        if(y==sliceY-1) {
          m.getI().push_back(Vec3i(GI(x,y),GI(x+1,y),0));
        } else {
          m.getI().push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
          m.getI().push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
        }
    //hook
    sizeType vOff=(sizeType)mm.getV().size();
    for(sizeType x=0; x<slice; x++)
      ret.push_back(vOff+1+x);
  }
  mm.addMesh(m);
  return ret;
#undef GI
}
std::pair<std::vector<sizeType>,std::vector<sizeType> > MakeMesh::addSpherePatch(int dim,ObjMesh& mm,const Vec2& n1,const Vec2& n2,const Vec3& nn1,const Vec3& nn2,sizeType slice,const Vec2& off)
{
  ObjMesh m;
  scalar rad=n1.norm();
  m.getDim()=mm.getDim();
  std::vector<sizeType> id1;
  std::vector<sizeType> id2;
  Vec3 axis=Vec3(n1[0],n1[1],0).cross(Vec3(n2[0],n2[1],0));
  scalar angle=acos(Vec3(n1[0],n1[1],0).dot(Vec3(n2[0],n2[1],0))/rad/rad);
  if(axis[2]<0)
    angle=2*M_PI-angle;
  if(dim==2) {
    scalar res=rad*M_PI/(scalar)slice;
    sizeType sliceY=(sizeType)std::ceil(rad*angle/res);
    for(sizeType iy=0; iy<=sliceY; iy++) {
      scalar theta0=angle*(scalar)iy/(scalar)sliceY;
      m.getV().push_back(Eigen::AngleAxis<scalar>(theta0,Vec3::UnitZ()).toRotationMatrix()*Vec3(n1[0],n1[1],0));
      m.getV().back()+=Vec3(off[0],off[1],0);
      if(iy<sliceY)
        m.getI().push_back(Vec3i(iy,iy+1,-1));
    }
    //hook
    sizeType vOff=(sizeType)mm.getV().size();
    id1.push_back(vOff);
    id2.push_back(vOff+sliceY);
  } else {
    Vec3 nn;
    scalar alpha,beta;
    sizeType lastSliceX=-1;
    sizeType lastVOff=-1;
    for(sizeType iy=0; iy<=slice*2; iy++) {
      //indices
      sizeType sliceX;
      sizeType vOff=(sizeType)m.getV().size();
      if(iy<=slice) {
        sliceX=iy;
        alpha=scalar(iy)/scalar(slice);
        nn=nn1;
      } else {
        sliceX=slice*2-iy;
        alpha=scalar(slice*2-iy)/scalar(slice);
        nn=nn2;
      }
      if(lastSliceX>=0) {
        sizeType two=lastSliceX>sliceX?0:1;
        sizeType indexLast=0,indexCurr=0;
        while(true) {
          if(two==0) {
            if(indexLast+1>lastSliceX)
              break;
            m.getI().push_back(Vec3i(lastVOff+indexLast+0,lastVOff+indexLast+1,vOff+indexCurr));
            indexLast++;
          } else {
            if(indexCurr+1>sliceX)
              break;
            m.getI().push_back(Vec3i(vOff+indexCurr+0,vOff+indexCurr+1,lastVOff+indexLast));
            indexCurr++;
          }
          two=1-two;
        }
      }
      lastSliceX=sliceX;
      lastVOff=vOff;
      //vertices
      for(sizeType ix=0; ix<=sliceX; ix++) {
        beta=scalar(ix)/std::max(scalar(sliceX),scalar(1));
        Vec3 v=Eigen::AngleAxis<scalar>(angle*beta,Vec3::UnitZ()).toRotationMatrix()*Vec3(n1[0],n1[1],0);
        Eigen::AngleAxis<scalar> aa(ScalarUtil<scalar>::ScalarQuat::FromTwoVectors(nn,v));
        Eigen::AngleAxis<scalar> bb(aa.angle()*alpha,aa.axis());
        v=bb.toRotationMatrix()*nn;
        if(ix==0)
          id1.push_back((sizeType)(m.getV().size()+mm.getV().size()));
        if(ix==sliceX)
          id2.push_back((sizeType)(m.getV().size()+mm.getV().size()));
        m.getV().push_back(v+Vec3(off[0],off[1],0));
      }
    }
  }
  mm.addMesh(m);
  return std::make_pair(id1,id2);
}
sizeType MakeMesh::GI(sizeType i,sizeType j,sizeType si,sizeType sj,const std::map<sizeType,sizeType>& vMap)
{
  std::map<sizeType,sizeType>::const_iterator it=vMap.find(GI(i,j,si,sj));
  if(it == vMap.end())
    return -1;
  else return it->second;
}
sizeType MakeMesh::GI(sizeType i,sizeType j,sizeType si,sizeType sj)
{
  return (i%si)*sj+(j%sj);
}
