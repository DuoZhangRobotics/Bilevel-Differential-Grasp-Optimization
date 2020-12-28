#include "Environment.h"
#include <Environment/ObjMeshGeomCellExact.h>
#include <Environment/MPQZIO.h>
#include <Utils/Utils.h>
#include <Utils/DebugGradient.h>
#include <CommonFile/MarchingCube3D.h>
#include <Articulated/ArticulatedBody.h>
#include <CommonFile/GridOp.h>

USE_PRJ_NAMESPACE

//Environment
template <typename T>
ObjMesh Environment<T>::getMeshProj2D() const
{
  ObjMesh m3D=getMesh(),ret;
  ret.setDim(2);
  for(sizeType i=0; i<(sizeType)m3D.getV().size(); i++)
    ret.getV().push_back(Vec3(m3D.getV(i)[0],m3D.getV(i)[1],0));
  for(sizeType i=0; i<(sizeType)m3D.getI().size(); i+=2)
    ret.getI().push_back(Vec3i(m3D.getI(i)[0],m3D.getI(i)[1],-1));
  return ret;
}
template <typename T>
void Environment<T>::debug(sizeType nrIter)
{
  Mat3T h;
  Vec3T pt,g,g2;
  BBox<scalarD> bb=getBB();
  DEFINE_NUMERIC_DELTA_T(T)
  for(sizeType i=0; i<nrIter; i++) {
    for(sizeType r=0; r<3; r++) {
      T alpha=RandEngine::randR(0,1);
      pt[r]=bb._minC[r]*(1-alpha)+bb._maxC[r]*alpha;
    }
    Vec3T delta=Vec3T::Random();
    T p=phi(pt,&g),p2=phi(pt+delta*DELTA);
    if(g.dot(delta)!=0) {
      DEBUG_GRADIENT("phiGrad",g.dot(delta),g.dot(delta)-(p2-p)/DELTA)
    }
    g=phiGrad(pt,&h),g2=phiGrad(pt+delta*DELTA);
    if(std::sqrt((h*delta).squaredNorm())!=0) {
      DEBUG_GRADIENT("phiHess",std::sqrt((h*delta).squaredNorm()),std::sqrt((h*delta-(g2-g)/DELTA).squaredNorm()))
    }
  }
}
//EnvironmentExact
template <typename T>
EnvironmentExact<T>::EnvironmentExact() {}
template <typename T>
EnvironmentExact<T>::EnvironmentExact(const ArticulatedBody& body)
{
  ObjMesh m,env;
  const StaticGeom& g=body.getGeomEnv();
  for(sizeType i=0; i<g.nrG(); i++) {
    g.getG(i).getMesh(m);
    env.addMesh(m);
  }
  ObjMeshGeomCell cell(Mat4::Identity(),env,0,true,false);
  _obj.reset(new ObjMeshGeomCellExact(cell));
}
template <typename T>
EnvironmentExact<T>::EnvironmentExact(const ObjMeshGeomCellExact& obj):_obj(new ObjMeshGeomCellExact(obj)) {}
template <typename T>
bool EnvironmentExact<T>::read(std::istream& is,IOData* dat)
{
  registerType<ObjMeshGeomCellExact>(dat);
  readBinaryData(_obj,is,dat);
  return is.good();
}
template <typename T>
bool EnvironmentExact<T>::write(std::ostream& os,IOData* dat) const
{
  registerType<ObjMeshGeomCellExact>(dat);
  writeBinaryData(_obj,os,dat);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> EnvironmentExact<T>::copy() const
{
  return std::shared_ptr<SerializableBase>(new EnvironmentExact);
}
template <typename T>
std::string EnvironmentExact<T>::type() const
{
  return typeid(EnvironmentExact).name();
}
template <typename T>
T EnvironmentExact<T>::phi(const Vec3T& x,Vec3T* g) const
{
  Mat3T hessian;
  Vec3T n,normal;
  Vec2i feat;
  T ret=_obj->closest<T>(x,n,normal,hessian,feat);
  if(g)
    *g=normal;
  return ret;
}
template <typename T>
typename EnvironmentExact<T>::Vec3T EnvironmentExact<T>::phiGrad(const Vec3T& x,Mat3T* h) const
{
  Mat3T hessian;
  Vec3T n,normal;
  Vec2i feat;
  _obj->closest<T>(x,n,normal,hessian,feat);
  if(h)
    *h=hessian;
  return normal;
}
template <typename T>
BBox<scalarD> EnvironmentExact<T>::getBB() const
{
  BBox<scalarD> bb;
  bb._minC=castRational<Vec3d,TriangleExact::PT>(_obj->getBB()._minC);
  bb._maxC=castRational<Vec3d,TriangleExact::PT>(_obj->getBB()._maxC);
  return bb;
}
template <typename T>
ObjMesh EnvironmentExact<T>::getMesh() const
{
  return getObj().getMesh();
}
template <typename T>
void EnvironmentExact<T>::debugVTK(const std::string& path,sizeType res) const
{
  recreate(path);
  _obj->getMesh().writeVTK(path+"/mesh.vtk",true);
  _obj->writePointDistVTK(path+"/dist.vtk",res);
}
template <typename T>
ObjMeshGeomCell EnvironmentExact<T>::createStair(scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n)
{
  ObjMesh m;
  m.getV().push_back(Vec3(0,-y,0));
  m.getV().push_back(Vec3(0, y,0));
  m.getV().push_back(Vec3(x,-y,0));
  m.getV().push_back(Vec3(x, y,0));
  for(sizeType i=1; i<=n; i++) {
    m.getV().push_back(Vec3(x+x0*(i-1)+slope*x0,-y,i*z0));
    m.getV().push_back(Vec3(x+x0*(i-1)+slope*x0, y,i*z0));
    m.getV().push_back(Vec3(x+x0*i,-y,i*z0));
    m.getV().push_back(Vec3(x+x0*i, y,i*z0));
  }
  for(sizeType i=0,off=0; i<n*2+1; i++,off+=2) {
    m.getI().push_back(Vec3i(0,3,1)+Vec3i::Constant(off));
    m.getI().push_back(Vec3i(0,2,3)+Vec3i::Constant(off));
  }
  m.smooth();
  ObjMeshGeomCell cell(Mat4::Identity(),m,0,true);
  _obj.reset(new ObjMeshGeomCellExact(cell));
  return cell;
}
template <typename T>
ObjMeshGeomCell EnvironmentExact<T>::createHills(scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res)
{
#define GI(X,Y) (X)*(resY+1)+(Y)
  ObjMesh m;
  sizeType resX=std::ceil(x*2*res);
  sizeType resY=std::ceil(y*2*res);
  m.getV().resize((resX+1)*(resY+1));
  for(sizeType r=0; r<=resX; r++)
    for(sizeType c=0; c<=resY; c++) {
      scalarD xx=-x+2*x*r/resX,yy=-y+2*y*c/resY;
      m.getV(GI(r,c))=Vec3(xx,yy,h(xx,yy));
      if(r<resX && c<resY) {
        m.getI().push_back(Vec3i(GI(r,c),GI(r+1,c),GI(r+1,c+1)));
        m.getI().push_back(Vec3i(GI(r,c),GI(r+1,c+1),GI(r,c+1)));
      }
    }
  m.smooth();
  ObjMeshGeomCell cell(Mat4::Identity(),m,0,true);
  _obj.reset(new ObjMeshGeomCellExact(cell));
  return cell;
#undef GI
}
template <typename T>
ObjMeshGeomCell EnvironmentExact<T>::createZigZag(scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z)
{
  ObjMesh m;
  for(sizeType pass=0; pass<2; pass++) {
    for(sizeType i=0; i<n; i++) {
      if(i==0)
        m.getV().push_back(Vec3(0,0,pass==0?-z:z));
      m.getV().push_back(Vec3(len*(i+0.5),off,pass==0?-z:z));
      m.getV().push_back(Vec3(len*(i+1.0),0  ,pass==0?-z:z));
    }
    for(sizeType i=0; i<n; i++) {
      if(i==0)
        m.getV().push_back(Vec3(len*n,wid,pass==0?-z:z));
      m.getV().push_back(Vec3(len*(n-i-0.5),wid+off,pass==0?-z:z));
      m.getV().push_back(Vec3(len*(n-i-1.0),wid+0  ,pass==0?-z:z));
    }
  }
  sizeType half=(sizeType)m.getV().size()/2;
  for(sizeType i=0; i<half; i++) {
    sizeType j=(i+1)%half;
    m.getI().push_back(Vec3i(i,j,j+half));
    m.getI().push_back(Vec3i(i,j+half,i+half));
  }
  m.smooth();
  ObjMeshGeomCell cell(Mat4::Identity(),m,0,true);
  _obj.reset(new ObjMeshGeomCellExact(cell));
  return cell;
}
template <typename T>
ObjMeshGeomCell EnvironmentExact<T>::createFloor(scalarD x,scalarD y,scalarD z)
{
  ObjMesh m;
  m.getV().push_back(Vec3(0,-y,0));
  m.getV().push_back(Vec3(0, y,0));
  m.getV().push_back(Vec3(x,-y,z));
  m.getV().push_back(Vec3(x, y,z));
  m.getI().push_back(Vec3i(0,3,1));
  m.getI().push_back(Vec3i(0,2,3));
  m.smooth();
  ObjMeshGeomCell cell(Mat4::Identity(),m,0,true);
  _obj.reset(new ObjMeshGeomCellExact(cell));
  return cell;
}
template <typename T>
ObjMeshGeomCell EnvironmentExact<T>::createFloor(const Vec4d& plane)
{
  Mat3d frm;
  sizeType id=-1;
  frm.col(0)=plane.template segment<3>(0);
  frm.col(0)/=std::sqrt(frm.col(0).squaredNorm());
  for(sizeType i=0; i<3; i++)
    if(id==-1 || std::abs(frm.col(0)[i])<std::abs(frm.col(0)[id]))
      id=i;
  frm.col(1)=Vec3d::Unit(id).cross(frm.col(0));
  frm.col(1)/=std::sqrt(frm.col(1).squaredNorm());
  frm.col(2)=frm.col(0).cross(frm.col(1));

  scalarD sz=std::sqrt(plane.template segment<3>(0).squaredNorm());
  Vec3d pt=-frm.col(0)*(plane[3]/sz);

  ObjMesh m;
  m.getV().push_back((pt-frm.col(1)*sz-frm.col(2)*sz).template cast<scalar>());
  m.getV().push_back((pt+frm.col(1)*sz-frm.col(2)*sz).template cast<scalar>());
  m.getV().push_back((pt+frm.col(1)*sz+frm.col(2)*sz).template cast<scalar>());
  m.getV().push_back((pt-frm.col(1)*sz+frm.col(2)*sz).template cast<scalar>());
  m.getI().push_back(Vec3i(0,1,2));
  m.getI().push_back(Vec3i(0,2,3));
  m.smooth();
  ObjMeshGeomCell cell(Mat4::Identity(),m,0,true);
  _obj.reset(new ObjMeshGeomCellExact(cell));
  return cell;
}
template <typename T>
const ObjMeshGeomCellExact& EnvironmentExact<T>::getObj() const
{
  return *_obj;
}
//EnvironmentCubic
template <typename T>
EnvironmentCubic<T>::EnvironmentCubic(scalarD dx):_dx(dx) {}
template <typename T>
EnvironmentCubic<T>::EnvironmentCubic(const ArticulatedBody& body,scalarD dx):_dx(dx)
{
  EnvironmentExact<T> env(body);
  buildDist(env);
}
template <typename T>
EnvironmentCubic<T>::EnvironmentCubic(const ObjMeshGeomCellExact& obj,scalarD dx):_dx(dx)
{
  EnvironmentExact<T> env(obj);
  buildDist(env);
}
template <typename T>
bool EnvironmentCubic<T>::read(std::istream& is,IOData* dat)
{
  _dist.read(is,dat);
  _mesh.readBinary(is);
  readBinaryData(_dx,is,dat);
  return is.good();
}
template <typename T>
bool EnvironmentCubic<T>::write(std::ostream& os,IOData* dat) const
{
  _dist.write(os,dat);
  _mesh.writeBinary(os);
  writeBinaryData(_dx,os,dat);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> EnvironmentCubic<T>::copy() const
{
  return std::shared_ptr<SerializableBase>(new EnvironmentCubic(_dx));
}
template <typename T>
std::string EnvironmentCubic<T>::type() const
{
  return typeid(EnvironmentCubic).name();
}
template <typename T>
T EnvironmentCubic<T>::phi(const Vec3T& x,Vec3T* g) const
{
  Vec3T pos;
  Vec3i posId=floorV(_dist.getIndexFrac(Vec3d(std::to_double(x[0]),std::to_double(x[1]),std::to_double(x[2]))));
  for(sizeType i=0; i<3; i++) {
    posId[i]=std::max<sizeType>(posId[i],1);
    posId[i]=std::min<sizeType>(posId[i],_dist.getNrPoint()[i]-3);
    pos[i]=(x[i]-_dist.getBB()._minC[i])*_dist.getInvCellSize()[i]-(_dist.isCenter()?0.5f:0.0f)-posId[i];
  }
  T val[4][4][4];
  for(sizeType x=0; x<4; x++)
    for(sizeType y=0; y<4; y++)
      for(sizeType z=0; z<4; z++)
        val[x][y][z]=_dist.get(posId+Vec3i(x-1,y-1,z-1));
  //gradient
  if(g) {
    g->coeffRef(0)=interp1D(pos[2],[&](sizeType zId) {
      return interp1D(pos[1],[&](sizeType yId) {
        return interp1DDiff(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[0];
    g->coeffRef(1)=interp1D(pos[2],[&](sizeType zId) {
      return interp1DDiff(pos[1],[&](sizeType yId) {
        return interp1D(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[1];
    g->coeffRef(2)=interp1DDiff(pos[2],[&](sizeType zId) {
      return interp1D(pos[1],[&](sizeType yId) {
        return interp1D(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[2];
  }
  //value
  return interp1D(pos[2],[&](sizeType zId) {
    return interp1D(pos[1],[&](sizeType yId) {
      return interp1D(pos[0],[&](sizeType xId) {
        return val[xId][yId][zId];
      });
    });
  });
}
template <typename T>
typename EnvironmentCubic<T>::Vec3T EnvironmentCubic<T>::phiGrad(const Vec3T& x,Mat3T* h) const
{
  Vec3T pos;
  Vec3i posId=floorV(_dist.getIndexFrac(Vec3d(std::to_double(x[0]),std::to_double(x[1]),std::to_double(x[2]))));
  for(sizeType i=0; i<3; i++) {
    posId[i]=std::max<sizeType>(posId[i],1);
    posId[i]=std::min<sizeType>(posId[i],_dist.getNrPoint()[i]-3);
    pos[i]=(x[i]-_dist.getBB()._minC[i])*_dist.getInvCellSize()[i]-(_dist.isCenter()?0.5f:0.0f)-posId[i];
  }
  T val[4][4][4];
  for(sizeType x=0; x<4; x++)
    for(sizeType y=0; y<4; y++)
      for(sizeType z=0; z<4; z++)
        val[x][y][z]=_dist.get(posId+Vec3i(x-1,y-1,z-1));
  //gradient
  Vec3T g;
  g[0]=interp1D(pos[2],[&](sizeType zId) {
    return interp1D(pos[1],[&](sizeType yId) {
      return interp1DDiff(pos[0],[&](sizeType xId) {
        return val[xId][yId][zId];
      });
    });
  })*_dist.getInvCellSize()[0];
  g[1]=interp1D(pos[2],[&](sizeType zId) {
    return interp1DDiff(pos[1],[&](sizeType yId) {
      return interp1D(pos[0],[&](sizeType xId) {
        return val[xId][yId][zId];
      });
    });
  })*_dist.getInvCellSize()[1];
  g[2]=interp1DDiff(pos[2],[&](sizeType zId) {
    return interp1D(pos[1],[&](sizeType yId) {
      return interp1D(pos[0],[&](sizeType xId) {
        return val[xId][yId][zId];
      });
    });
  })*_dist.getInvCellSize()[2];
  //hessian
  if(h) {
    //diag
    h->coeffRef(0,0)=interp1D(pos[2],[&](sizeType zId) {
      return interp1D(pos[1],[&](sizeType yId) {
        return interp1DDDiff(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[0]*_dist.getInvCellSize()[0];
    h->coeffRef(1,1)=interp1D(pos[2],[&](sizeType zId) {
      return interp1DDDiff(pos[1],[&](sizeType yId) {
        return interp1D(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[1]*_dist.getInvCellSize()[1];
    h->coeffRef(2,2)=interp1DDDiff(pos[2],[&](sizeType zId) {
      return interp1D(pos[1],[&](sizeType yId) {
        return interp1D(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[2]*_dist.getInvCellSize()[2];
    //offdiag
    h->coeffRef(0,1)=h->coeffRef(1,0)=interp1D(pos[2],[&](sizeType zId) {
      return interp1DDiff(pos[1],[&](sizeType yId) {
        return interp1DDiff(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[0]*_dist.getInvCellSize()[1];
    h->coeffRef(0,2)=h->coeffRef(2,0)=interp1DDiff(pos[2],[&](sizeType zId) {
      return interp1D(pos[1],[&](sizeType yId) {
        return interp1DDiff(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[0]*_dist.getInvCellSize()[2];
    h->coeffRef(1,2)=h->coeffRef(2,1)=interp1DDiff(pos[2],[&](sizeType zId) {
      return interp1DDiff(pos[1],[&](sizeType yId) {
        return interp1D(pos[0],[&](sizeType xId) {
          return val[xId][yId][zId];
        });
      });
    })*_dist.getInvCellSize()[1]*_dist.getInvCellSize()[2];
  }
  return g;
}
template <typename T>
BBox<scalarD> EnvironmentCubic<T>::getBB() const
{
  return _dist.getBB();
}
template <typename T>
ObjMesh EnvironmentCubic<T>::getMesh() const
{
#if 0
  ObjMeshD m;
  MarchingCube3D<scalarD> mc3D(m);
  mc3D.solve(_dist,0);

  ObjMesh mF;
  m.cast<scalar>(mF);
  mF.smooth();
  mF.insideOut();
  mF.smooth();
  return mF;
#else
  return _mesh;
#endif
}
template <typename T>
void EnvironmentCubic<T>::debugVTK(const std::string& path) const
{
  recreate(path);
  getMesh().writeVTK(path+"/mesh.vtk",true);
}
template <typename T>
void EnvironmentCubic<T>::writeDistVTK(const std::string& path) const
{
  GridOp<scalarD,scalarD>::write3DScalarGridVTK(path,_dist);
}
template <typename T>
void EnvironmentCubic<T>::createStair(scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n)
{
  EnvironmentExact<T> env;
  env.createStair(x,y,x0,z0,slope,n);
  buildDist(env);
}
template <typename T>
void EnvironmentCubic<T>::createHills(scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res)
{
  EnvironmentExact<T> env;
  env.createHills(x,y,h,res);
  buildDist(env);
}
template <typename T>
void EnvironmentCubic<T>::createZigZag(scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z)
{
  EnvironmentExact<T> env;
  env.createZigZag(wid,len,off,n,z);
  buildDist(env);
}
template <typename T>
void EnvironmentCubic<T>::createFloor(scalarD x,scalarD y,scalarD z)
{
  EnvironmentExact<T> env;
  env.createFloor(x,y,z);
  buildDist(env);
}
template <typename T>
void EnvironmentCubic<T>::createFloor(const Vec4d& plane)
{
  EnvironmentExact<T> env;
  env.createFloor(plane);
  buildDist(env);
}
//helper
template <typename T>
void EnvironmentCubic<T>::buildDist(const EnvironmentExact<T>& env)
{
  BBox<scalarD> bb;
  bb._minC=castRational<Vec3d,TriangleExact::PT>(env.getObj().getBB()._minC);
  bb._maxC=castRational<Vec3d,TriangleExact::PT>(env.getObj().getBB()._maxC);
  _dist.reset(ceilV(bb.getExtent()/_dx),bb,0,false);
  _dist.expand(Vec3i::Constant(5),0);
  _mesh=env.getMesh();

  sizeType nrLayer=0;
  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<_dist.getNrPoint()[0]; x++) {
    OMP_ATOMIC_
    nrLayer++;
    INFOV("Building %ld/%ld layer!",nrLayer,_dist.getNrPoint()[0])
    for(sizeType y=0; y<_dist.getNrPoint()[1]; y++)
      for(sizeType z=0; z<_dist.getNrPoint()[2]; z++)
        _dist.get(Vec3i(x,y,z))=std::to_double(env.phi(_dist.getPt(Vec3i(x,y,z)).template cast<T>()));
  }
}
template <typename T>
T EnvironmentCubic<T>::interp1D(T t,std::function<T(sizeType)> f)
{
  T t2=t*t,t3=t2*t;
  T C0=-(t-2*t2+t3)/2;
  T C1=(2-5*t2+3*t3)/2;
  T C2=(t+4*t2-3*t3)/2;
  T C3=(t3-t2)/2;
  return C0*f(0)+C1*f(1)+C2*f(2)+C3*f(3);
}
template <typename T>
T EnvironmentCubic<T>::interp1DDiff(T t,std::function<T(sizeType)> f)
{
  T t2=t*t;
  T C0=-(1-4*t+3*t2)/2;
  T C1=(9*t2-10*t)/2;
  T C2=(1+8*t-9*t2)/2;
  T C3=(3*t2-2*t)/2;
  return C0*f(0)+C1*f(1)+C2*f(2)+C3*f(3);
}
template <typename T>
T EnvironmentCubic<T>::interp1DDDiff(T t,std::function<T(sizeType)> f)
{
  T C0=2-3*t;
  T C1=9*t-5;
  T C2=4-9*t;
  T C3=3*t-1;
  return C0*f(0)+C1*f(1)+C2*f(2)+C3*f(3);
}
//EnvironmentHeight
template <typename T>
EnvironmentHeight<T>::EnvironmentHeight(scalarD dx):EnvironmentCubic<T>(dx) {}
template <typename T>
EnvironmentHeight<T>::EnvironmentHeight(const std::string& path,bool is2D,scalar dxMul):EnvironmentCubic<T>(0)
{
  ASSERT_MSGV(exists(path),"Cannot find: %s!",path.c_str())
  INFOV("Creating terrain from: %s, is2D=%d, dxMul=%f!",path.c_str(),is2D,dxMul)
  ObjMesh mesh;
  std::string tp;
  std::ifstream is(path);
  std::vector<sizeType> index;
  index.push_back(0);
  //vertex
  while(std::getline(is,tp)) {
    sizeType xpos=tp.find_first_of(' ');
    scalar x=std::atof(tp.substr(0,xpos).c_str());
    tp=tp.substr(xpos+1);
    scalar ypos=tp.find_first_of(' ');
    scalar y=std::atof(tp.substr(0,ypos).c_str());
    scalar z=std::atof(tp.substr(ypos+1).c_str());
    if(!mesh.getV().empty() && x!=mesh.getV().back()[0]) {
      if(is2D) {
        mesh.getV()[index.back()+1]=mesh.getV().back();
        mesh.getV().resize(index.back()+2);
      }
      index.push_back(mesh.getV().size());
      _dx=std::abs(x-mesh.getV().back()[0])*dxMul;
    }
    mesh.getV().push_back(Vec3(x,y,z));
  }
  //index
  for(sizeType i=0; i<(sizeType)index.size()-1; i++) {
    sizeType i0=index[i+0];
    sizeType i1=index[i+1];
    while(i0<index[i+1]-1) {
      mesh.getI().push_back(Vec3i(i0,i0+1,i1+1));
      mesh.getI().push_back(Vec3i(i0,i1+1,i1));
      i0++;
      i1++;
    }
  }
  mesh.smooth();
  ObjMeshGeomCell cell(Mat4::Identity(),mesh,0,true);
  buildDist(cell);
}
template <typename T>
EnvironmentHeight<T>::EnvironmentHeight(const ArticulatedBody& body,scalarD dx):EnvironmentCubic<T>(dx)
{
  ObjMesh m,env;
  const StaticGeom& g=body.getGeomEnv();
  for(sizeType i=0; i<g.nrG(); i++) {
    g.getG(i).getMesh(m);
    env.addMesh(m);
  }
  ObjMeshGeomCell cell(Mat4::Identity(),env,0,true,false);
  buildDist(cell);
}
template <typename T>
EnvironmentHeight<T>::EnvironmentHeight(const ObjMeshGeomCell& obj,scalarD dx):EnvironmentCubic<T>(dx)
{
  buildDist(obj);
}
template <typename T>
T EnvironmentHeight<T>::phi(const Vec3T& x,Vec3T* g) const
{
  Vec3T pos;
  Vec3i posId=floorV(_dist.getIndexFrac(Vec3d(std::to_double(x[0]),std::to_double(x[1]),std::to_double(x[2]))));
  for(sizeType i=0; i<2; i++) {
    posId[i]=std::max<sizeType>(posId[i],1);
    posId[i]=std::min<sizeType>(posId[i],_dist.getNrPoint()[i]-3);
    pos[i]=(x[i]-_dist.getBB()._minC[i])*_dist.getInvCellSize()[i]-(_dist.isCenter()?0.5f:0.0f)-posId[i];
  }
  T val[4][4];
  for(sizeType x=0; x<4; x++)
    for(sizeType y=0; y<4; y++)
      val[x][y]=_dist.get(posId+Vec3i(x-1,y-1,0));
  //gradient
  if(g) {
    g->coeffRef(0)=-interp1D(pos[1],[&](sizeType yId) {
      return interp1DDiff(pos[0],[&](sizeType xId) {
        return val[xId][yId];
      });
    })*_dist.getInvCellSize()[0];
    g->coeffRef(1)=-interp1DDiff(pos[1],[&](sizeType yId) {
      return interp1D(pos[0],[&](sizeType xId) {
        return val[xId][yId];
      });
    })*_dist.getInvCellSize()[1];
    g->coeffRef(2)=1;
  }
  //value
  return x[2]-interp1D(pos[1],[&](sizeType yId) {
    return interp1D(pos[0],[&](sizeType xId) {
      return val[xId][yId];
    });
  });
}
template <typename T>
typename EnvironmentHeight<T>::Vec3T EnvironmentHeight<T>::phiGrad(const Vec3T& x,Mat3T* h) const
{
  Vec3T pos;
  Vec3i posId=floorV(_dist.getIndexFrac(Vec3d(std::to_double(x[0]),std::to_double(x[1]),std::to_double(x[2]))));
  for(sizeType i=0; i<2; i++) {
    posId[i]=std::max<sizeType>(posId[i],1);
    posId[i]=std::min<sizeType>(posId[i],_dist.getNrPoint()[i]-3);
    pos[i]=(x[i]-_dist.getBB()._minC[i])*_dist.getInvCellSize()[i]-(_dist.isCenter()?0.5f:0.0f)-posId[i];
  }
  T val[4][4];
  for(sizeType x=0; x<4; x++)
    for(sizeType y=0; y<4; y++)
      val[x][y]=_dist.get(posId+Vec3i(x-1,y-1,0));
  //gradient
  Vec3T g;
  g[0]=interp1D(pos[1],[&](sizeType yId) {
    return interp1DDiff(pos[0],[&](sizeType xId) {
      return val[xId][yId];
    });
  })*_dist.getInvCellSize()[0];
  g[1]=interp1DDiff(pos[1],[&](sizeType yId) {
    return interp1D(pos[0],[&](sizeType xId) {
      return val[xId][yId];
    });
  })*_dist.getInvCellSize()[1];
  g[2]=1;
  T len=std::sqrt(g.squaredNorm()),gradLen=0;
  //hessian
  if(h) {
    Mat2T hessian;
    hessian(0,0)=interp1D(pos[1],[&](sizeType yId) {
      return interp1DDDiff(pos[0],[&](sizeType xId) {
        return val[xId][yId];
      });
    })*_dist.getInvCellSize()[0]*_dist.getInvCellSize()[0];
    hessian(1,1)=interp1DDDiff(pos[1],[&](sizeType yId) {
      return interp1D(pos[0],[&](sizeType xId) {
        return val[xId][yId];
      });
    })*_dist.getInvCellSize()[1]*_dist.getInvCellSize()[1];
    hessian(0,1)=hessian(1,0)=interp1DDiff(pos[1],[&](sizeType yId) {
      return interp1DDiff(pos[0],[&](sizeType xId) {
        return val[xId][yId];
      });
    })*_dist.getInvCellSize()[0]*_dist.getInvCellSize()[1];
    //fill
    h->setZero();
    for(sizeType d=0; d<2; d++) {
      gradLen=2*(g[0]*hessian(0,d)+g[1]*hessian(1,d));
      h->col(d)=Vec3T(-hessian(0,d),-hessian(1,d),0)/len;
      h->col(d)-=Vec3T(-g[0],-g[1],1)*(0.5*gradLen/len/len/len);
    }
  }
  return Vec3T(-g[0],-g[1],1)/len;
}
template <typename T>
BBox<scalarD> EnvironmentHeight<T>::getBB() const
{
  BBox<scalarD> bb=_dist.getBB();
  _dist.minMax(bb._minC[2],bb._maxC[2]);
  return bb;
}
template <typename T>
ObjMesh EnvironmentHeight<T>::getMesh() const
{
#if 1
#define GI(X,Y) (X)*_dist.getNrPoint()[1]+(Y)
  ObjMeshD m;
  for(sizeType x=0; x<_dist.getNrPoint()[0]; x++)
    for(sizeType y=0; y<_dist.getNrPoint()[1]; y++) {
      Vec3d pt=_dist.getPt(Vec3i(x,y,0));
      pt[2]=_dist.get(Vec3i(x,y,0));
      m.getV().push_back(pt);
    }
  for(sizeType x=0; x<_dist.getNrPoint()[0]-1; x++)
    for(sizeType y=0; y<_dist.getNrPoint()[1]-1; y++) {
      m.getI().push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x+1,y)));
      m.getI().push_back(Vec3i(GI(x,y),GI(x,y+1),GI(x+1,y+1)));
    }
  m.smooth();

  ObjMesh mF;
  m.cast<scalar>(mF);
  mF.insideOut();
  return mF;
#undef GI
#else
  return _mesh;
#endif
}
template <typename T>
void EnvironmentHeight<T>::createStair(scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n)
{
  EnvironmentExact<T> env;
  ObjMeshGeomCell cell=env.createStair(x,y,x0,z0,slope,n);
  buildDist(cell);
}
template <typename T>
void EnvironmentHeight<T>::createHills(scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res)
{
  EnvironmentExact<T> env;
  ObjMeshGeomCell cell=env.createHills(x,y,h,res);
  buildDist(cell);
}
template <typename T>
void EnvironmentHeight<T>::createZigZag(scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z)
{
  EnvironmentExact<T> env;
  ObjMeshGeomCell cell=env.createZigZag(wid,len,off,n,z);
  buildDist(cell);
}
template <typename T>
void EnvironmentHeight<T>::createFloor(scalarD x,scalarD y,scalarD z)
{
  EnvironmentExact<T> env;
  ObjMeshGeomCell cell=env.createFloor(x,y,z);
  buildDist(cell);
}
template <typename T>
void EnvironmentHeight<T>::createFloor(const Vec4d& plane)
{
  EnvironmentExact<T> env;
  ObjMeshGeomCell cell=env.createFloor(plane);
  buildDist(cell);
}
//helper
template <typename T>
void EnvironmentHeight<T>::buildDist(const ObjMeshGeomCell& cell)
{
  BBox<scalar> bb=cell.getBB();
  bb._minC[2]=bb._maxC[2]=0;
  _dist.reset(floorV(bb.getExtent()/_dx),bb,0,true);
  cell.getMesh(_mesh);

  bb=cell.getBB();
  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<_dist.getNrPoint()[0]; x++)
    for(sizeType y=0; y<_dist.getNrPoint()[1]; y++) {
      Vec3d pt=_dist.getPt(Vec3i(x,y,0));
      Vec3 x0(pt[0],pt[1],bb._minC[2]);
      Vec3 dir=Vec3::UnitZ()*bb.getExtent()[2];
      scalar s=cell.rayQuery(x0,dir);
      _dist.get(Vec3i(x,y,0))=x0[2]+dir[2]*s;
    }
}
//instance
PRJ_BEGIN
#define INSTANTIATE(T)  \
template class EnvironmentExact<T>; \
template class EnvironmentCubic<T>;   \
template class EnvironmentHeight<T>;
INSTANTIATE(double)
#ifdef ALL_TYPES
INSTANTIATE(__float128)
INSTANTIATE(mpfr::mpreal)
#endif
PRJ_END
