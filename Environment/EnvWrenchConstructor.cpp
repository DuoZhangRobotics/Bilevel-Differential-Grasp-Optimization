#include "EnvWrenchConstructor.h"
#ifdef TRAJ_OPT_SUPPORT
#include <TrajOpt/PBDDynamicsSequence.h>
#include <TrajOpt/NEDynamicsSequence.h>
#endif
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Articulated/NEArticulatedGradientInfo.h>
#include <Utils/CrossSpatialUtil.h>
#include <CommonFile/MakeMesh.h>

USE_PRJ_NAMESPACE

//EnvWrenchConstructor
template <typename T>
EnvWrenchConstructor<T>::EnvWrenchConstructor(std::shared_ptr<Environment<T>> env):_env(env) {}
#ifdef TRAJ_OPT_SUPPORT
template <typename T>
void EnvWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,sizeType id,const PBDDynamicsSequence<T,3>& dyn,bool grad)
{
  operator()(externalWrench,[&](sizeType JID) {
    return dyn.getJointTrans(id,JID);
  },grad);
}
template <typename T>
void EnvWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,sizeType id,const PBDDynamicsSequence<T,2>& dyn,bool grad)
{
  operator()(externalWrench,[&](sizeType JID) {
    return dyn.getJointTrans(id,JID);
  },grad);
}
#endif
template <typename T>
void EnvWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,const PBDArticulatedGradientInfo<T>& info,bool grad)
{
  operator()(externalWrench,[&](sizeType JID) {
    return TRANSI(info._TM,JID);
  },grad);
}
template <typename T>
void EnvWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,const NEArticulatedGradientInfo<T>& info,bool grad)
{
  operator()(externalWrench,[&](sizeType JID) {
    return fromSpatial<T>(spatialInv<T>(TRANSI6(info._invTM,JID)));
  },grad);
}
template <typename T>
void EnvWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>&,std::function<Mat3X4T(sizeType)>,bool)
{
  FUNCTION_NOT_IMPLEMENTED
}
//C0FloorWrenchConstructor
template <typename T>
C0EnvWrenchConstructor<T>::C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n,const Vec3T& g,sizeType nrDir,T mu,T coef)
  :C0EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentExact<T>()),g,nrDir,mu,coef)
{
  std::shared_ptr<EnvironmentExact<T>> env=std::dynamic_pointer_cast<EnvironmentExact<T>>(_env);
  ASSERT(env)
  env->createStair(x,y,x0,z0,slope,n);
}
template <typename T>
C0EnvWrenchConstructor<T>::C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res,const Vec3T& g,sizeType nrDir,T mu,T coef)
  :C0EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentExact<T>()),g,nrDir,mu,coef)
{
  std::shared_ptr<EnvironmentExact<T>> env=std::dynamic_pointer_cast<EnvironmentExact<T>>(_env);
  ASSERT(env)
  env->createHills(x,y,h,res);
}
template <typename T>
C0EnvWrenchConstructor<T>::C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z,const Vec3T& g,sizeType nrDir,T mu,T coef)
  :C0EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentExact<T>()),g,nrDir,mu,coef)
{
  std::shared_ptr<EnvironmentExact<T>> env=std::dynamic_pointer_cast<EnvironmentExact<T>>(_env);
  ASSERT(env)
  env->createZigZag(wid,len,off,n,z);
}
template <typename T>
C0EnvWrenchConstructor<T>::C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD z,const Vec3T& g,sizeType nrDir,T mu,T coef)
  :C0EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentExact<T>()),g,nrDir,mu,coef)
{
  std::shared_ptr<EnvironmentExact<T>> env=std::dynamic_pointer_cast<EnvironmentExact<T>>(_env);
  ASSERT(env)
  env->createFloor(x,y,z);
}
template <typename T>
C0EnvWrenchConstructor<T>::C0EnvWrenchConstructor(const ArticulatedBody& body,const Vec4d& plane,const Vec3T& g,sizeType nrDir,T mu,T coef)
  :C0EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentExact<T>()),g,nrDir,mu,coef)
{
  std::shared_ptr<EnvironmentExact<T>> env=std::dynamic_pointer_cast<EnvironmentExact<T>>(_env);
  ASSERT(env)
  env->createFloor(plane);
}
template <typename T>
C0EnvWrenchConstructor<T>::C0EnvWrenchConstructor(const ArticulatedBody& body,const std::string& path,bool is2D,scalarD dxMul,const Vec3T& g,sizeType nrDir,T mu,T coef)
  :C0EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(path,is2D,dxMul)),g,nrDir,mu,coef) {}
template <typename T>
C0EnvWrenchConstructor<T>::C0EnvWrenchConstructor(const ArticulatedBody& body,std::shared_ptr<Environment<T>> env,const Vec3T& g,sizeType nrDir,T mu,T coef):EnvWrenchConstructor<T>(env),_body(body)
{
  sizeType id=-1;
  _frm.col(0)=-g;
  _frm.col(0)/=std::sqrt(_frm.col(0).squaredNorm());
  for(sizeType i=0; i<3; i++)
    if(id==-1 || std::abs(_frm.col(0)[i])<std::abs(_frm.col(0)[id]))
      id=i;
  _frm.col(1)=Vec3T::Unit(id).cross(_frm.col(0));
  _frm.col(1)/=std::sqrt(_frm.col(1).squaredNorm());
  _frm.col(2)=_frm.col(0).cross(_frm.col(1));

  _B.resize(3,nrDir);
  for(sizeType i=0; i<nrDir; i++) {
    T angle=M_PI*2*i/nrDir;
    _B.col(i)=_frm*Vec3T(1,std::cos(angle)*mu,std::sin(angle)*mu)*coef;
  }
}
template <typename T>
void C0EnvWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool)
{
  //INFO("Using C0Floor!")
  externalWrench.clear();
  for(sizeType i=0; i<(sizeType)_externalForces.size(); i++) {
    const EndEffectorBounds& ee=_externalForces[i];
    Mat3X4T t=tfunc(ee._JID.back());
    Vec3T pt=ROT(t)*ee._localPos.template cast<T>()+CTR(t),Rb,n;
    Mat3T R;
    T phi=_env->phi(pt);
    ExternalWrench<T> eeW;
    for(sizeType r=0; r<3; r++)
      for(sizeType c=0; c<4; c++)
        eeW._DBDX[r][c].setZero(6,_B.cols());
    if(phi<ee._phi0) {
      R=rotateVec(_frm.col(0),n=_env->phiGrad(pt));
      eeW._jointId=ee._JID.back();
      eeW._B.resize(6,_B.cols());
      for(sizeType b=0; b<_B.cols(); b++) {
        Rb=R*_B.col(b);
        eeW._B.col(b)=concat<Vec>((pt-n*ee._phi0).cross(Rb),Rb);
      }
    } else {
      eeW._jointId=ee._JID.back();
      eeW._B.setZero(6,_B.cols());
    }
    externalWrench.push_back(eeW);
  }
}
template <typename T>
typename C0EnvWrenchConstructor<T>::Mat3T C0EnvWrenchConstructor<T>::rotateVec(const Vec3T a,const Vec3T& b,Mat3T* dRdX,Mat3T* dRdY,Mat3T* dRdZ)
{
  //input
  //T ax;
  //T ay;
  //T az;
  //T bx;
  //T by;
  //T bz;
  Mat3T R;

  //temp
  T tt1;
  T tt2;
  T tt3;
  T tt4;
  T tt5;
  T tt6;
  T tt7;
  T tt8;
  T tt9;
  T tt10;
  T tt11;
  T tt12;
  T tt13;
  T tt14;
  T tt15;
  T tt16;
  T tt17;
  T tt18;
  T tt19;
  T tt20;
  T tt21;
  T tt22;
  T tt23;
  T tt24;
  T tt25;
  T tt26;
  T tt27;
  T tt28;
  T tt29;
  T tt30;
  T tt31;
  T tt32;
  T tt33;
  T tt34;
  T tt35;
  T tt36;
  T tt37;
  T tt38;
  T tt39;
  T tt40;
  T tt41;
  T tt42;

  tt1=a[2]*b[2]+a[1]*b[1]+a[0]*b[0]+1;
  tt2=1/tt1;
  tt3=a[1]*b[0];
  tt4=-a[0]*b[1];
  tt5=tt4+tt3;
  tt6=-a[1]*b[0];
  tt7=a[0]*b[1];
  tt8=tt7+tt6;
  tt9=tt5*tt8;
  tt10=a[2]*b[0];
  tt11=-a[0]*b[2];
  tt12=tt11+tt10;
  tt13=-a[2]*b[0];
  tt14=a[0]*b[2];
  tt15=tt14+tt13;
  tt16=tt12*tt15;
  tt17=tt16+tt9;
  tt18=-a[2]*b[1];
  tt19=a[1]*b[2];
  tt20=tt19+tt18;
  tt21=a[2]*b[1];
  tt22=-a[1]*b[2];
  tt23=tt22+tt21;
  tt24=tt23*tt20;
  tt25=tt24+tt9;
  tt26=tt24+tt16;
  tt27=-a[1]*tt5;
  tt28=a[1]*tt8;
  tt29=-a[2]*tt12;
  tt30=a[2]*tt15;
  tt31=1/pow(tt1,2);
  tt32=-a[1];
  tt33=-a[2];
  tt34=a[0]*tt5;
  tt35=-a[0]*tt8;
  tt36=-a[0];
  tt37=-a[2]*tt23;
  tt38=a[2]*tt20;
  tt39=a[0]*tt12;
  tt40=-a[0]*tt15;
  tt41=a[1]*tt23;
  tt42=-a[1]*tt20;
  R(0,0)=tt2*tt17+1;
  R(0,1)=tt12*tt20*tt2+tt4+tt3;
  R(0,2)=tt5*tt23*tt2+tt11+tt10;
  R(1,0)=tt15*tt23*tt2+tt7+tt6;
  R(1,1)=tt2*tt25+1;
  R(1,2)=tt8*tt12*tt2+tt22+tt21;
  R(2,0)=tt8*tt20*tt2+tt14+tt13;
  R(2,1)=tt5*tt15*tt2+tt19+tt18;
  R(2,2)=tt2*tt26+1;
  if(dRdX) {
    (*dRdX)(0,0)=tt2*(tt27+tt28+tt29+tt30)-a[0]*tt31*tt17;
    (*dRdX)(0,1)=a[2]*tt20*tt2-a[0]*tt12*tt20*tt31+a[1];
    (*dRdX)(0,2)=a[1]*tt23*tt2-a[0]*tt5*tt23*tt31+a[2];
    (*dRdX)(1,0)=(-a[2]*tt23*tt2)-a[0]*tt15*tt23*tt31+tt32;
    (*dRdX)(1,1)=(tt27+tt28)*tt2-a[0]*tt31*tt25;
    (*dRdX)(1,2)=(-a[1]*tt12*tt2)+a[2]*tt8*tt2-a[0]*tt8*tt12*tt31;
    (*dRdX)(2,0)=(-a[1]*tt20*tt2)-a[0]*tt8*tt20*tt31+tt33;
    (*dRdX)(2,1)=a[1]*tt15*tt2-a[2]*tt5*tt2-a[0]*tt5*tt15*tt31;
    (*dRdX)(2,2)=tt2*(tt29+tt30)-a[0]*tt31*tt26;
  }
  if(dRdY) {
    (*dRdY)(0,0)=(tt34+tt35)*tt2-a[1]*tt31*tt17;
    (*dRdY)(0,1)=(-a[2]*tt12*tt2)-a[1]*tt12*tt20*tt31+tt36;
    (*dRdY)(0,2)=(-a[0]*tt23*tt2)+a[2]*tt5*tt2-a[1]*tt5*tt23*tt31;
    (*dRdY)(1,0)=a[2]*tt15*tt2-a[1]*tt15*tt23*tt31+a[0];
    (*dRdY)(1,1)=tt2*(tt34+tt35+tt37+tt38)-a[1]*tt31*tt25;
    (*dRdY)(1,2)=a[0]*tt12*tt2-a[1]*tt8*tt12*tt31+a[2];
    (*dRdY)(2,0)=a[0]*tt20*tt2-a[2]*tt8*tt2-a[1]*tt8*tt20*tt31;
    (*dRdY)(2,1)=(-a[0]*tt15*tt2)-a[1]*tt5*tt15*tt31+tt33;
    (*dRdY)(2,2)=tt2*(tt37+tt38)-a[1]*tt31*tt26;
  }
  if(dRdZ) {
    (*dRdZ)(0,0)=tt2*(tt39+tt40)-a[2]*tt31*tt17;
    (*dRdZ)(0,1)=(-a[0]*tt20*tt2)+a[1]*tt12*tt2-a[2]*tt12*tt20*tt31;
    (*dRdZ)(0,2)=(-a[1]*tt5*tt2)-a[2]*tt5*tt23*tt31+tt36;
    (*dRdZ)(1,0)=a[0]*tt23*tt2-a[1]*tt15*tt2-a[2]*tt15*tt23*tt31;
    (*dRdZ)(1,1)=tt2*(tt41+tt42)-a[2]*tt31*tt25;
    (*dRdZ)(1,2)=(-a[0]*tt8*tt2)-a[2]*tt8*tt12*tt31+tt32;
    (*dRdZ)(2,0)=a[1]*tt8*tt2-a[2]*tt8*tt20*tt31+a[0];
    (*dRdZ)(2,1)=a[0]*tt5*tt2-a[2]*tt5*tt15*tt31+a[1];
    (*dRdZ)(2,2)=tt2*(tt39+tt40+tt41+tt42)-a[2]*tt31*tt26;
  }
  return R;
}
template <typename T>
void C0EnvWrenchConstructor<T>::debugRotateVec(sizeType nrIter)
{
  DEFINE_NUMERIC_DELTA_T(T)
  for(sizeType it=0; it<nrIter; it++) {
    Vec3T a=Vec3T::Random();
    a/=std::sqrt(a.squaredNorm());
    Vec3T b=Vec3T::Random();
    b/=std::sqrt(b.squaredNorm());
    Vec3T db=Vec3T::Random();
    Mat3T dRdX,dRdY,dRdZ,R=rotateVec(a,b,&dRdX,&dRdY,&dRdZ);
    DEBUG_GRADIENT("rotateVec",1,std::sqrt((R*a-b).squaredNorm()))
    Mat3T R2=rotateVec(a,b+db*DELTA);
    Mat3T dR=dRdX*db[0]+dRdY*db[1]+dRdZ*db[2];
    DEBUG_GRADIENT("dR",std::sqrt(dR.squaredNorm()),std::sqrt((dR-(R2-R)/DELTA).squaredNorm()))
  }
}
template <typename T>
void C0EnvWrenchConstructor<T>::writeEndEffectorVTK(const std::string& path) const
{
  ObjMesh eeMesh;
  PBDArticulatedGradientInfo<T> info(_body,Vec::Zero(_body.nrDOF()));
  for(sizeType i=0; i<(sizeType)_externalForces.size(); i++) {
    const EndEffectorBounds& ee=_externalForces[i];
    //draw EE
    ObjMesh s;
    MakeMesh::makeSphere3D(s,ee._phi0,16);
    PBDArticulatedGradientInfo<scalarD> info(_body,Cold::Zero(_body.nrDOF()));
    s.getPos()=(ROTI(info._TM,ee._JID.back())*ee._localPos+CTRI(info._TM,ee._JID.back())).cast<scalar>();
    s.applyTrans();
    eeMesh.addMesh(s);
  }
  eeMesh.writeVTK(path,true);
}
template <typename T>
void C0EnvWrenchConstructor<T>::writeVTK(const std::string& path) const
{
  _env->getMesh().writeVTK(path,true);
}
//C2FloorWrenchConstructor
template <typename T>
C2EnvWrenchConstructor<T>::C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n,const Vec3T& g,sizeType nrDir,T mu,T coef,T relax)
  :C2EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),g,nrDir,mu,coef,relax)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createStair(x,y,x0,z0,slope,n);
}
template <typename T>
C2EnvWrenchConstructor<T>::C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res,const Vec3T& g,sizeType nrDir,T mu,T coef,T relax)
  :C2EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),g,nrDir,mu,coef,relax)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createHills(x,y,h,res);
}
template <typename T>
C2EnvWrenchConstructor<T>::C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z,const Vec3T& g,sizeType nrDir,T mu,T coef,T relax)
  :C2EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),g,nrDir,mu,coef,relax)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createZigZag(wid,len,off,n,z);
}
template <typename T>
C2EnvWrenchConstructor<T>::C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD z,const Vec3T& g,sizeType nrDir,T mu,T coef,T relax)
  :C2EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),g,nrDir,mu,coef,relax)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createFloor(x,y,z);
}
template <typename T>
C2EnvWrenchConstructor<T>::C2EnvWrenchConstructor(const ArticulatedBody& body,const Vec4d& plane,const Vec3T& g,sizeType nrDir,T mu,T coef,T relax)
  :C2EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),g,nrDir,mu,coef,relax)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createFloor(plane);
}
template <typename T>
C2EnvWrenchConstructor<T>::C2EnvWrenchConstructor(const ArticulatedBody& body,const std::string& path,bool is2D,scalarD dxMul,const Vec3T& g,sizeType nrDir,T mu,T coef,T relax)
  :C2EnvWrenchConstructor<T>(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(path,is2D,dxMul)),g,nrDir,mu,coef,relax) {}
template <typename T>
C2EnvWrenchConstructor<T>::C2EnvWrenchConstructor(const ArticulatedBody& body,std::shared_ptr<Environment<T>> env,const Vec3T& g,sizeType nrDir,T mu,T coef,T relax)
  :C0EnvWrenchConstructor<T>(body,env,g,nrDir,mu,1),_coef(coef),_relax(relax)
{
  std::ostringstream oss;
  oss << "C2EnvWrenchConstructor: mu=" << mu << ", coef=" << coef << "!";
  INFO(oss.str().c_str())
  _crossB.resize(3,_B.cols()*3);
  for(sizeType i=0; i<_B.cols(); i++)
    _crossB.template block<3,3>(0,i*3)=cross<T>(_B.col(i));
}
template <typename T>
void C2EnvWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool grad)
{
  //INFO("Using C2Floor!")
  externalWrench.clear();
  for(sizeType i=0; i<(sizeType)_externalForces.size(); i++) {
    const EndEffectorBounds& ee=_externalForces[i];
    Vec4T localPosHomo(ee._localPos[0],ee._localPos[1],ee._localPos[2],1);
    Mat3X4T t=tfunc(ee._JID.back());
    Vec3T pt=t*localPosHomo,Rb,n,DRDptb,DPhiDPt;
    Mat3T R,DRDn[3],DnDpt;
    Vec6T DBDpt;
    T phi=_env->phi(pt,&DPhiDPt);
    ExternalWrench<T> eeW;
    for(sizeType r=0; r<3; r++)
      for(sizeType c=0; c<4; c++)
        eeW._DBDX[r][c].setZero(6,_B.cols());
    if(phi<ee._phi0 || _relax>0) {
      T depth=std::min<T>(phi-ee._phi0,0.);
      T coef=(depth*depth+_relax)*_coef;
      T DCoef=2*depth*_coef;
      R=rotateVec(_frm.col(0),n=_env->phiGrad(pt,&DnDpt),DRDn+0,DRDn+1,DRDn+2);
      //B
      eeW._jointId=ee._JID.back();
      eeW._B.resize(6,_B.cols());
      //std::cout << "Wrench: " << _coef << " " << depth << " " << phi << " " << ee._phi0 << std::endl;
      for(sizeType b=0; b<_B.cols(); b++) {
        Rb=R*_B.col(b);
        eeW._B.col(b)=concat<Vec>((pt-n*ee._phi0).cross(Rb),Rb);
        //DBDX
        if(grad)
          for(sizeType r=0; r<3; r++) {
            DRDptb=(DRDn[0]*DnDpt(0,r)+DRDn[1]*DnDpt(1,r)+DRDn[2]*DnDpt(2,r))*_B.col(b);
            DBDpt=concat<Vec>((pt-n*ee._phi0).cross(DRDptb),DRDptb);
            DBDpt.template segment<3>(0)+=(Vec3T::Unit(r)-DnDpt.col(r)*ee._phi0).cross(Rb);
            for(sizeType c=0; c<4; c++)
              eeW._DBDX[r][c].col(b)=(eeW._B.col(b)*DCoef*DPhiDPt[r]+DBDpt*coef)*localPosHomo[c];
          }
        eeW._B.col(b)*=coef;
      }
    } else {
      eeW._jointId=ee._JID.back();
      eeW._B.setZero(6,_B.cols());
    }
    externalWrench.push_back(eeW);
  }
}

PRJ_BEGIN
#define INSTANTIATE(T)  \
template class EnvWrenchConstructor<T>; \
template class C0EnvWrenchConstructor<T>;   \
template class C2EnvWrenchConstructor<T>;
INSTANTIATE(double)
#ifdef ALL_TYPES
INSTANTIATE(__float128)
INSTANTIATE(mpfr::mpreal)
#endif
PRJ_END
