#include "GranularWrenchConstructor.h"
#include <CommonFile/MakeMesh.h>
#include <Utils/CrossSpatialUtil.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

template <typename T>
C2GranularWrenchConstructor<T>::C2GranularWrenchConstructor
(const ArticulatedBody& body,scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n,const std::string& pathInput,const std::string& pathWeight)
  :C2GranularWrenchConstructor(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),pathInput,pathWeight)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createStair(x,y,x0,z0,slope,n);
}
template <typename T>
C2GranularWrenchConstructor<T>::C2GranularWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res,const std::string& pathInput,const std::string& pathWeight)
  :C2GranularWrenchConstructor(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),pathInput,pathWeight)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createHills(x,y,h,res);
}
template <typename T>
C2GranularWrenchConstructor<T>::C2GranularWrenchConstructor(const ArticulatedBody& body,scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z,const std::string& pathInput,const std::string& pathWeight)
  :C2GranularWrenchConstructor(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),pathInput,pathWeight)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createZigZag(wid,len,off,n,z);
}
template <typename T>
C2GranularWrenchConstructor<T>::C2GranularWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD z,const std::string& pathInput,const std::string& pathWeight)
  :C2GranularWrenchConstructor(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),pathInput,pathWeight)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createFloor(x,y,z);
}
template <typename T>
C2GranularWrenchConstructor<T>::C2GranularWrenchConstructor(const ArticulatedBody& body,const Vec4d& plane,const std::string& pathInput,const std::string& pathWeight)
  :C2GranularWrenchConstructor(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(0.1)),pathInput,pathWeight)
{
  std::shared_ptr<EnvironmentHeight<T>> env=std::dynamic_pointer_cast<EnvironmentHeight<T>>(_env);
  ASSERT(env)
  env->createFloor(plane);
}
template <typename T>
C2GranularWrenchConstructor<T>::C2GranularWrenchConstructor(const ArticulatedBody& body,const std::string& path,bool is2D,scalarD dxMul,const std::string& pathInput,const std::string& pathWeight)
  :C2GranularWrenchConstructor(body,std::shared_ptr<Environment<T>>(new EnvironmentHeight<T>(path,is2D,dxMul)),pathInput,pathWeight) {}
template <typename T>
C2GranularWrenchConstructor<T>::C2GranularWrenchConstructor
(const ArticulatedBody& body,std::shared_ptr<Environment<T>> env,const std::string& pathInput,const std::string& pathWeight)
  :EnvWrenchConstructor<T>(env),_body(body)
{
  ASSERT_MSGV(exists(pathInput),"Input file (%s) does not exist!",pathInput.c_str())
  ASSERT_MSGV(exists(pathWeight),"Weight file (%s) does not exist!",pathWeight.c_str())
  std::vector<double> vars;
  {
    vars.clear();
    std::string tp;
    std::fstream newfile;
    newfile.open(pathInput,std::ios::in);
    if(newfile.is_open()) {
      while(getline(newfile,tp)) {
        vars.resize(vars.size()+2);
        std::istringstream(tp) >> vars[vars.size()-2] >> vars[vars.size()-1];
      }
      newfile.close();
    }
    _input=Eigen::Map<Eigen::Matrix<double,2,-1>>(&vars[0],2,vars.size()/2).template cast<T>();
  }
  //weight is ordered as: fx,fz,wy
  {
    vars.clear();
    std::string tp;
    std::fstream newfile;
    newfile.open(pathWeight,std::ios::in);
    if(newfile.is_open()) {
      while(getline(newfile,tp)) {
        vars.resize(vars.size()+3);
        std::istringstream(tp) >> vars[vars.size()-3] >> vars[vars.size()-2] >> vars[vars.size()-1];
      }
      newfile.close();
    }
    _weight=Eigen::Map<Eigen::Matrix<double,3,-1>>(&vars[0],3,vars.size()/3).template cast<T>();
  }
  //detect external force
  //this assumes that the collision mesh is used,
  //and this assumes that we are using: robosimian_caesar_new_all_active.urdf
  for(sizeType i=0; i<body.nrJ(); i++)
    if(body.children(i,true).empty()) {
      Vec3d zRange;
      EndEffectorBounds ee(i);
      if(beginsWith(body.joint(i)._name,"limb1"))
        zRange=Vec3d(-0.001f,0,0);
      else if(beginsWith(body.joint(i)._name,"limb2"))
        zRange=Vec3d( 0.001f,0,0);
      else if(beginsWith(body.joint(i)._name,"limb3"))
        zRange=Vec3d( 0.001f,0,0);
      else if(beginsWith(body.joint(i)._name,"limb4"))
        zRange=Vec3d(-0.001f,0,0);
      SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0,zRange);
      //In collision mesh we assume cylinder length is 0.15m, but Yifan assumes this length is 0.1572m
      ee._localPos.x()+=0.0f;
      ee._phi0=0.01f;
      _externalForces.push_back(ee);
      //Torque center is moved up by 0.15m
      _torqueCenters.push_back(ee._localPos.template cast<T>()-Vec3T(0.1572f,0,0));
      //Upward normal position
      _normalDirs.push_back(Vec3T(-1,0,0));
    }
}
template <typename T>
void C2GranularWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool grad)
{
  Vec2T w(5,1);
  operator()(externalWrench,tfunc,[&](const Vec2T& r,Vec2T* diff) {
    Vec2T wr=w.asDiagonal()*r;
    T rSqr=wr.squaredNorm();
    if(diff)
      *diff=w.asDiagonal()*wr*(std::log(rSqr)+1);
    return rSqr*std::log(rSqr)/2;
  },grad);
}
template <typename T>
void C2GranularWrenchConstructor<T>::operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,RBF kernel,bool grad)
{
  sizeType nVert=_weight.cols()/_input.cols();
  externalWrench.resize(_externalForces.size());
  for(sizeType i=0; i<(sizeType)_externalForces.size(); i++) {
    const EndEffectorBounds& ee=_externalForces[i];
    Mat3X4T t=tfunc(ee._JID.back()),DwDX,DPhiDX,DtDX,DtcDX[3];
    Vec4T localPosHomo(ee._localPos[0],ee._localPos[1],ee._localPos[2],1),torqueCenterHomo(_torqueCenters[i][0],_torqueCenters[i][1],_torqueCenters[i][2],1);
    Vec3T pt=t*localPosHomo,tc=t*torqueCenterHomo,n,DPhiDPt,Dt0,Dt1;
    Mat3T ctc=cross<T>(tc),R,DnDpt,dR[3];
    Vec2T DwDf;
    T phi=_env->phi(pt,&DPhiDPt);
    //create wrench
    ExternalWrench<T>& W=externalWrench[i];
    W._B.setZero(6,nVert);
    for(sizeType r=0; r<3; r++)
      for(sizeType c=0; c<4; c++)
        W._DBDX[r][c].setZero(6,nVert);
    W._jointId=ee.jointId();
    if(phi<0.1f) {
      n=_env->phiGrad(pt,&DnDpt);
      //compute feature
      Vec2T feature;
      feature[0]=phi;
      feature[1]=theta(n,ROT(t)*_normalDirs[i],&Dt0,&Dt1);
      DPhiDX=DPhiDPt*localPosHomo.transpose();
      DtDX=(DnDpt.transpose()*Dt0)*localPosHomo.transpose();
      ROT(DtDX)+=Dt1*_normalDirs[i].transpose();
      //fill wrench
      for(sizeType ptid=0; ptid<_input.cols(); ptid++) {
        T w=kernel(feature-_input.col(ptid),grad?&DwDf:NULL);
        if(grad)
          DwDX=DwDf[0]*DPhiDX+DwDf[1]*DtDX;
        for(sizeType vid=0; vid<nVert; vid++) {
          W._B(3,vid)+=_weight(0,vid*_input.cols()+ptid)*w;
          W._B(5,vid)+=_weight(1,vid*_input.cols()+ptid)*w;
          W._B(1,vid)+=_weight(2,vid*_input.cols()+ptid)*w;
          if(grad)
            for(sizeType r=0; r<3; r++)
              for(sizeType c=0; c<4; c++) {
                W._DBDX[r][c](3,vid)+=_weight(0,vid*_input.cols()+ptid)*DwDX(r,c);
                W._DBDX[r][c](5,vid)+=_weight(1,vid*_input.cols()+ptid)*DwDX(r,c);
                W._DBDX[r][c](1,vid)+=_weight(2,vid*_input.cols()+ptid)*DwDX(r,c);
              }
        }
      }
      //move torque center
      if(grad)
        for(sizeType vid=0; vid<nVert; vid++) {
          R=cross<T>(W._B.template block<3,1>(3,vid));
          DtcDX[0]=-R.row(0).transpose()*torqueCenterHomo.transpose();
          DtcDX[1]=-R.row(1).transpose()*torqueCenterHomo.transpose();
          DtcDX[2]=-R.row(2).transpose()*torqueCenterHomo.transpose();
          for(sizeType r=0; r<3; r++)
            for(sizeType c=0; c<4; c++) {
              W._DBDX[r][c](0,vid)+=ctc(0,0)*W._DBDX[r][c](3,vid)+ctc(0,1)*W._DBDX[r][c](4,vid)+ctc(0,2)*W._DBDX[r][c](5,vid)+DtcDX[0](r,c);
              W._DBDX[r][c](1,vid)+=ctc(1,0)*W._DBDX[r][c](3,vid)+ctc(1,1)*W._DBDX[r][c](4,vid)+ctc(1,2)*W._DBDX[r][c](5,vid)+DtcDX[1](r,c);
              W._DBDX[r][c](2,vid)+=ctc(2,0)*W._DBDX[r][c](3,vid)+ctc(2,1)*W._DBDX[r][c](4,vid)+ctc(2,2)*W._DBDX[r][c](5,vid)+DtcDX[2](r,c);
            }
        }
      W._B.block(0,0,3,nVert)+=ctc*W._B.block(3,0,3,nVert);
      //rotate according to slope of normal
      R=rotateVecPlanar(Vec3T::UnitZ(),n,NULL,dR);
      if(grad)
        for(sizeType vid=0; vid<nVert; vid++) {
          ctc.col(0)=dR[0]*W._B.template block<3,1>(3,vid);
          ctc.col(1)=dR[1]*W._B.template block<3,1>(3,vid);
          ctc.col(2)=dR[2]*W._B.template block<3,1>(3,vid);
          ctc*=DnDpt;
          DtcDX[0]=ctc.row(0).transpose()*localPosHomo.transpose();
          DtcDX[1]=ctc.row(1).transpose()*localPosHomo.transpose();
          DtcDX[2]=ctc.row(2).transpose()*localPosHomo.transpose();
          for(sizeType r=0; r<3; r++)
            for(sizeType c=0; c<4; c++) {
              DPhiDPt=Vec3T(W._DBDX[r][c](3,vid),W._DBDX[r][c](4,vid),W._DBDX[r][c](5,vid));
              Dt0=R*DPhiDPt;
              W._DBDX[r][c](3,vid)=Dt0[0]+DtcDX[0](r,c);
              W._DBDX[r][c](4,vid)=Dt0[1]+DtcDX[1](r,c);
              W._DBDX[r][c](5,vid)=Dt0[2]+DtcDX[2](r,c);
            }
        }
      W._B.block(3,0,3,nVert)=(R*W._B.block(3,0,3,nVert)).eval();
    }
  }
}
template <typename T>
void C2GranularWrenchConstructor<T>::writeEndEffectorVTK(const std::string& path,bool torqueCenter) const
{
  ObjMesh eeMesh;
  PBDArticulatedGradientInfo<T> info(_body,Vec::Zero(_body.nrDOF()));
  for(sizeType i=0; i<(sizeType)_externalForces.size(); i++) {
    const EndEffectorBounds& ee=_externalForces[i];
    //draw EE
    ObjMesh s;
    Vec3d local;
    MakeMesh::makeSphere3D(s,ee._phi0,16);
    PBDArticulatedGradientInfo<scalarD> info(_body,Cold::Zero(_body.nrDOF()));
    if(torqueCenter)
      local=_torqueCenters[i].unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    });
    else local=ee._localPos;
    s.getPos()=(ROTI(info._TM,ee._JID.back())*local+CTRI(info._TM,ee._JID.back())).cast<scalar>();
    s.applyTrans();
    eeMesh.addMesh(s);
  }
  eeMesh.writeVTK(path,true);
}
template <typename T>
void C2GranularWrenchConstructor<T>::writeVTK(const std::string& path) const
{
  _env->getMesh().writeVTK(path,true);
}
template <typename T>
typename C2GranularWrenchConstructor<T>::Mat3T C2GranularWrenchConstructor<T>::rotateVecPlanar(const Vec3T d0,const Vec3T& d1,Mat3T dR0[3],Mat3T dR1[3])
{
  T sinVal=d0[2]*d1[0]-d0[0]*d1[2];
  T cosVal=d0[0]*d1[0]+d0[2]*d1[2];
  Mat3T R=Mat3T::Identity();
  R(0,0)=R(2,2)=cosVal;
  R(0,2)=sinVal;
  R(2,0)=-sinVal;
  if(dR0) {
    dR0[0].setZero();
    sinVal=-d1[2];
    cosVal=d1[0];
    dR0[0](0,0)=dR0[0](2,2)=cosVal;
    dR0[0](0,2)=sinVal;
    dR0[0](2,0)=-sinVal;

    dR0[1].setZero();

    dR0[2].setZero();
    sinVal=d1[0];
    cosVal=d1[2];
    dR0[2](0,0)=dR0[2](2,2)=cosVal;
    dR0[2](0,2)=sinVal;
    dR0[2](2,0)=-sinVal;
  }
  if(dR1) {
    dR1[0].setZero();
    sinVal=d0[2];
    cosVal=d0[0];
    dR1[0](0,0)=dR1[0](2,2)=cosVal;
    dR1[0](0,2)=sinVal;
    dR1[0](2,0)=-sinVal;

    dR1[1].setZero();

    dR1[2].setZero();
    sinVal=-d0[0];
    cosVal=d0[2];
    dR1[2](0,0)=dR1[2](2,2)=cosVal;
    dR1[2](0,2)=sinVal;
    dR1[2](2,0)=-sinVal;
  }
  return R;
}
template <typename T>
T C2GranularWrenchConstructor<T>::theta(const Vec3T& d0,const Vec3T& d1,Vec3T* dT0,Vec3T* dT1)
{
  T crossVal=d0[2]*d1[0]-d0[0]*d1[2];
  T ret=std::asin(crossVal);
  if(dT0) {
    T coef=1/std::sqrt(std::max<T>(1-crossVal*crossVal,0));
    dT0->setZero();
    (*dT0)[0]=-d1[2]*coef;
    (*dT0)[2]= d1[0]*coef;
  }
  if(dT1) {
    T coef=1/std::sqrt(std::max<T>(1-crossVal*crossVal,0));
    dT1->setZero();
    (*dT1)[0]= d0[2]*coef;
    (*dT1)[2]=-d0[0]*coef;
  }
  return ret;
}
template <typename T>
void C2GranularWrenchConstructor<T>::debugTheta(sizeType nrIter)
{
  DEFINE_NUMERIC_DELTA_T(T)
  for(sizeType it=0; it<nrIter; it++) {
    Vec3T a=Vec3T::Random();
    //a[1]=0;
    a/=std::sqrt(a.squaredNorm());
    Vec3T b=Vec3T::Random();
    //b[1]=0;
    b/=std::sqrt(b.squaredNorm());
    Vec3T dab=Vec3T::Random();
    T ret,ret2;

    Vec3T dT0,dT1;
    ret=theta(a,b,&dT0,&dT1);
    ret2=theta(a+dab*DELTA,b);
    DEBUG_GRADIENT("Dtheta0",dT0.dot(dab),dT0.dot(dab)-(ret2-ret)/DELTA)
    ret2=theta(a,b+dab*DELTA);
    DEBUG_GRADIENT("Dtheta1",dT1.dot(dab),dT1.dot(dab)-(ret2-ret)/DELTA)
  }
  for(sizeType it=0; it<nrIter; it++) {
    Vec3T a=Vec3T::Random();
    a[1]=0;
    a/=std::sqrt(a.squaredNorm());
    Vec3T b=Vec3T::Random();
    b[1]=0;
    b/=std::sqrt(b.squaredNorm());
    Vec3T dab=Vec3T::Random();
    Mat3T ret,ret2;

    Mat3T dT0[3],dT1[3];
    ret=rotateVecPlanar(a,b,dT0,dT1);
    DEBUG_GRADIENT("Rab",std::sqrt((ret*a).squaredNorm()),std::sqrt((ret*a-b).squaredNorm()))
    ret2=rotateVecPlanar(a+dab*DELTA,b);
    DEBUG_GRADIENT("DR0",std::sqrt((dT0[0]*dab[0]+dT0[1]*dab[1]+dT0[2]*dab[2]).squaredNorm()),
                   std::sqrt((dT0[0]*dab[0]+dT0[1]*dab[1]+dT0[2]*dab[2]-(ret2-ret)/DELTA).squaredNorm()))
    ret2=rotateVecPlanar(a,b+dab*DELTA);
    DEBUG_GRADIENT("DR1",std::sqrt((dT1[0]*dab[0]+dT1[1]*dab[1]+dT1[2]*dab[2]).squaredNorm()),
                   std::sqrt((dT1[0]*dab[0]+dT1[1]*dab[1]+dT1[2]*dab[2]-(ret2-ret)/DELTA).squaredNorm()))
  }
}

PRJ_BEGIN
template class C2GranularWrenchConstructor<double>;
#ifdef ALL_TYPES
template class C2GranularWrenchConstructor<__float128>;
template class C2GranularWrenchConstructor<mpfr::mpreal>;
#endif
PRJ_END
