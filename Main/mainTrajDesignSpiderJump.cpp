#include <Articulated/PBDSimulator.h>
#include <Articulated/SimplifiedDynamics.h>
#include <TrajOpt/Environment/EnvWrenchConstructor.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/RotationUtil.h>
#include <Articulated/PDTarget.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

typedef double T;
DECL_MAP_TYPES_T
Vec3T getEEPos(const ArticulatedBody& body,std::shared_ptr<C2EnvWrenchConstructor<T>> wrench,std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>>* pos,bool* even,const Vec& DOF)
{
  Vec3T ctr=Vec3T::Zero();
  PBDArticulatedGradientInfo<T> info;
  info.reset(body,DOF);
  if(pos)
    pos->resize(wrench->_externalForces.size());
  for(sizeType i=0; i<(sizeType)wrench->_externalForces.size(); i++) {
    const EndEffectorBounds& ee=wrench->_externalForces[i];
    Vec3T EEPos=ROTI(info._TM,ee._JID[0])*ee._localPos.template cast<T>()+CTRI(info._TM,ee._JID[0]);
    if(pos)
      (*pos)[i]=EEPos;
    ctr+=EEPos;
  }
  ctr/=(sizeType)wrench->_externalForces.size();
  for(sizeType i=0; i<(sizeType)wrench->_externalForces.size(); i++) {
    const EndEffectorBounds& ee=wrench->_externalForces[i];
    Vec3T EEPos=ROTI(info._TM,ee._JID[0])*ee._localPos.template cast<T>()+CTRI(info._TM,ee._JID[0]);
    if(even)
      even[i]=((EEPos-ctr)[0]>0)^((EEPos-ctr)[1]>0);
  }
  return ctr;
}
void trajDesignSpider(const std::string& path,sizeType step=2,T dt=1.0f,T dt2=0.2f,Vec3T dirZ=Vec3T(0,0,0.4f),T dirX=0.0f)
{
  ArticulatedBody body(ArticulatedLoader::createSpider(0));
  ROT(body.joint(0)._trans)=expWGradV<scalarD,Vec3d>(Vec3d::UnitX()*M_PI/2);
  ArticulatedUtils(body).addBase(2,Vec3d::UnitY());
  std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(body,Vec4d(0,0,1,2)*20,Vec3T(0,0,-9.81f),6,0.7f,100000));
  for(sizeType i=0; i<body.nrJ(); i++)
    if(body.children(i,true).empty()) {
      EndEffectorBounds ee(i);
      SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0=0);
      wrench->_externalForces.push_back(ee);
    }
  wrench->writeEndEffectorVTK(std::experimental::filesystem::v1::path(path).parent_path().string()+"/EndEffector.vtk");
  body.debugBase(std::experimental::filesystem::v1::path(path).parent_path().string()+"/Spider-DebugBase",1);
  //if(exists(path+"-PDTarget/PDTarget.dat")) {
  //  PDTarget PDT;
  //  PDT.SerializableBase::read(path+"-PDTarget/PDTarget.dat");
  //  PDT.writeVTKSeq(body,path+"-PDTargetReplay",0.05f);
  //  return;
  //}

  //init
  Vec s=Vec::Zero(body.nrDOF());
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    std::cout << "Joint" << i << ": name=" << J._name << " parent=" << J._parent << std::endl;
  }

  //ctr
  bool even[4];
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> pos;
  Vec3T ctr=getEEPos(body,wrench,NULL,even,s);

  //walk
  T t=dt;
  Vec init=s;
  Mat3T H=Vec3T(1,1,1).asDiagonal();
  std::vector<std::tuple<scalarD,Vec,Vec>> PD;
  PD.push_back(std::make_tuple(0,s,Vec::Zero(body.nrDOF())));
  for(sizeType i=0; i<step; i++) {
    Vec s=std::get<1>(PD.back());
    getEEPos(body,wrench,&pos,NULL,s);
    //stage1
    for(sizeType ie=0; ie<4; ie++)
      pos[ie]+=dirZ/2;
    s.segment(3,s.size()-3)=init.segment(3,init.size()-3);
    s=SimplifiedDynamics::inverseKinematics(body,wrench->_externalForces,pos,H,&s,false,false);
    PD.push_back(std::make_tuple(t,s,Vec::Zero(body.nrDOF())));
    t+=dt2;
    if(s.size()==0) {
      INFO("Failed!")
      return;
    }
    //stage2
    for(sizeType ie=0; ie<4; ie++) {
      pos[ie]-=dirZ;
      if(pos[ie].x()<ctr.x())
        pos[ie].x()+=dirX;
      else pos[ie].x()-=dirX;
    }
    s.segment(3,s.size()-3)=init.segment(3,init.size()-3);
    s=SimplifiedDynamics::inverseKinematics(body,wrench->_externalForces,pos,H,&s,false,false);
    PD.push_back(std::make_tuple(t,s,Vec::Zero(body.nrDOF())));
    t+=dt;
    if(s.size()==0) {
      INFO("Failed!")
      return;
    }
    //stage3
    for(sizeType ie=0; ie<4; ie++) {
      pos[ie]+=dirZ/2;
      if(pos[ie].x()<ctr.x())
        pos[ie].x()-=dirX;
      else pos[ie].x()+=dirX;
    }
    s.segment(3,s.size()-3)=init.segment(3,init.size()-3);
    s=SimplifiedDynamics::inverseKinematics(body,wrench->_externalForces,pos,H,&s,false,false);
    PD.push_back(std::make_tuple(t,s,Vec::Zero(body.nrDOF())));
    t+=dt;
    if(s.size()==0) {
      INFO("Failed!")
      return;
    }
  }

  //compile
  {
    Vec PCoef=Vec::Zero(body.nrDOF()),DCoef=Vec::Zero(body.nrDOF());
    for(sizeType i=0; i<body.nrJ(); i++)
      if(!body.joint(i).isRoot(body)) {
        PCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(1000);
        DCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(1.0);
      }
    PDTarget PDT(PCoef,DCoef,PD);
    PDT.writeVTKSeq(body,path+"-PDTarget",0.05f);
    PDT.SerializableBase::write(path+"-PDTarget/PDTarget.dat");
  }
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  create("Spider");
  RandEngine::seed(0);
  trajDesignSpider("Spider/SpiderJump",10);
  trajDesignSpider("Spider/SpiderJumpHigh",20,1.0f,0.2f,Vec3T(0,0,0.6f),0.2f);
  return 0;
}

