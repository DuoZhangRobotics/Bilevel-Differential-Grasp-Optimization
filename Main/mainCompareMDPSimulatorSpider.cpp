#include <Articulated/MDPSimulator.h>
#include <Articulated/ArticulatedBodyPragma.h>
#include <TrajOpt/Environment/EnvWrenchConstructor.h>
#include <Articulated/SpatialRotationUtil.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Articulated/PDTarget.h>
#include <Utils/Scalar.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

typedef double T;
#define debugOnly false
#define VISUAL_MESH true
void debugGradientInfoSpider(const std::string& PDTargetPath,const std::string& path,T dt,T dtw,MDP_SIMULATOR_MODE mode,bool floor=true)
{
  DECL_MAP_TYPES_T
  ArticulatedBody body(ArticulatedLoader::createSpider(0));
  ROT(body.joint(0)._trans)=expWGradV<scalarD,Vec3d>(Vec3d::UnitX()*M_PI/2);
  ArticulatedUtils(body).addBase(2,Vec3d::UnitY());
  ArticulatedUtils(body).scaleMass(100/ArticulatedUtils(body).totalMass());
  std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(body,Vec4d(0,0,1,0.43)*20,Vec3T(0,0,-9.81f),6,1.5f,1000000));
  wrench->writeVTK(path+"/Spider-MDPFloor.vtk");
  if(floor)
    for(sizeType i=0; i<body.nrJ(); i++)
      if(body.children(i,true).empty()) {
        EndEffectorBounds ee(i);
        SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0=0);
        wrench->_externalForces.push_back(ee);
      }
  body.debugBase(std::experimental::filesystem::v1::path(path).parent_path().string()+"/Spider-DebugBase",1);

  Options ops;
  MDPSimulator<T> sim1(body,ops,Vec3T(0,0,-9.81f),mode);
  ops.setOptions<MultiPrecisionLQP<T>,bool>("callback",false);
  ops.setOptions<MultiPrecisionLQP<T>,bool>("highPrec",false);
  sim1.reset(ops);
  sim1.setWrenchConstructor(wrench);
  if(debugOnly) {
    MDP<T>::debugWrenchConstructor(body,*wrench,10);
    sim1.debugNMDP(0.01f,10,1000);
    //sim1.debug(0.01f,10,0);
    //sim1.debug(0.01f,10,1000);
    return;
  }
  sizeType N=body.nrDOF();
  Vec s=Vec::Zero(N*2);
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=body.joint(i);
    std::cout << "Joint" << i << ": name=" << J._name << " parent=" << J._parent;
    if(J._typeJoint==Joint::HINGE_JOINT)
      std::cout << " lower=" << J._limits(0,0) << " upper=" << J._limits(1,0);
    std::cout << " type=" << Joint::typeToString((Joint::JOINT_TYPE)J._typeJoint) << std::endl;
  }
  sim1.setState(s.segment(0,N),s.segment(0,N));
  std::shared_ptr<PDTarget> PD;
  if(exists(PDTargetPath)) {
    PD.reset(new PDTarget);
    PD->SerializableBase::read(PDTargetPath);
    INFOV("Using PDTarget: %s!",PDTargetPath.c_str())
  } else {
    Vec PCoef=Vec::Zero(N),DCoef=Vec::Zero(N);
    for(sizeType i=0; i<body.nrJ(); i++)
      if(!body.joint(i).isRoot(body)) {
        PCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(10);
        DCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(0.1);
      }
    std::cout << PCoef.transpose() << std::endl << DCoef.transpose() << std::endl;
    PD.reset(new PDTarget(PCoef.template cast<double>(),DCoef.template cast<double>(),s.template cast<double>()));
  }
  sim1.writeVTKSeq(path,20,dtw,dt,PD);
  wrench->writeEndEffectorVTK(path+"/Spider-EndEffector.vtk");
  wrench->writeVTK(path+"/Spider-MDPFloor.vtk");
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  create("Spider");
  RandEngine::seed(0);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.005",0.005,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.010",0.010,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.015",0.015,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.020",0.020,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.025",0.025,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.030",0.030,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.035",0.035,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.040",0.040,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.045",0.045,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalk-PDTarget/PDTarget.dat","Spider/SpiderWalk-MDP-FORWARD_RK1F-0.050",0.050,0.05,FORWARD_RK1F);
  RandEngine::seed(0);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.005",0.005,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.010",0.010,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.015",0.015,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.020",0.020,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.025",0.025,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.030",0.030,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.035",0.035,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.040",0.040,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.045",0.045,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderWalkFast-PDTarget/PDTarget.dat","Spider/SpiderWalkFast-MDP-FORWARD_RK1F-0.050",0.050,0.05,FORWARD_RK1F);
  RandEngine::seed(0);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.005",0.005,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.010",0.010,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.015",0.015,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.020",0.020,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.025",0.025,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.030",0.030,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.035",0.035,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.040",0.040,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.045",0.045,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJump-PDTarget/PDTarget.dat","Spider/SpiderJump-MDP-FORWARD_RK1F-0.050",0.050,0.05,FORWARD_RK1F);
  RandEngine::seed(0);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.005",0.005,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.010",0.010,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.015",0.015,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.020",0.020,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.025",0.025,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.030",0.030,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.035",0.035,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.040",0.040,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.045",0.045,10000,FORWARD_RK1F);
  debugGradientInfoSpider("Spider/SpiderJumpHigh-PDTarget/PDTarget.dat","Spider/SpiderJumpHigh-MDP-FORWARD_RK1F-0.050",0.050,0.05,FORWARD_RK1F);
  return 0;
}

