#include <Articulated/PBDSimulator.h>
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
void debugGradientInfoRobosimian(const std::string& PDTargetPath,const std::string& path,T dt,T dtw,PBD_SIMULATOR_MODE mode,bool floor=true)
{
  DECL_MAP_TYPES_T
  ArticulatedBody body=ArticulatedLoader::readURDF("data/Robosimian/robosimian_caesar_new_all_active.urdf",false,VISUAL_MESH);
  body.addBase(2,Vec3d::UnitY());
  std::vector<std::string> unusedJoints;
  unusedJoints.push_back("limb1_link0");
  unusedJoints.push_back("limb1_link1");
  unusedJoints.push_back("limb1_link3");
  unusedJoints.push_back("limb1_link5");
  unusedJoints.push_back("limb1_link7");
  unusedJoints.push_back("limb2_link0");
  unusedJoints.push_back("limb2_link1");
  unusedJoints.push_back("limb2_link3");
  unusedJoints.push_back("limb2_link5");
  unusedJoints.push_back("limb2_link7");
  unusedJoints.push_back("limb3_link0");
  unusedJoints.push_back("limb3_link1");
  unusedJoints.push_back("limb3_link3");
  unusedJoints.push_back("limb3_link5");
  unusedJoints.push_back("limb3_link7");
  unusedJoints.push_back("limb4_link0");
  unusedJoints.push_back("limb4_link1");
  unusedJoints.push_back("limb4_link3");
  unusedJoints.push_back("limb4_link5");
  unusedJoints.push_back("limb4_link7");
  body.eliminateJoint(unusedJoints,Cold::Zero(body.nrDOF()),10);
  std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(body,Vec4d(0,0,1,0.8)*20,Vec3T(0,0,-9.81f),6,0.7f,1000000,1E-4f));
  if(floor)
    for(sizeType i=0; i<body.nrJ(); i++)
      if(body.children(i,true).empty()) {
        Vec3d zRange;
        EndEffectorBounds ee(i);
        if(beginsWith(body.joint(i)._name,"limb1"))
          zRange=Vec3d(-0.1f,0,0);
        else if(beginsWith(body.joint(i)._name,"limb2"))
          zRange=Vec3d( 0.1f,0,0);
        else if(beginsWith(body.joint(i)._name,"limb3"))
          zRange=Vec3d( 0.1f,0,0);
        else if(beginsWith(body.joint(i)._name,"limb4"))
          zRange=Vec3d(-0.1f,0,0);
        SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0,zRange);
        wrench->_externalForces.push_back(ee);
      }
  body.debugBase(std::experimental::filesystem::v1::path(path).parent_path().string()+"/Robosimian-DebugBase",1);
  std::cout << ArticulatedUtils(body).totalMass() << std::endl;
  //exit(-1);

  Options ops;
  PBDSimulator<T> sim1(body,ops,Vec3T(0,0,-9.81f),mode);
  ops.setOptions<PBDSimulator<T>,bool>("callback",false);
  sim1.reset(ops);
  sim1.setWrenchConstructor(wrench);
  if(debugOnly) {
    sim1.debug(0.01f,10,1000);
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
    if(J._name=="limb1_link2+limb1_link3")
      s[J._offDOF]= M_PI*0.9f/2;
    if(J._name=="limb1_link4+limb1_link5")
      s[J._offDOF]=-M_PI/2;
    if(J._name=="limb1_link6+limb1_link7")
      s[J._offDOF]= M_PI*1.1f/2;

    if(J._name=="limb2_link2+limb2_link3")
      s[J._offDOF]=-M_PI*0.9f/2;
    if(J._name=="limb2_link4+limb2_link5")
      s[J._offDOF]= M_PI/2;
    if(J._name=="limb2_link6+limb2_link7")
      s[J._offDOF]=-M_PI*1.1f/2;

    if(J._name=="limb3_link2+limb3_link3")
      s[J._offDOF]= M_PI*0.9f/2;
    if(J._name=="limb3_link4+limb3_link5")
      s[J._offDOF]=-M_PI/2;
    if(J._name=="limb3_link6+limb3_link7")
      s[J._offDOF]= M_PI*1.1f/2;

    if(J._name=="limb4_link2+limb4_link3")
      s[J._offDOF]=-M_PI*0.9f/2;
    if(J._name=="limb4_link4+limb4_link5")
      s[J._offDOF]= M_PI/2;
    if(J._name=="limb4_link6+limb4_link7")
      s[J._offDOF]=-M_PI*1.1f/2;
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
        PCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(1000);
        DCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(1.0);
      }
    //std::cout << PCoef.transpose() << std::endl << DCoef.transpose() << std::endl;
    PD.reset(new PDTarget(PCoef.template cast<double>(),DCoef.template cast<double>(),s.template cast<double>()));
  }
  sim1.writeVTKSeq(path,10,dtw,dt,PD);
  wrench->writeEndEffectorVTK(path+"/Robosimian-EndEffector.vtk");
  wrench->writeVTK(path+"/Robosimian-MDPFloor.vtk");
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  create("Robosimian");
  /*RandEngine::seed(0);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.005",0.005,0.05,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.010",0.010,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.015",0.015,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.020",0.020,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.025",0.025,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.030",0.030,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.035",0.035,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.040",0.040,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.045",0.045,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalk-PDTarget/PDTarget.dat","Robosimian/RobosimianWalk-PBD-ZOPGM-0.050",0.050,10000,NPBD_ZOPGM);
  RandEngine::seed(0);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.005",0.005,0.05,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.010",0.010,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.015",0.015,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.020",0.020,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.025",0.025,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.030",0.030,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.035",0.035,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.040",0.040,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.045",0.045,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianWalkFast-PDTarget/PDTarget.dat","Robosimian/RobosimianWalkFast-PBD-ZOPGM-0.050",0.050,10000,NPBD_ZOPGM);
  RandEngine::seed(0);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.005",0.005,0.05,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.010",0.010,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.015",0.015,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.020",0.020,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.025",0.025,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.030",0.030,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.035",0.035,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.040",0.040,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.045",0.045,10000,NPBD_ZOPGM);
  debugGradientInfoRobosimian("Robosimian/RobosimianJump-PDTarget/PDTarget.dat","Robosimian/RobosimianJump-PBD-ZOPGM-0.050",0.050,10000,NPBD_ZOPGM);*/
  return 0;
}

