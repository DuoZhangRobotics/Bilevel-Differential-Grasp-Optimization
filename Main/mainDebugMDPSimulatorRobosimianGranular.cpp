#include <Articulated/MDPSimulator.h>
#include <Articulated/ArticulatedBodyPragma.h>
#include <TrajOpt/Environment/GranularWrenchConstructor.h>
#include <TrajOpt/Environment/Environment.h>
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
void debugGradientInfoRobosimian(const std::string& PDTargetPath,const std::string& path,T dt,T dtw,MDP_SIMULATOR_MODE mode,const std::string& pathTerrain="")
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
  std::shared_ptr<Environment<T>> env;
  if(!pathTerrain.empty()) {
    env.reset(new EnvironmentHeight<T>(pathTerrain,true,4));
    env->getMesh().writeVTK(std::experimental::filesystem::v1::path(pathTerrain).replace_extension(".vtk"),true);
  } else {
    env.reset(new EnvironmentExact<T>());
    if(debugOnly)
      std::dynamic_pointer_cast<EnvironmentExact<T>>(env)->createHills(5,5,[&](scalarD x,scalarD y) {
      return std::sin(x)*std::sin(y);
    },16);
    else std::dynamic_pointer_cast<EnvironmentExact<T>>(env)->createFloor(Vec4d(0,0,1,0.8)*20);
  }
  std::shared_ptr<C2GranularWrenchConstructor<T>> wrench(new C2GranularWrenchConstructor<T>(body,env));
  //ArticulatedBody bodyVisual=ArticulatedLoader::readURDF("Robosimian/Robosimian/robosimian_caesar_new_all_active.urdf",false,true);
  //bodyVisual.writeVTK(NEArticulatedGradientInfo<scalarD>(bodyVisual,Cold::Zero(bodyVisual.nrDOF())).getTrans(),"Robosimian/Robosimian-Rest.vtk",Joint::MESH);
  body.debugBase(std::experimental::filesystem::v1::path(path).parent_path().string()+"/Robosimian-DebugBase",1);

  Options ops;
  MDPSimulator<T> sim1(body,ops,Vec3T(0,0,-9.81f),mode);
  ops.setOptions<MultiPrecisionLQP<T>,bool>("callback",false);
  ops.setOptions<MultiPrecisionLQP<T>,bool>("highPrec",false);
  ops.setOptions<MultiPrecisionLQP<T>,T>("muFinal",1e-5f);
  sim1.reset(ops);
  sim1.setWrenchConstructor(wrench);
  if(debugOnly) {
    MDP<T>::debugWrenchConstructor(body,*wrench,10);
    wrench->debugTheta();
    sim1.debugNMDP(0.01f,10,1000);
    //sim1.debug(0.01f,10,0);
    //sim1.debug(0.01f,10,1000);
    return;
  }
  sizeType N=body.nrDOF();
  Vec s=concat<Vec>(sim1._info._qM,Vec::Zero(N));
  //s[1]=-1.5;
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
  sim1.begLog("log.dat");
  sim1.setState(s.segment(0,N),s.segment(N,N));
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
    std::cout << PCoef.transpose() << std::endl << DCoef.transpose() << std::endl;
    PD.reset(new PDTarget(PCoef.template cast<double>(),DCoef.template cast<double>(),s.template cast<double>()));
  }
  std::map<std::string,std::set<sizeType>> jointMask;
  for(sizeType i=0; i<body.nrJ(); i++)
    if(body.children(i,true).size()==4) {
      for(sizeType j=i; j>=0; j=body.joint(j)._parent)
        jointMask["TORSO"].insert(j);
      for(sizeType j:body.children(i,true))
        jointMask["LEG"+std::to_string(j)]=body.children(j);
      break;
    }
  sim1.writeVTKSeq(path,10,dtw,dt,PD,&jointMask);
  wrench->writeEndEffectorVTK(path+"/EndEffector.vtk",false);
  wrench->writeEndEffectorVTK(path+"/TorqueCenter.vtk",true);
  wrench->writeVTK(path+"/MDPFloor.vtk");
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  std::string path="/media/zherong/Extreme SSD/NMDP_Video/RobosimianMDPGranular";
  create("Robosimian");
  RandEngine::seed(0);
  //debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.001",0.001,0.1,FORWARD_RK1F,"terrain1.txt");
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.001",0.001,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.002",0.002,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.003",0.003,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.004",0.004,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.005",0.005,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.006",0.006,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.007",0.007,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.008",0.008,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.009",0.009,0.05,FORWARD_RK1F);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.010",0.010,0.05,FORWARD_RK1F);
  //debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-FORWARD_RK1F-0.050",0.050,0.05,FORWARD_RK1F);
  RandEngine::seed(0);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.005",0.005,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.010",0.010,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.015",0.015,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.020",0.020,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.025",0.025,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.030",0.030,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.035",0.035,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.040",0.040,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.045",0.045,0.05,NMDP_PGM);
  debugGradientInfoRobosimian(path+"/RobosimianWalkFast-PDTarget/PDTarget.dat",path+"/RobosimianGranularWalkFast-MDP-PGM-0.050",0.050,0.05,NMDP_PGM);
  return 0;
}

