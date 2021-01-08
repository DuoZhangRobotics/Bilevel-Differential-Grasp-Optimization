#include <Articulated/ArticulatedLoader.h>
#include <TrajOpt/Environment/Environment.h>
#include <CommonFile/Timing.h>
#include <TrajOpt/TrajOpt.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

#define STAGE_II
//#define PLANAR_ONLY
template <typename T>
void trajOpt(int argc,char** argv,bool debugDynamics)
{
  disableTiming();
  std::shared_ptr<EnvironmentExact<T>> env(new EnvironmentExact<T>);
  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody(ArticulatedLoader::readURDF("Robosimian/robosimian_caesar_new_all_active.urdf",false,false)));
  INFOV("Total mass=%f",body->totalMass())
  body->scaleMass(1/body->totalMass());
#ifdef PLANAR_ONLY
  body->addBase(2,Vec3d::UnitY());
#else
  body->addBase(3,Vec3d::Zero());
#endif
  std::vector<std::string> unusedJoints;
#ifdef PLANAR_ONLY
  unusedJoints.push_back("limb1_link0");
  unusedJoints.push_back("limb1_link1");
  unusedJoints.push_back("limb2_link0");
  unusedJoints.push_back("limb2_link1");
  unusedJoints.push_back("limb3_link0");
  unusedJoints.push_back("limb3_link1");
  unusedJoints.push_back("limb4_link0");
  unusedJoints.push_back("limb4_link1");
#endif
  unusedJoints.push_back("limb1_link3");
  unusedJoints.push_back("limb1_link5");
  unusedJoints.push_back("limb1_link7");
  unusedJoints.push_back("limb2_link3");
  unusedJoints.push_back("limb2_link5");
  unusedJoints.push_back("limb2_link7");
  unusedJoints.push_back("limb3_link3");
  unusedJoints.push_back("limb3_link5");
  unusedJoints.push_back("limb3_link7");
  unusedJoints.push_back("limb4_link3");
  unusedJoints.push_back("limb4_link5");
  unusedJoints.push_back("limb4_link7");
  body->eliminateJoint(unusedJoints,Cold::Zero(body->nrDOF()),10);
  env->createFloor(20,5,0);
  //init pose
  Cold pose=Cold::Zero(body->nrDOF());
  for(sizeType i=0; i<body->nrJ(); i++) {
    const Joint& J=body->joint(i);
    if(J._typeJoint==Joint::TRANS_3D)
      pose.segment<3>(J._offDOF)=Vec3d(10,0,0.75);
    if(J._typeJoint==Joint::TRANS_2D)
      pose.segment<2>(J._offDOF)=Vec2d(10,-0.75);
    if(J._name=="limb1_link2+limb1_link3")
      pose[J._offDOF]= M_PI*0.9f/2;
    if(J._name=="limb1_link4+limb1_link5")
      pose[J._offDOF]=-M_PI/2;
    if(J._name=="limb1_link6+limb1_link7")
      pose[J._offDOF]= M_PI*1.1f/2;

    if(J._name=="limb2_link2+limb2_link3")
      pose[J._offDOF]=-M_PI*0.9f/2;
    if(J._name=="limb2_link4+limb2_link5")
      pose[J._offDOF]= M_PI/2;
    if(J._name=="limb2_link6+limb2_link7")
      pose[J._offDOF]=-M_PI*1.1f/2;

    if(J._name=="limb3_link2+limb3_link3")
      pose[J._offDOF]= M_PI*0.9f/2;
    if(J._name=="limb3_link4+limb3_link5")
      pose[J._offDOF]=-M_PI/2;
    if(J._name=="limb3_link6+limb3_link7")
      pose[J._offDOF]= M_PI*1.1f/2;

    if(J._name=="limb4_link2+limb4_link3")
      pose[J._offDOF]=-M_PI*0.9f/2;
    if(J._name=="limb4_link4+limb4_link5")
      pose[J._offDOF]= M_PI/2;
    if(J._name=="limb4_link6+limb4_link7")
      pose[J._offDOF]=-M_PI*1.1f/2;
  }
  PBDArticulatedGradientInfo<T> poseInfo(*body,pose.template cast<T>());
  //Stage I
  TrajOpt<T,3> optStage1(env,body);
  optStage1.parsePtree(argc,argv);
  optStage1.formulation(TrajOpt<T,3>::STATIC_PBTO_LQP);
  optStage1.solver(TrajOpt<T,3>::KNITRO);
  optStage1.derivCheck(0);
  optStage1.PBTOPose(true);
  optStage1.savePath("stage1Solution.dat");
  optStage1.buildProblem([&](sizeType jid,EndEffectorBounds& ee) {
    Vec3d zRange;
    if(beginsWith(body->joint(jid)._name,"limb1"))
      zRange=Vec3d(-0.1f,0,0);
    else if(beginsWith(body->joint(jid)._name,"limb2"))
      zRange=Vec3d( 0.1f,0,0);
    else if(beginsWith(body->joint(jid)._name,"limb3"))
      zRange=Vec3d( 0.1f,0,0);
    else if(beginsWith(body->joint(jid)._name,"limb4"))
      zRange=Vec3d(-0.1f,0,0);
    else {
      ASSERT_MSG(false,"Unknown joint!")
    }
    ee=EndEffectorBounds(jid);
    SimplifiedDynamics::detectEndEffector(*body,jid,ee._localPos,ee._phi0,zRange);
    bool largerX=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[0]>10;
    bool largerY=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[1]>0;
    return largerX!=largerY;
  });
  optStage1.debug(10,0,[&](const std::string& str) {
    return debugDynamics ? (str=="PBDStaticSequence" || str=="PBDDynamicsSequence") : false;
  });
  optStage1.setInitPose(pose.template cast<T>());
  optStage1.solve();
  optStage1.writeVTK("solutionPBTOStage1");
  //Stage II
#ifdef STAGE_II
  TrajOpt<T,3> optStage2(env,body);
  optStage2.parsePtree(argc,argv);
  optStage2.formulation(TrajOpt<T,3>::PBTO_LQP);
  optStage2.solver(TrajOpt<T,3>::KNITRO);
  optStage2.useDirectOrCG(true);
  optStage2.betaPBTO(2000);
  optStage2.maxIter(500);
  optStage2.mu(0.2);
  optStage2.subStage(5);
  optStage2.callback(5);
  optStage2.derivCheck(0);
  optStage2.coefSmallForce(0);
  optStage2.coefShuffleAvoidance(1);
  optStage2.period(40);
  optStage2.coefPeriod(1);
  optStage2.QPMollifier(true);
  optStage2.PBTOPose(true);
  optStage2.savePath("stage2Solution.dat");
  optStage2.buildProblem([&](sizeType jid,EndEffectorBounds& ee) {
    Vec3d zRange;
    if(beginsWith(body->joint(jid)._name,"limb1"))
      zRange=Vec3d(-0.1f,0,0);
    else if(beginsWith(body->joint(jid)._name,"limb2"))
      zRange=Vec3d( 0.1f,0,0);
    else if(beginsWith(body->joint(jid)._name,"limb3"))
      zRange=Vec3d( 0.1f,0,0);
    else if(beginsWith(body->joint(jid)._name,"limb4"))
      zRange=Vec3d(-0.1f,0,0);
    else {
      ASSERT_MSG(false,"Unknown joint!")
    }
    ee=EndEffectorBounds(jid);
    SimplifiedDynamics::detectEndEffector(*body,jid,ee._localPos,ee._phi0,zRange);
    bool largerX=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[0]>10;
    bool largerY=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[1]>0;
    return largerX!=largerY;
  });
  //movement objective
  optStage2.setInit(optStage1);
  //optStage2.constrainPosition(Vec2i(0,1));
  optStage2.constrainRootVelocity(Coli::LinSpaced(197,3,199),Vec3d(0.3,0,0),1);
  optStage2.debug(10,0,[&](const std::string& str) {
    return debugDynamics ? (str=="PBDStaticSequence" || str=="PBDDynamicsSequence" || endsWith(str,"ShuffleAvoidanceObj")) : false;
  });
  optStage2.solve();
#endif
  //Visualization
#ifdef STAGE_II
  TrajOpt<T,3>& finalStage=optStage2;
#else
  TrajOpt<T,3>& finalStage=optStage1;
#endif
  finalStage.writeVTK("solutionPBTO");
}
int main(int argc,char** argv)
{
#ifdef MPFR
  OmpSettings::getOmpSettingsNonConst().setNrThreads(1);
  mpfr_set_default_prec(1024U);
  trajOpt<mpfr::mpreal>(argc,argv,false);
#else
  trajOpt<scalarD>(argc,argv,false);
#endif
  return 0;
}
