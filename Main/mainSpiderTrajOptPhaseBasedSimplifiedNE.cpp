#include <Articulated/ArticulatedLoader.h>
#include <TrajOpt/Environment/Environment.h>
#include <CommonFile/Timing.h>
#include <TrajOpt/TrajOpt.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

//#define MPFR
#define STAGE_II
template <typename T>
void trajOpt(int argc,char** argv,bool debugDynamics)
{
  disableTiming();
  std::shared_ptr<EnvironmentExact<T>> env(new EnvironmentExact<T>);
  std::shared_ptr<SimplifiedDynamics> simplified(new SimplifiedDynamics);
  simplified->SerializableBase::read("spider/SimplifiedDynamics.dat");
  std::shared_ptr<ArticulatedBody> body=simplified->getBody();
  INFOV("Total mass=%f",body->totalMass())
  //body->scaleMass(0.15f/body->totalMass());
  env->createFloor(20,5,0);
  //init pose
  Cold pose=Cold::Zero(body->nrDOF());
  for(sizeType i=0; i<body->nrJ(); i++)
    if(body->joint(i)._typeJoint==Joint::TRANS_3D)
      pose.segment<3>(body->joint(i)._offDOF)=Vec3d(10,0,1);
    else if(body->joint(i)._typeJoint==Joint::ROT_3D_EXP)
      pose.segment<3>(body->joint(i)._offDOF)=Vec3d(M_PI/2,0,0);
  PBDArticulatedGradientInfo<T> poseInfo(*body,pose.template cast<T>());
  //Stage I
  TrajOpt<T,3> optStage1(env,simplified);
  optStage1.parsePtree(argc,argv);
  optStage1.formulation(TrajOpt<T,3>::STATIC);
  optStage1.solver(TrajOpt<T,3>::KNITRO);
  optStage1.usePBDOrNewtonEuler(false);
  optStage1.savePath("stage1Solution.dat");
  optStage1.buildProblem([&](sizeType jid,const EndEffectorBounds& ee) {
    bool largerX=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[0]>10;
    bool largerY=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[1]>0;
    return largerX!=largerY;
  },Vec3d::Zero());
  optStage1.debug(10,0,[&](const std::string& str) {
    return debugDynamics ? (str=="NESimplifiedDynamicsSequence" || str=="NESimplifiedStaticSequence") : false;
  });
  optStage1.setInitPose(pose.template cast<T>());
  optStage1.solve();
  optStage1.writeVTK("solutionPhaseBasedSimplifiedStage1");
  //Stage II
#ifdef STAGE_II
  TrajOpt<T,3> optStage2(env,simplified);
  optStage1.coefSmallForce(1E-6f);
  optStage2.parsePtree(argc,argv);
  optStage2.formulation(TrajOpt<T,3>::PHASE_BASED);
  optStage2.solver(TrajOpt<T,3>::KNITRO);
  optStage2.usePBDOrNewtonEuler(false);
  optStage2.useDirectOrCG(true);
  optStage2.maxIter(500);
  optStage2.callback(5);
  optStage2.derivCheck(0);
  optStage2.coefSmallForce(0);
  optStage2.savePath("stage2Solution.dat");
  optStage2.buildProblem([&](sizeType jid,const EndEffectorBounds& ee) {
    bool largerX=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[0]>10;
    bool largerY=(ROTI(poseInfo._TM,jid)*ee._localPos.template cast<T>()+CTRI(poseInfo._TM,jid))[1]>0;
    return largerX!=largerY;
  },Vec3d::Zero());
  //movement objective
  optStage2.setInit(optStage1);
  //optStage2.constrainPosition(Vec2i(0,1));
  optStage2.constrainRootVelocity(Coli::LinSpaced(197,3,199),Vec3d(0.5,0,0));
  optStage2.debug(10,0,[&](const std::string& str) {
    return debugDynamics ? (str=="NESimplifiedDynamicsSequence" || str=="NESimplifiedStaticSequence") : false;
  });
  optStage2.solve();
#endif
  //Visualization
#ifdef STAGE_II
  TrajOpt<T,3>& finalStage=optStage2;
#else
  TrajOpt<T,3>& finalStage=optStage1;
#endif
  finalStage.writeVTK("solutionPhaseBasedSimplifiedNE");
}
int main(int argc,char** argv)
{
#ifdef MPFR
  OmpSettings::getOmpSettingsNonConst().setNrThreads(1);
  mpfr_set_default_prec(1024U);
  trajOpt<mpfr::mpreal>(argc,argv,false);
#else
  trajOpt<scalarD>(argc,argv,true);
#endif
  return 0;
}
