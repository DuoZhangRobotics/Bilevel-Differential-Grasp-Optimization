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
void debugGradientInfoBird(const std::string& PDTargetPath,const std::string& path,T dt,T dtw,PBD_SIMULATOR_MODE mode,bool floor=true)
{
  DECL_MAP_TYPES_T
  scalarD rad=0.01f;
  ArticulatedBody body(ArticulatedLoader::createBird(0));
  ROT(body.joint(0)._trans)=expWGradV<scalarD,Vec3d>(Vec3d::UnitX()*M_PI/2);
  ArticulatedUtils(body).addBase(2,Vec3d::UnitY());
  ArticulatedUtils(body).scaleMass(100/ArticulatedUtils(body).totalMass());
  std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(body,Vec4d(0,0,1,0.8)*20,Vec3T(0,0,-9.81f),6,2.0f,100000000));
  wrench->writeVTK(path+"/Bird-MDPFloor.vtk");
  if(floor)
    for(sizeType i=0; i<body.nrJ(); i++)
      if(body.children(i,true).empty()) {
        if(body.joint(i).getMesh(Joint::MESH).getV().size()>100)
          continue;
        Vec3d zRange;
        EndEffectorBounds ee(i);
        zRange=Vec3d( 1, 1,0)*0.001f;
        SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0=0,zRange);
        ee._localPos-=Vec3d::UnitY()*(ee._phi0-rad);
        ee._phi0=rad;
        wrench->_externalForces.push_back(ee);

        zRange=Vec3d(-1,1,0)*0.001f;
        SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0=0,zRange);
        ee._localPos-=Vec3d::UnitY()*(ee._phi0-rad);
        ee._phi0=rad;
        wrench->_externalForces.push_back(ee);

        zRange=Vec3d( 1,-1,0)*0.001f;
        SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0=0,zRange);
        ee._localPos-=Vec3d::UnitY()*(ee._phi0-rad);
        ee._phi0=rad;
        wrench->_externalForces.push_back(ee);

        zRange=Vec3d(-1,-1,0)*0.001f;
        SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0=0,zRange);
        ee._localPos-=Vec3d::UnitY()*(ee._phi0-rad);
        ee._phi0=rad;
        wrench->_externalForces.push_back(ee);
      }
  body.debugBase(std::experimental::filesystem::v1::path(path).parent_path().string()+"/Bird-DebugBase",1);

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
        PCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(0);
        DCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(0);
      }
    std::cout << PCoef.transpose() << std::endl << DCoef.transpose() << std::endl;
    PD.reset(new PDTarget(PCoef.template cast<double>(),DCoef.template cast<double>(),s.template cast<double>()));
  }
  sim1.writeVTKSeq(path,10.0,dtw,dt,PD);
  wrench->writeEndEffectorVTK(path+"/Bird-EndEffector.vtk");
  wrench->writeVTK(path+"/Bird-MDPFloor.vtk");
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  create("Bird");
  RandEngine::seed(0);
  debugGradientInfoBird("Bird/BirdWalk-PDTarget/PDTarget.dat","Bird/BirdWalk-PBD-PGM-0.050",0.001,0.05,NPBD_PGM);
  RandEngine::seed(0);
  debugGradientInfoBird("Bird/BirdWalkFast-PDTarget/PDTarget.dat","Bird/BirdWalkFast-PBD-PGM-0.050",0.001,0.05,NPBD_PGM);
  return 0;
}

