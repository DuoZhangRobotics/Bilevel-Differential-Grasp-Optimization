#include <Articulated/MDPSimulator.h>
#include <Utils/ArticulatedBodyPragma.h>
#include <Environment/EnvWrenchConstructor.h>
#include <Articulated/ArticulatedLoader.h>
#include <Utils/SpatialRotationUtil.h>
#include <Utils/Scalar.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

#ifdef ENVIRONMENT_SUPPORT
#define SIN
typedef double T;
#define debugOnly false
void debugGradientInfoChain(const std::string& path,T dt,T dtw,MDP_SIMULATOR_MODE mode,bool floor=true,sizeType root=Joint::BALL_JOINT,sizeType nrLink=10,sizeType GD=1)
{
  DECL_MAP_TYPES_T
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),root,nrLink,2,0.2f,M_PI/4,0,0,0,GD);
  ArticulatedBody body(*(pt.RootElement()));

#if defined(CONE)
  std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(body,20,20,[](scalarD x,scalarD y) {
    return -(x*x+y*y)*0.05f-(debugOnly?0:3);
  },5,Vec3T(0,0,-9.81f)));
#elif defined(SIN)
  std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(body,20,20,[](scalarD x,scalarD y) {
    return -std::sin(x/2)*std::cos(y/2)*2-(debugOnly?0:3);
  },5,Vec3T(0,0,-9.81f)));
#endif
  if(floor)
    for(sizeType i=0; i<body.nrJ(); i++) {
      EndEffectorBounds ee(i);
      SimplifiedDynamics::detectEndEffector(body,i,ee._localPos,ee._phi0);
      wrench->_externalForces.push_back(ee);
    }

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
  sim1.writeVTKSeq(path,10,dtw,dt);
  wrench->writeVTK(path+"/floor.vtk");
}
#endif
int main()
{
#ifdef ENVIRONMENT_SUPPORT
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debugGradientInfoChain("chain-MDP-PGM-0.025",0.025,0.025,NMDP_PGM);
#endif
  return 0;
}

