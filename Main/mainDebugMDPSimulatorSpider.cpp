#include <Articulated/MDPSimulator.h>
#include <Utils/ArticulatedBodyPragma.h>
#include <Environment/EnvWrenchConstructor.h>
#include <Articulated/ArticulatedLoader.h>
#include <Utils/SpatialRotationUtil.h>
#include <Utils/Scalar.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

#ifdef ENVIRONMENT_SUPPORT
typedef double T;
#define debugOnly false
void debugGradientInfoSpider(const std::string& path,T dt,T dtw,MDP_SIMULATOR_MODE mode,bool floor=true,sizeType root=Joint::TRANS_3D|Joint::ROT_3D_XYZ,sizeType lmt=30)
{
  DECL_MAP_TYPES_T
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  scalar footLen=0.2f*sqrt(2.0f)+0.16f;
  ArticulatedLoader::createSpider(*(pt.RootElement()),root,0.25f,footLen,0.08f,D2R(10),D2R(lmt),D2R(lmt),true);
  ArticulatedBody body(*(pt.RootElement()));
  std::shared_ptr<C2EnvWrenchConstructor<T>> wrench(new C2EnvWrenchConstructor<T>(body,Vec4d(0,0,1,2)*20,Vec3T(0,0,-9.81f)));
  if(floor)
    for(sizeType i=0; i<body.nrJ(); i++)
      if(body.children(i,true).empty()) {
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
    sim1.debug(0.01f,10,0);
    sim1.debug(0.01f,10,1000);
    return;
  }
  sizeType N=body.nrDOF();
  Vec s=concat<Vec>(sim1._info._qM,Vec::Zero(N));
  s[3]=M_PI/2;
  Vec PCoef=Vec::Zero(N),DCoef=Vec::Zero(N);
  for(sizeType i=0; i<body.nrJ(); i++)
    if(!body.joint(i).isRoot(body)) {
      PCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(1);
      DCoef.segment(body.joint(i)._offDOF,body.joint(i).nrDOF()).setConstant(0.1);
    }
  std::cout << PCoef.transpose() << std::endl << DCoef.transpose() << std::endl;
  sim1.setState(s.segment(0,N),s.segment(N,N));
  sim1.writeVTKSeq(path,10,dtw,dt);
  wrench->writeVTK(path+"/MDPFloor.vtk");
}
#endif
int main()
{
#ifdef OPTIMIZER_SUPPORT
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debugGradientInfoSpider("Spider-MDP-PGM-0.025",0.025,0.025,NMDP_PGM);
#endif
  return 0;
}

