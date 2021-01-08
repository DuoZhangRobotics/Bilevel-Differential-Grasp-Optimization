#ifdef ENVIRONMENT_SUPPORT
#include <Articulated/MDP.h>
#include <Environment/EnvWrenchConstructor.h>
#include <Articulated/ArticulatedLoader.h>
#include <Utils/SpatialRotationUtil.h>
#include <Utils/Scalar.h>

USE_PRJ_NAMESPACE

typedef double T;
#define debugMDPSwitch true
void debugGradientInfoChain(sizeType root,sizeType nrLink=10,sizeType GD=1)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  C0EnvWrenchConstructor<T>::debugRotateVec();
  ArticulatedLoader::createChain(*(pt.RootElement()),root,nrLink,2,0.2f,M_PI/4,M_PI/10.0f,0.1f,0.2f,GD);
  ArticulatedBody body(*(pt.RootElement()));
  MDP<T>::debug(body,10);
}
void debugGradientInfoSpider(sizeType root,sizeType lmt=30)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  C0EnvWrenchConstructor<T>::debugRotateVec();
  scalar footLen=0.2f*sqrt(2.0f)+0.16f;
  ArticulatedLoader::createSpider(*(pt.RootElement()),root,0.25f,footLen,0.08f,D2R(10),D2R(lmt),D2R(lmt),true);
  ArticulatedBody body(*(pt.RootElement()));
  MDP<T>::debug(body,10);
}
#endif
int main()
{
#ifdef ENVIRONMENT_SUPPORT
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debugGradientInfoChain(Joint::TRANS_3D|Joint::ROT_3D_EXP);
  debugGradientInfoChain(Joint::TRANS_2D|Joint::HINGE_JOINT);
  debugGradientInfoSpider(Joint::TRANS_3D|Joint::ROT_3D_XYZ);
#endif
  return 0;
}

