#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Articulated/ArticulatedLoader.h>
#include <Utils/Scalar.h>

USE_PRJ_NAMESPACE

typedef double T;
void debugGradientInfoChain(sizeType root,sizeType nrLink=10,sizeType GD=1)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  ArticulatedLoader::createChain(*(pt.RootElement()),root,nrLink,2,0.2f,M_PI/4,M_PI/10.0f,0.1f,0.2f,GD);
  ArticulatedBody body(*(pt.RootElement()));
  PBDArticulatedGradientInfo<T>::debug(body);
}
void debugGradientInfoSpider(sizeType root,sizeType lmt=30)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  scalar footLen=0.2f*sqrt(2.0f)+0.16f;
  ArticulatedLoader::createSpider(*(pt.RootElement()),root,0.25f,footLen,0.08f,D2R(10),D2R(lmt),D2R(lmt),true);
  ArticulatedBody body(*(pt.RootElement()));
  PBDArticulatedGradientInfo<T>::debug(body);
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debugGradientInfoChain(Joint::TRANS_3D|Joint::ROT_3D_EXP);
  debugGradientInfoChain(Joint::TRANS_2D|Joint::HINGE_JOINT);
  debugGradientInfoSpider(Joint::TRANS_3D|Joint::ROT_3D_XYZ);
  return 0;
}

