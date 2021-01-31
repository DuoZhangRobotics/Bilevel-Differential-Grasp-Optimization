#include <Quasistatic/PhysicsRegistration.h>
#include <Quasistatic/PointCloudRegistrationEnergy.h>
#include <Quasistatic/ContactIndexConstraint.h>
#include <Quasistatic/GravityEnergy.h>
#include <CommonFile/MakeMesh.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

typedef double T;
typedef PointCloudObject<T>::Vec Vec;
typedef PointCloudObject<T>::Vec3T Vec3T;
typedef PointCloudObject<T>::Mat3XT Mat3XT;
int main(int argn,char** argc)
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);

  ASSERT_MSG(argn==3,"mainDebugRegistration: [density] [sample density]")
  sizeType density=std::atoi(argc[1]);
  std::string pathObj(argc[2]);

  //test objective
  PointCloudObject<T> object;
  object.SerializableBase::read(pathObj);
  object.writeVTK("pointCloud",1);

  //load objects
  std::vector<ObjMesh> objs(2);
  MakeMesh::makeSphere3D(objs[0],0.05f,16);
  MakeMesh::makeBox3D(objs[1],Vec3::Constant(0.05f));
  Options ops;
  PhysicsRegistrationParameter param(ops);
  PhysicsRegistration<T> planner;
  planner.reset(objs,1.0f/density);
  Vec x0=Vec::Random(planner.body().nrDOF());
  planner.writeVTK(x0,"objects",1);
  planner.debugPenetration(10,0.1f);
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    GravityEnergy<T>(obj,info,planner,object,Vec3T::Random(),10).debug(1,0,&x0,true);
  }
  for(sizeType pass=0; pass<2; pass++) {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    PointCloudRegistrationEnergy<T> pe(obj,info,planner,object,10);
    for(sizeType i=0; i<object.pss().cols(); i++)
      pe.addPointCloudMap(i,i>object.pss().cols()/2?0:1,pass==0?Vec3T(Vec3T::Constant(ScalarUtil<T>::scalar_nanq())):Vec3T(Vec3T::Random()));
    pe.debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    ContactIndexConstraint<T> cc(obj,info,planner,object,Vec3T::Random(),1.0f);
    for(sizeType i=0; i<100; i++)
      if(i<100/4) {
        const Mat3XT& pss=planner.pnss()[2].first;
        cc.addContactIndexMap(0,1,RandEngine::randR(0,pss.cols()-1),0.7f,6);
      } else if(i<100*2/4) {
        const Mat3XT& pss=planner.pnss()[4].first;
        cc.addContactIndexMap(1,0,RandEngine::randR(0,pss.cols()-1),0.7f,6);
      } else if(i<100*3/4) {
        const Mat3XT& pss=object.pss();
        cc.addContactIndexMap(-1,0,RandEngine::randR(0,pss.cols()-1),0.7f,6);
      } else {
        const Mat3XT& pss=object.pss();
        cc.addContactIndexMap(-1,1,RandEngine::randR(0,pss.cols()-1),0.7f,6);
      }
    cc.debug(1,0,&x0,true);
  }
  return 0;
}
