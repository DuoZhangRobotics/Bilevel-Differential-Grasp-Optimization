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

  ASSERT_MSG(argn==3,"mainObjectSettle: [density] [sample density]")
  sizeType density=std::atoi(argc[1]);
  std::string pathObj(argc[2]);

  //load objects
  PhysicsRegistration<T> planner;
  if(exists("objects.dat")) {
    planner.SerializableBase::read("objects.dat");
  } else {
    std::vector<ObjMesh> objs(2);
    MakeMesh::makeSphere3D(objs[0],0.1f,16);
    MakeMesh::makeBox3D(objs[1],Vec3(0.2f,0.2f,0.1f));
    planner.reset(objs,1.0f/density,false,0.005f);
    planner.SerializableBase::write("objects.dat");
  }

  //load env
  PointCloudObject<T> env;
  env.SerializableBase::read(pathObj);
  env.writeVTK("env",1);

  sizeType itAll=0;
  recreate("objectSettle");
  Vec x0=Vec::Zero(planner.body().nrDOF());
  x0[6*0+2]=0.5f;
  x0[6*1+2]=0.2f;
  planner.writeVTK(x0,"objects",1);
  planner.setIndexModifier([&](sizeType,const Vec& x) {
    //if(x.size()>60)
    //  return;
    //write
    PBDArticulatedGradientInfo<T> info(planner.body(),x);
    planner.body().writeVTK(info._TM,"objectSettle/itConfig"+std::to_string(itAll)+".vtk",Joint::MESH);
    planner.writeContactVTK(x,"objectSettle/itContact"+std::to_string(itAll)+".vtk");
    itAll++;
    //add index
    std::set<PhysicsRegistration<T>::Penetration> pss;
    planner.getDeepestPenetration(pss);
    for(const PhysicsRegistration<T>::Penetration& p:pss) {
      std::tuple<sizeType,Vec3T,T> t=p._deepestPenetration;
      if(std::get<2>(t)<0 && !planner.existContactIndexMap(p._oid,p._oidOther,std::get<0>(t)))
        planner.addContactIndexMap(p._oid,p._oidOther,std::get<0>(t),0.7f,6);
    }
  });
  Options ops;
  PhysicsRegistrationParameter param(ops);
  param._coefPotential=1;
  param._useAugLag=false;
  planner.optimize(false,x0,env,NULL,param);
  return 0;
}
