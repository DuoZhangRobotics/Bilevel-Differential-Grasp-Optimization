#include <Quasistatic/PhysicsRegistration.h>
#include <Quasistatic/PointCloudRegistrationEnergy.h>
#include <Quasistatic/ContactIndexConstraint.h>
#include <Quasistatic/GravityEnergy.h>
#include <CommonFile/CameraModel.h>
#include <CommonFile/MakeMesh.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

#define NR_DIR 3
typedef double T;
typedef PointCloudObject<T>::Vec Vec;
typedef PointCloudObject<T>::Vec3T Vec3T;
typedef PointCloudObject<T>::Mat3XT Mat3XT;
int main(int argn,char** argc)
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);

  ASSERT_MSG(argn==3,"mainObjectRegister: [density] [sample density]")
  sizeType density=std::atoi(argc[1]);
  std::string pathObj(argc[2]);
  Vec x0;

  //load objects
  PhysicsRegistration<T> planner;
  if(exists("objects.dat")) {
    planner.SerializableBase::read("objects.dat");
  } else {
    std::vector<ObjMesh> objs(2);
    MakeMesh::makeSphere3D(objs[0],0.1f,16);
    MakeMesh::makeBox3D(objs[1],Vec3(0.2f,0.2f,0.1f));
    planner.reset(objs,1.0f/density,false,0.005f,0.2f,true);
    planner.SerializableBase::write("objects.dat");
  }

  //load point cloud
  PointCloudObject<T> object;
  Mat4 prj=getProjection<scalar>(90,1,0.1,10);
  Mat4 m=lookAt<scalar>(Vec3(1,1,1),Vec3(0,0,0),Vec3(0,0,1));
  x0.setZero(planner.body().nrDOF());
  x0[6*0+0]+=0.1f;
  x0[6*1+0]+=0.1f;
  x0[6*0+1]+=0.1f;
  x0[6*1+1]+=0.1f;
  x0[6*0+2]=0.3f;
  x0[6*1+2]=0.1f;
  object.resetPointCloud(planner.body(),x0,m,prj,Vec2i(512,512));
  object.writeVTK("pointCloud",0.01f);

  //load env
  PointCloudObject<T> env;
  env.SerializableBase::read(pathObj);
  env.writeVTK("env",1);

  sizeType itAll=0;
  recreate("objectRegister");
  x0.setZero(planner.body().nrDOF());
  x0[6*0+2]=0.3f;
  x0[6*1+2]=0.1f;
  planner.writeVTK(x0,"objects",1);
  planner.setIndexModifier([&](sizeType,const Vec& x) {
    //if(x.size()>60)
    //  return;
    //write
    PBDArticulatedGradientInfo<T> info(planner.body(),x);
    planner.body().writeVTK(info._TM,"objectRegister/itConfig"+std::to_string(itAll)+".vtk",Joint::MESH);
    planner.writeContactVTK(x,"objectRegister/itContact"+std::to_string(itAll)+".vtk");
    itAll++;
    //add index
    std::set<PhysicsRegistration<T>::Penetration> pss;
    planner.getDeepestPenetration(pss);
    for(const PhysicsRegistration<T>::Penetration& p:pss) {
      std::tuple<sizeType,Vec3T,T> t=p._deepestPenetration;
      if(std::get<2>(t)<0 && !planner.existContactIndexMap(p._oid,p._oidOther,std::get<0>(t)))
        planner.addContactIndexMap(p._oid,p._oidOther,std::get<0>(t),0.7f,NR_DIR);
    }
  });
  Options ops;
  PhysicsRegistrationParameter param(ops);
  param._coefPotential=1;
  //param._useAugLag=false;
  param._g=Vec3d(0,-1,-9.81f);
  planner.optimize(false,x0,env,&object,param);
  return 0;
}
