#include <Environment/ObjMeshGeomCellExact.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/Environment.h>
#include <Articulated/ConvexHull.h>
#include <CommonFile/MakeMesh.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

#ifdef ENVIRONMENT_SUPPORT
sizeType res=10;
#define FORCE_CREATE false
template <typename T>
void debugEnvExact()
{
  std::shared_ptr<EnvironmentExact<T>> env;
  if(exists("envExact.dat") && !FORCE_CREATE) {
    env.reset(new EnvironmentExact<T>());
    env->SerializableBase::read("envExact.dat");
  } else {
    env.reset(new EnvironmentExact<T>());
    env->createStair(5,5,5,2.5,0.1,10);
    env->getMesh().write("Stair.obj");
    env->createHills(5,5,[&](scalarD x,scalarD y) {
      return std::sin(x)*std::sin(y);
    },16);
    env->getMesh().write("Hills.obj");
    env->SerializableBase::write("envExact.dat");
  }
  env->debugVTK("envExact",res);
  env->debug();
}
template <typename T>
void debugEnvCubic()
{
  std::shared_ptr<EnvironmentCubic<T>> env;
  if(exists("envCubic.dat") && !FORCE_CREATE) {
    env.reset(new EnvironmentCubic<T>(0.1f));
    env->SerializableBase::read("envCubic.dat");
  } else {
    env.reset(new EnvironmentCubic<T>(0.1f));
    //env->createStair(5,5,5,2.5,0.1,10);
    env->createHills(5,5,[&](scalarD x,scalarD y) {
      return std::sin(x)*std::sin(y);
    },16);
    env->SerializableBase::write("envCubic.dat");
  }
  env->debugVTK("envCubic");
  env->debug();
}
template <typename T>
void debugEnvHeight()
{
  std::shared_ptr<EnvironmentHeight<T>> env;
  if(exists("envHeight.dat") && !FORCE_CREATE) {
    env.reset(new EnvironmentHeight<T>(0.1f));
    env->SerializableBase::read("envHeight.dat");
  } else {
    env.reset(new EnvironmentHeight<T>(0.1f));
    //env->createStair(5,5,5,2.5,0.1,10);
    env->createHills(5,5,[&](scalarD x,scalarD y) {
      return std::sin(x)*std::sin(y);
    },16);
    env->SerializableBase::write("envHeight.dat");
  }
  env->debugVTK("envHeight");
  env->debug();
}
template <typename T>
void debugConvexHull()
{
  ObjMesh m;
  MakeMesh::makeSphere3D(m,1,8);
  //MakeMesh::makeBox3D(m,Vec3::Ones());
  m=makeConvex(m);
  ObjMeshGeomCell cell(Mat4::Identity(),m,0,true);

  ConvexHullExact cExact(cell);
  ObjMeshGeomCellExact mExact(cell);
  cExact.writePointDistVTK("dist.vtk");
  cExact.writeHistoryDistVTK("history");
  cExact.getMesh().writeVTK("mesh.vtk",true);
  cExact.compare(mExact);
}
#endif
int main()
{
#ifdef ENVIRONMENT_SUPPORT
  mpfr_set_default_prec(1024U);
  //debugEnvExact<double>();
  //debugEnvCubic<double>();
  //debugEnvHeight<double>();
  debugConvexHull<double>();
#endif
  return 0;
}
