#include <Grasp/GraspPlanner.h>
#include <Grasp/CentroidClosednessEnergy.h>
#include <Grasp/ObjectClosednessEnergy.h>
#include <Grasp/LogBarrierObjEnergy.h>
#include <Grasp/LogBarrierSelfEnergy.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

typedef double T;
typedef GraspQualityMetric<T>::Vec Vec;
typedef GraspQualityMetric<T>::Vec3T Vec3T;
int main(int argn,char** argc)
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);

  ASSERT_MSG(argn==4,"mainGraspPlan: [urdf path] [sample density] [obj path]")
  std::string path(argc[1]);
  sizeType density=std::atoi(argc[2]);
  std::string pathObj(argc[3]);

  //load hand
  std::experimental::filesystem::v1::path pathIO(path);
  pathIO.replace_extension("");
  pathIO.replace_filename(pathIO.filename().string()+"_"+std::to_string(density));
  pathIO.replace_extension(".dat");
  GraspPlanner<T> planner;
  if(exists(pathIO.string())) {
    planner.SerializableBase::read(pathIO.string());
  } else {
    planner.reset(path,1.0f/density);
    planner.fliterSample([&](sizeType lid,const Vec3T&,const Vec3T& n)->bool{
      if(lid==1)
        return n.dot(Vec3T(0,0,1))>0.9f;
      else if(lid==2)
        return n.dot(Vec3T(0,0,-1))>0.9f;
      else if(lid==5)
        return n.dot(Vec3T(0,0,1))>0.9f;
      else return n.dot(Vec3T(0,1,0))>0.9f;
    });
    planner.SerializableBase::write(pathIO.string());
  }

  //test objective
  GraspQualityMetric<T> obj;
  obj.SerializableBase::read(pathObj);
  Vec x0=Vec::Zero(planner.body().nrDOF());
  x0.template segment<3>(0)=Vec3T(0,0,-0.2f);

  pathIO=path;
  pathIO.replace_extension("");
  recreate(pathIO.filename().string());
  planner.writeVTK(x0,pathIO.filename().string(),1);
  planner.writeLocalVTK(pathIO.filename().string(),1);
  planner.writeLimitsVTK("limits");
  Vec x1=planner.optimize(x0,obj,3,20,Q_INF_MEAN,-1000,0,0,1e-2f,0);
  planner.writeVTK(x0,"beforeOptimize",1);
  planner.writeVTK(x1,"afterOptimize",1);
  obj.writeVTK("object",1);
  return 0;
}
