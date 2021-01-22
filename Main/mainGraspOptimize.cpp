#include <Grasp/GraspPlanner.h>
#include <Grasp/CentroidClosednessEnergy.h>
#include <Grasp/ObjectClosednessEnergy.h>
#include <Grasp/LogBarrierObjEnergy.h>
#include <Grasp/LogBarrierSelfEnergy.h>
#include <Utils/Utils.h>
#include <string>
#include <fstream>

USE_PRJ_NAMESPACE

typedef double T;
typedef GraspQualityMetric<T>::Vec Vec;
typedef GraspQualityMetric<T>::Vec3T Vec3T;

Vec initializeParams(std::string path, Vec &x)
{
    std::string line;
    int i = 0;
    std::ifstream fin;
    fin.open((path));
    float a;
    std::cout << "Initial parameters are: ";
    while (fin >> a)
    {
      x[i] = a;
      std::cout << x[i] << " ";
      i++;
    }
    std::cout << std::endl;
    return x;
}

int main(int argn,char** argc)
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);

  ASSERT_MSG(argn>=6,"mainGraspPlan: [urdf path] [sample density] [obj path] [obj name] [obj scale] [initial parameters]")
  std::string path(argc[1]);
  sizeType density=std::atoi(argc[2]);
  std::string pathObj(argc[3]);
  std::string objName(argc[4]);
  std::string objScale(argc[5]);
  std::string initParamsPath;
  std::cout << "argn = " << argn << std::endl;
  if (argn == 7)
  {
      std::string initParamPath(argc[6]);
      initParamsPath = initParamPath;
  }
  else initParamsPath = "";
  std::cout << initParamsPath << std::endl;
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
  if (initParamsPath != "")
  {
    x0 = initializeParams(initParamsPath, x0);
  }
  else
  {
    x0.template segment<3>(0)=Vec3T(0.0f, -0.23f, -0.19f);
    x0[5]=M_PI/2;
    x0[6]=0.5f;
    x0[9]=0.5f;

  }
  x0[1] -= 0.1;
  pathIO=path;
  pathIO.replace_extension("");
  recreate(pathIO.filename().string());
  planner.writeVTK(x0,pathIO.filename().string(),1);
  planner.writeLocalVTK(pathIO.filename().string(),1);
  planner.writeLimitsVTK("limits");
  std::string beforeOptimizeFileName = "beforeOptimize_" + objName + "_" + objScale;
  planner.writeVTK(x0, beforeOptimizeFileName, 1);
  std::ofstream initialParameters(beforeOptimizeFileName + "/initialParameters.txt");
  for (const auto &e : x0) initialParameters << e << " ";
  Options ops;
  GraspPlannerParameter param(ops);
  param._normalExtrude=10;
  param._maxIter=15000;
  x0=planner.optimize(false,x0,obj,param);
  param._normalExtrude=2;
  param._maxIter=15000;
  x0=planner.optimize(false,x0,obj,param);
  std::string afterOptimizeFileName = "afterOptimize_" + objName + "_" + objScale;
  planner.writeVTK(x0, afterOptimizeFileName, 1);
  obj.writeVTK("object",1,planner.rad()*param._normalExtrude);
  std::ofstream afterOptimizeFile(afterOptimizeFileName + "/parameters.txt");
  for (const auto &e : x0) afterOptimizeFile << e << " ";

  return 0;
}
