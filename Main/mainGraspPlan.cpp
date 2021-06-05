#include <Quasistatic/GraspPlanner.h>
#include <Quasistatic/CentroidClosednessEnergy.h>
#include <Quasistatic/ObjectClosednessEnergy.h>
#include <Quasistatic/LogBarrierObjEnergy.h>
#include <Quasistatic/ConvexLogBarrierSelfEnergy.h>
#include <Utils/Utils.h>
#include <string>
#include <fstream>
#include <chrono>

USE_PRJ_NAMESPACE

typedef double T;
typedef PointCloudObject<T>::Vec Vec;
typedef PointCloudObject<T>::Vec3T Vec3T;
Vec initializeParams(std::string path, Vec &x)
{
  std::cout << "Input path is: " <<  path << std::endl;
  std::string line;
  int i = 0;
  std::ifstream fin;
  fin.open((path));
  float a;
  std::cout << "Initial parameters are: ";
  while (fin >> a) {
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

  ASSERT_MSG(argn>=7,"mainGraspPlan: [urdf path] [sample density] [obj path] [obj name] [obj scale] [use_FGT] [max_iters] [saving dir] [initial parameters]")
  std::string path(argc[1]);
  sizeType density=std::atoi(argc[2]);
  std::string pathObj(argc[3]);
  std::cout << "Obj Path is: "<< pathObj << std::endl;
  std::string objName(argc[4]);
  std::string objScale(argc[5]);
  int useFGT = std::atoi(argc[6]);
  int max_iters = std::atoi(argc[7]);
  std::string initParamsPath;
  std::cout << "argn = " << argn << std::endl;
  std::string savingDir;
  if(argn>=9){
      std::string tmp(argc[8]);
      savingDir = tmp;
  }
  else{
      savingDir = "./";
  }
  if(argn==10) {
    std::string initParamPath(argc[9]);
    initParamsPath=initParamPath;
  }
  else initParamsPath="";
  std::cout << initParamsPath << std::endl;

  //load hand
  std::experimental::filesystem::v1::path pathIO(path);
  std::cout << "Path is: " << path << std::endl;
  pathIO.replace_extension("");
  pathIO.replace_filename(pathIO.filename().string()+"_"+std::to_string(density));
  pathIO.replace_extension(".dat");
  GraspPlanner<T> planner;
  ASSERT_MSG(exists(pathIO.string()),"Use mainGripper to create gripper first")
  planner.SerializableBase::read(pathIO.string());

  //test objective
  PointCloudObject<T> obj;
  obj.SerializableBase::read(pathObj);
  Vec x0=Vec::Zero(planner.body().nrDOF());
  std::string handName=" ";
  if(useFGT==1)
    handName="Q_INF_CONSTRAINT_FGT";
  else if(useFGT==2)
  {

      handName="Q_INF";
  }
  else if(useFGT==3)
  {

      handName="No_Metric_OC";
  }
  else{

      handName="Q_1";
  }
  Options ops;
  std::string type;
  GraspPlannerParameter param(ops);
  if(initParamsPath!="") {
    x0=initializeParams(initParamsPath, x0);
    if(pathIO.string().find("BarrettHand")!=std::string::npos) {
        handName+="_BarrettHand";
        type = "BarrettHand";
        param._coefM = -100;
        param._coefO=10;
        param._coefS=1;
    }
    else if(pathIO.string().find("ShadowHand")!=std::string::npos) {
        handName += "_Shadowhand";
        type = "Shadowhand";
        param._coefM=-1;
        param._coefO=100;
        param._coefS=1;
    }

  } else if(pathIO.string().find("BarrettHand")!=std::string::npos) {
      std::cout<< "find barretthand" << std::endl;
    x0.template segment<3>(0)=Vec3T(0,-0.2f,-0.2f);
    x0[5]=M_PI/2;
    x0[6]=0.5f;
    x0[9]=0.5f;
    handName+="_BarrettHand";
  } else if(pathIO.string().find("ShadowHand")!=std::string::npos) {
    x0.template segment<3>(0)=Vec3T(-0.1f,0.02f,-0.1f);
    x0[5]=-M_PI/2;

//     x0[3]=-M_PI/2;
//          x0[4]=-M_PI/2;
    handName += "_Shadowhand";
  }
   else{
      std::cout<< "wrong" << std::endl;
  }
//  pathIO=path;
//  pathIO.replace_extension("");
//  recreate(pathIO.filename().string());
//  planner.writeVTK(x0,pathIO.filename().string(),1);
  // planner.writeLocalVTK(pathIO.filename().string(),1);
  // planner.writeLimitsVTK("limits");
  std::string beforeOptimizeFileName= savingDir+"beforeOptimize_"+handName+ "_"+ objName+"_"+objScale;
  std::cout << "Initial parameters saved at: "<< beforeOptimizeFileName<< std::endl;
  planner.writeVTK(x0, beforeOptimizeFileName,1);
  std::ofstream initialParameters(beforeOptimizeFileName+"/initialParameters.txt");
  for(const auto &e:x0)
    initialParameters << e << " ";


  if(useFGT==1)
    param._metric=Q_INF_CONSTRAINT_FGT;
  else if(useFGT==2)
  {
      std::cout << "using Q_Inf" << std::endl;
      param._metric=Q_INF;

  }
  else if(useFGT==3)
  {
      std::cout << "using No_metric OC" << std::endl;
      param._metric=NO_METRIC;
      if(type=="Shadowhand")
        param._coefOC = 1;
      else if(type=="BarrettHand")
          param._coefOC = 100;
  }
  else{
      std::cout << "using Q_1" << std::endl;
      param._metric=Q_1;
  }
  if(max_iters==1)
  {
       param._normalExtrude=2;
       planner.evaluateQInf(x0, obj, param);
   }
  else
  {
      param._normalExtrude=10;
      param._maxIter=max_iters;
      x0=planner.optimize(false,x0,obj,param);
      auto start=std::chrono::high_resolution_clock::now();

      param._normalExtrude=2;
      param._maxIter=max_iters;
      x0=planner.optimize(false,x0,obj,param);
      auto stop=std::chrono::high_resolution_clock::now();
      auto duration=std::chrono::duration_cast<std::chrono::seconds>(stop - start);
      std::cout << "Optimization time using " << (useFGT?"FGT": "No FGT") << " is: " << duration.count() << std::endl;

      std::string afterOptimizeFileName=savingDir+"afterOptimize_"+handName+ "_" + objName+"_"+objScale;
      std::cout << "Output paramters saved at: " << afterOptimizeFileName << std::endl;
      planner.writeVTK(x0,afterOptimizeFileName, 1);
      obj.writeVTK("object",1,planner.rad()*param._normalExtrude);
      std::ofstream afterOptimizeFile(afterOptimizeFileName + "/parameters.txt");
      for(const auto &e:x0)
        afterOptimizeFile << e << " ";
  }

  return 0;

}
