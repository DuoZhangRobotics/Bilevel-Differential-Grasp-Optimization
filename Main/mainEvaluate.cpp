#include <Quasistatic/GraspPlanner.h>
#include <Quasistatic/CentroidClosednessEnergy.h>
#include <Quasistatic/ObjectClosednessEnergy.h>
#include <Quasistatic/LogBarrierObjEnergy.h>
#include <Quasistatic/ConvexLogBarrierSelfEnergy.h>
#include <Utils/Utils.h>
#include <string>
#include <fstream>

USE_PRJ_NAMESPACE

typedef double T;
typedef PointCloudObject<T>::Vec Vec;
typedef PointCloudObject<T>::Vec3T Vec3T;
Vec initializeParams(std::string path, Vec &x)
{
  // std::cout << "Input path is: " <<  path << std::endl;
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
  // std::cout << "Obj Path is: "<< pathObj << std::endl;
  std::string objName(argc[4]);
  std::string objScale(argc[5]);
  int useFGT = std::atoi(argc[6]);
//  int max_iters = std::atoi(argc[7]);
  std::string initParamsPath;
  // std::cout << "argn=" << argn << std::endl;
  std::string savingDir;
  if(argn>=9) {
    std::string tmp(argc[8]);
    savingDir = tmp;
  } else {
    savingDir="./";
  }
  if(argn>=10) {
    std::string initParamPath(argc[9]);
    initParamsPath=initParamPath;
  } else initParamsPath="";
  // std::cout << initParamsPath << std::endl;

  //load hand
  std::experimental::filesystem::v1::path pathIO(path);
  // std::cout << "Path is: " << path << std::endl;
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
  std::cout << "dofnum = " << planner.body().nrDOF() << std::endl;
  std::string handName=" ";
  if(useFGT==1)
    handName="Q_INF_CONSTRAINT_FGT";
  else if(useFGT==2)
    handName="Q_INF";
  else if(useFGT==3)
    handName="No_Metric_OC";
  else
    handName="Q_1";
  Options ops;
  std::string type;
  GraspPlannerParameter param(ops);
  if(argn>=11) {
    param._FGTThres=std::atof(argc[10]);
    std::cout << "setting FGTThres=" << param._FGTThres << std::endl;
  }
  if(initParamsPath!="") {
    x0=initializeParams(initParamsPath, x0);
    if(pathIO.string().find("BarrettHand")!=std::string::npos) {
      handName+="_BarrettHand";
      type = "BarrettHand";
      param._coefM = -100;
      param._coefO=10;
      param._coefS=1;
    } else if(pathIO.string().find("ShadowHand")!=std::string::npos) {
      handName += "_Shadowhand";
      type = "Shadowhand";
      param._coefM=-1;
      param._coefO=100;
      param._coefS=1;
    }
  } else if(pathIO.string().find("BarrettHand")!=std::string::npos) {
    std::cout << "find barretthand" << std::endl;
    x0.template segment<3>(0)=Vec3T(0,-0.0f,-0.0f);
//    x0[5]=M_PI/2;
//    x0[4]=M_PI/2;
//    x0[6]=0.5f;
//    x0[9]=0.5f;
    handName+="_BarrettHand";
  } else if(pathIO.string().find("ShadowHand")!=std::string::npos) {
    x0.template segment<3>(0)=Vec3T(-0.1f,0.02f,-0.1f);
    x0[5]=-M_PI/2;

    // x0[3]=-M_PI/2;
    // x0[4]=-M_PI/2;
    handName += "_Shadowhand";
  } else {
    std::cout<< "wrong" << std::endl;
  }
  std::cout << "Hand volume = " << planner.body().totalMass() << std::endl;
  // pathIO=path;
  // pathIO.replace_extension("");
  // recreate(pathIO.filename().string());
  // planner.writeVTK(x0,pathIO.filename().string(),1);
  // planner.writeLocalVTK(pathIO.filename().string(),1);
  // planner.writeLimitsVTK("limits");
  
  // std::cout << "Initial parameters saved at: "<< beforeOptimizeFileName<< std::endl;
  for (double i = 0; i < 20; i++)
  {
    std::string beforeOptimizeFileName=savingDir+"TestRange"+std::to_string(i);
    x0[12] = i/10;
    planner.writeVTK(x0, beforeOptimizeFileName,1);
  }
  
  // planner.writeVTK(x0, beforeOptimizeFileName,1);
  // std::ofstream initialParameters(beforeOptimizeFileName+"/initialParameters.txt");
  // for(sizeType i=0; i<x0.size(); i++)
  //   initialParameters << x0[i] << " ";
  return 0;
}
