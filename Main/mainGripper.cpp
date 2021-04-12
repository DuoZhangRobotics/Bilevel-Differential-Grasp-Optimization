#include <Quasistatic/GraspPlanner.h>
#include <Quasistatic/PrimalDualQInfMetricEnergyFGT.h>
#include <Quasistatic/PrimalDualQInfMetricEnergy.h>
#include <Quasistatic/CentroidClosednessEnergy.h>
#include <Quasistatic/ObjectClosednessEnergy.h>
#include <Quasistatic/ConvexLogBarrierSelfEnergy.h>
#include <Quasistatic/LogBarrierObjEnergy.h>
#include <Quasistatic/MetricEnergy.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

typedef double T;
typedef PointCloudObject<T>::Vec Vec;
typedef PointCloudObject<T>::Vec3T Vec3T;
int main(int argn,char** argc)
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);

  ASSERT_MSG(argn==4,"mainGripper: [urdf path] [sample density] [obj path]")
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
    if(pathIO.string().find("BarrettHand")!=std::string::npos) {
      planner.fliterSample([&](sizeType lid,const Vec3T&,const Vec3T& n)->bool{
        if(lid==1)
          return n.dot(Vec3T(0,0,1))>0.9f;
        else if(lid==2)
          return n.dot(Vec3T(0,0,-1))>0.9f;
        else if(lid==5)
          return n.dot(Vec3T(0,0,1))>0.9f;
        else return n.dot(Vec3T(0,1,0))>0.9f;
      });
    } else if(pathIO.string().find("ShadowHand")!=std::string::npos) {
      planner.fliterSample([&](sizeType lid,const Vec3T&,const Vec3T& n)->bool{
        if(lid==20)
          return n.dot(Vec3T(0,0,1))>0.9f;
        else if(lid==22)
          return n.dot(Vec3T(-1,0,0))>0.9f;
        else return n.dot(Vec3T(0,-1,0))>0.9f;
      });
    }
    planner.SerializableBase::write(pathIO.string());
  }

  //test objective
  PointCloudObject<T> object;
  object.SerializableBase::read(pathObj);
  Vec x0=Vec::Zero(planner.body().nrDOF());
  if(pathIO.string().find("BarrettHand")!=std::string::npos)
    x0.template segment<3>(0)=Vec3T(0,0,-0.2f);
  else if(pathIO.string().find("ShadowHand")!=std::string::npos)
    x0.template segment<3>(0)=Vec3T(0,0.2f,0);
  
  pathIO=path;
  pathIO.replace_extension("");
  recreate(pathIO.filename().string());
  planner.writeVTK(x0,pathIO.filename().string(),1);
  planner.writeLocalVTK(pathIO.filename().string(),1);

  pathIO=pathObj;
  pathIO.replace_extension("");
  recreate(pathIO.filename().string());
  object.writeVTK(pathIO.filename().string(),1);

  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    MetricEnergy<T>(obj,info,planner,object,0,0.1f,10,Q_1,SQR_EXP_ACTIVATION,1).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    MetricEnergy<T>(obj,info,planner,object,0,0.1f,10,Q_INF,SQR_EXP_ACTIVATION,1).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    MetricEnergy<T>(obj,info,planner,object,0,0.1f,10,Q_INF_BARRIER,SQR_EXP_ACTIVATION,1).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    PrimalDualQInfMetricEnergy<T>(obj,info,planner,object,0.1f,10,SQR_EXP_ACTIVATION,1).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    PrimalDualQInfMetricEnergyFGT<T>(obj,info,planner,object,0.1f,10,1).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    CentroidClosednessEnergy<T>(obj,info,planner,object,10).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    ObjectClosednessEnergy<T>(obj,info,planner,object,10).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    LogBarrierObjEnergy<T>(obj,info,planner,object,0.1,10,false).debug(1,0,&x0,true);
  }
  {
    DSSQPObjectiveCompound<T> obj;
    PBDArticulatedGradientInfo<T> info;
    ConvexLogBarrierSelfEnergy<T> es(obj,info,planner,object,5,10,true);
    es.debug(1,0,&x0,true);
    es.updatePlanes();
  }
  

  return 0;
}
