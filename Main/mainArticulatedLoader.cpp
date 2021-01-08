#include <Articulated/ArticulatedUtils.h>
#include <Articulated/ArticulatedLoader.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <CommonFile/geom/StaticGeom.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

int main(int argc,char** argv)
{
  ASSERT_MSG(argc>=2,"Usage: mainArticulatedLoader [path of urdf]")
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  parseProps(argc,argv,pt);
  scalar rad=get<scalar>(pt,"rad",0.1f);
  sizeType res=get<sizeType>(pt,"res",2);
  sizeType nrLink=get<sizeType>(pt,"nrLink",10);
  char visualMesh=get<sizeType>(pt,"visualMesh",0);
  char convex=get<sizeType>(pt,"convex",0);
  char simplify=get<sizeType>(pt,"simplify",1);
  char addBase=get<sizeType>(pt,"addBase",1);
  ArticulatedBody body;
  std::string pathURDF=argv[1];
  if(std::experimental::filesystem::v1::is_directory(pathURDF)) {
    for(std::experimental::filesystem::v1::directory_iterator it(pathURDF); it!=std::experimental::filesystem::v1::directory_iterator(); ++it)
      if(endsWith(it->path().string(),"urdf")||endsWith(it->path().string(),"URDF")) {
        std::string filePathURDF=std::experimental::filesystem::v1::absolute(it->path()).string();
        if(endsWith(filePathURDF,"urdf")||endsWith(filePathURDF,"URDF")) {
          ASSERT_MSGV(exists(filePathURDF),"Cannot find: %s!",filePathURDF.c_str())
          ArticulatedLoader::compareVisualMeshURDF(filePathURDF);
          body=ArticulatedLoader::readURDF(filePathURDF,convex,visualMesh);
        }
        break;
      }
  } else if(std::experimental::filesystem::v1::is_regular_file(pathURDF)) {
    ASSERT(endsWith(pathURDF,"urdf") || endsWith(pathURDF,"URDF"))
    ArticulatedLoader::compareVisualMeshURDF(pathURDF);
    body=ArticulatedLoader::readURDF(pathURDF,convex,visualMesh);
    pathURDF=std::experimental::filesystem::v1::path(pathURDF).parent_path().string();
  }
  if(body.nrJ()==0) {
    create(pathURDF);
    if(pathURDF=="bird")
      body=ArticulatedLoader::createBird(Joint::JOINT_TYPE(Joint::TRANS_3D|Joint::ROT_3D_EXP));
    else if(pathURDF=="bipedal")
      body=ArticulatedLoader::createBipedal(Joint::JOINT_TYPE(Joint::TRANS_3D|Joint::ROT_3D_EXP),true);
    else if(pathURDF=="chain")
      body=ArticulatedLoader::createChain(Joint::JOINT_TYPE(Joint::TRANS_3D|Joint::ROT_3D_EXP),rad,nrLink,res);
    else if(pathURDF=="spider")
      body=ArticulatedLoader::createSpider(Joint::JOINT_TYPE(Joint::TRANS_3D|Joint::ROT_3D_EXP));
    else if(pathURDF=="arm")
      body=ArticulatedLoader::createArm(1,0.1);
    else {
      ASSERT_MSG(false,"Unknown articulated body type!")
    }
  }
  {
    if(simplify) {
      ArticulatedUtils utils(body);
      utils.simplify(10);
    }
    if(addBase>0) {
      ArticulatedUtils utils(body);
      utils.addBase(addBase,Vec3d::UnitZ());
    }
    ASSERT_MSGV(body.nrJ()>0,"Cannot find ArticulatedBody using path: %s!",pathURDF.c_str())
    ArticulatedLoader::visualizeJointLimitVTK(body,pathURDF+"/visualize",Joint::MESH);
    PBDArticulatedGradientInfo<scalarD> info(body,Cold::Zero(body.nrDOF()));
    body.writeVTK(info._TM,pathURDF+"/visualize/mesh.vtk",Joint::MESH);
    body.getGeomEnv().writeVTK(pathURDF+"/visualize/geom.vtk");
    body.SerializableBase::write(pathURDF+"/ArticulatedBody.dat");
  }
  return 0;
}
