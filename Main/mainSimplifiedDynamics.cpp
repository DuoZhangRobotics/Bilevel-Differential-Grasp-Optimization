#include <Articulated/ArticulatedBody.h>
#include <Articulated/SimplifiedDynamics.h>
#include <Utils/ArticulatedBodyPragma.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

int main(int argc,char** argv)
{
#ifdef OPTIMIZER_SUPPORT
  ASSERT_MSG(argc>=2,"Usage: mainSimplifiedDynamics [path of urdf]")
  std::string pathURDF=argv[1];
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  parseProps(argc,argv,pt);
  scalar res=get<scalar>(pt,"res",1.0f);
  sizeType symDir=get<sizeType>(pt,"symDir",-1);
  Vec3d zRange=parsePtreeDef<Vec3d>(pt,"zRange",Vec3d::Zero());
  scalarD zMargin=get<scalar>(pt,"zMargin",0.0f);

  std::shared_ptr<ArticulatedBody> body(new ArticulatedBody);
  body->SerializableBase::read(pathURDF+"/ArticulatedBody.dat");
  if(exists(pathURDF+"/SimplifiedDynamics.dat")) {
    SimplifiedDynamics simp;
    simp.SerializableBase::read(pathURDF+"/SimplifiedDynamics.dat");
    simp.cullSymmetry(symDir);
    simp.writeVTK(pathURDF);
    simp.debugDOF(pathURDF);
  } else {
    SimplifiedDynamics simp(body,res,symDir,zRange,zMargin);
    simp.SerializableBase::write(pathURDF+"/SimplifiedDynamics.dat");
    simp.writeVTK(pathURDF);
    simp.debugDOF(pathURDF);
  }
#endif
  return 0;
}
