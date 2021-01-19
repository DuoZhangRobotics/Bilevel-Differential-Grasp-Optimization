#include <Grasp/GraspQualityMetric.h>
#include <CommonFile/MakeMesh.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

typedef double T;
typedef GraspQualityMetric<T>::Vec Vec;
int main(int argn,char** argc)
{
  ASSERT_MSG(argn>=4,"mainGraspQuality: [ObjMesh path] [radius of disk] [scale] [scaleY]")
  std::string path(argc[1]);
  sizeType density=std::atoi(argc[2]);
  T scale=std::atof(argc[3]);
  T scaleY=std::atof(argn>4?argc[4]:argc[3]);

  ObjMesh m;
  if(beginsWith(path,"cube"))
    MakeMesh::makeBox3D(m,Vec3::Constant(scale));
  else if(beginsWith(path,"sphere"))
    MakeMesh::makeSphere3D(m,scale,16);
  else if(beginsWith(path,"cylinder"))
    MakeMesh::makeCylinder3D(m,scale,scaleY,16);
  else {
    std::ifstream is(path.c_str());
    m.read(is,false,false);
    scale/=m.getBB().getExtent().maxCoeff();
    m.getScale()=scale;
    m.applyTrans();
  }

  std::experimental::filesystem::v1::path pathIO(path);
  pathIO.replace_extension("");
  if (argn < 4){
      pathIO.replace_filename(pathIO.filename().string()+"_"+std::to_string(density) + "_" + std::to_string(scale) + ".tmp");
  }
  else{
      pathIO.replace_filename(pathIO.filename().string()+"_"+std::to_string(density) + "_" + std::to_string(scale) + "_" + std::to_string(scaleY) + ".tmp");
  }

  pathIO.replace_extension(".dat");
  GraspQualityMetric<T> q;
  if(exists(pathIO.string())) {
    q.SerializableBase::read(pathIO.string());
  } else {
    q.reset(m,1.0f/density);
    q.SerializableBase::write(pathIO.string());
  }
  q.debug(10);
  pathIO.replace_extension("");
  recreate(pathIO.filename().string());
  q.writeVTK(pathIO.filename().string(),1);
  return 0;
}
