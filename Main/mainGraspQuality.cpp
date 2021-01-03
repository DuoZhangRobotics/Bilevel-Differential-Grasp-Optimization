/**
 * * Get the target object prepared
*/
#include <Grasp/GraspQualityMetric.h>
#include <CommonFile/MakeMesh.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

typedef double T;
typedef GraspQualityMetric<T>::Vec Vec;
//TODO: fill up the scale argument
/**
    Function Command Line Input:
        @param argc[1]: Target object path //** [.obj file]
        @param argc[2]: density for sampling points on the surface of each link of the robot hand //** [how many points need to sample on the target object]
        @param argc[3]: target object file //*? Not Sure
    Function Output:
        ** Write processed target object to .dat file for optimization
*/
int main(int argn,char** argc)
{
  ASSERT_MSG(argn==4,"mainGraspQuality: [ObjMesh path] [sample density] [scale]")
  std::string path(argc[1]);
  sizeType density=std::atoi(argc[2]);
  T scale=std::atof(argc[3]);

  // Special case for cube
  ObjMesh m;
  if(beginsWith(path,"cube")) {
    MakeMesh::makeBox3D(m,Vec3::Constant(scale));
  } else {
    std::ifstream is(path.c_str());
    m.read(is,false,false);
    scale/=m.getBB().getExtent().maxCoeff();
    m.getScale()=scale;
    m.applyTrans();
  }
  
  // Read the obj file and write the result to .dat file
  std::experimental::filesystem::v1::path pathIO(path);
  pathIO.replace_extension("");
  pathIO.replace_filename(pathIO.filename().string()+"_"+std::to_string(density));
  pathIO.replace_extension(".dat");
  GraspQualityMetric<T> q;
  if(exists(pathIO.string())) {
    q.SerializableBase::read(pathIO.string());
  } else {
    q.reset(m, 1.0f/density);
    q.SerializableBase::write(pathIO.string());
  }
  q.debug(10);
  pathIO.replace_extension("");
  recreate(pathIO.filename().string());
  q.writeVTK(pathIO.filename().string(),1);
  return 0;
}
