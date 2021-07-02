#include <Quasistatic/PointCloudObject.h>
#include <CommonFile/MakeMesh.h>
#include <Utils/Utils.h>
#include <fstream>
USE_PRJ_NAMESPACE

typedef double T;
typedef PointCloudObject<T>::Vec Vec;
DECL_MAP_TYPES_T
Vec maxRange(PointCloudObject<T> q){
    Mat3XT vertices = q.pss();
    Vec x_ax = vertices.row(0);
    Vec y_ax = vertices.row(1);
    Vec z_ax = vertices.row(2);
    Vec x = Vec::Zero(3);
    x[0] = x_ax.maxCoeff();
    x[1] = y_ax.maxCoeff();
    x[2] = z_ax.maxCoeff();
//    std::cout << x<< std::endl;
    return x;
}

int main(int argn,char** argc)
{
  ASSERT_MSG(argn>=4,"mainPointCloudObject: [ObjMesh path] [radius of disk] [scale] [scaleY]")
  std::string path(argc[1]);
  sizeType density=std::atoi(argc[2]);
  T scale=std::atof(argc[3]);
  T scaleY=std::atof(argn>4?argc[4]:argc[3]);
  bool graspable=true;
  // bool graspable=false;
  ObjMesh m;
  if(beginsWith(path,"cube"))
    MakeMesh::makeBox3D(m,Vec3::Constant(scale));
  else if(beginsWith(path,"sphere"))
    MakeMesh::makeSphere3D(m,scale,16);
  else if(beginsWith(path,"cylinder"))
    MakeMesh::makeCylinder3D(m,scale,scaleY,16);
  else if(beginsWith(path,"plane")) {
    MakeMesh::makeGrid(m,Vec2i(10,10));
    m.getPos()=-Vec3(0.5f,0.5f,0);
    m.applyTrans();
    m.getScale()=scale;
    m.applyTrans();
    graspable=false;
  } else {
    std::ifstream is(path.c_str());
    m.read(is,false,false);
    scale/=m.getBB().getExtent().maxCoeff();

    std::cout << "scale1 = " << m.getScale()<< " " << m.getBB().getExtent().maxCoeff() << std::endl;

    m.getScale()=scale;
    std::cout << "x, y, z = "<< m.getVolumeCentroid()[0] << " " << m.getVolumeCentroid()[1]<< " "<< m.getVolumeCentroid()[2]<< std::endl;
    std::cout << "x, y, z = "<< m.getPos()[0] << " " << m.getPos()[1]<< " "<< m.getPos()[2]<< std::endl;
//    m.getPos()-= m.getVolumeCentroid();
    m.applyTrans();
    std::cout << "scale2 = " << m.getScale() << std::endl;
     std::cout << "x, y, z = "<< m.getVolumeCentroid()[0] << " " << m.getVolumeCentroid()[1]<< " "<< m.getVolumeCentroid()[2]<< std::endl;
     std::cout << "x, y, z = "<< m.getPos()[0] << " " << m.getPos()[1]<< " "<< m.getPos()[2]<< std::endl;
     std::cout << "scale = " << scale << std::endl;

  }

  std::experimental::filesystem::v1::path pathIO(path);
  pathIO.replace_extension("");
  pathIO.replace_filename(pathIO.filename().string()+"_"+std::to_string(density) + "_" + std::to_string(scaleY) + ".tmp");
  pathIO.replace_extension(".dat");
  PointCloudObject<T> q;
  if(exists(pathIO.string())) {
    q.SerializableBase::read(pathIO.string());
  } else {
    if(graspable) {
      q.resetGraspable(m,1.0f/density);
    } else {
      q.reset(m,1.0f/density);
    }
    q.SerializableBase::write(pathIO.string());
  }
  Vec max = maxRange(q);
  std::cout << max[0] <<" "<< max[1]<<" "<<max[2]<<std::endl;
  std::ofstream file("maxRange_Scale.txt");
  for(sizeType i=0; i<max.size(); i++){
    file << max[i] << " ";
}
  file << "\n";
  file<<scale;
  q.debug(10);
  pathIO.replace_extension("");
  recreate(pathIO.filename().string());
  q.writeVTK(pathIO.filename().string(),0);
  m.writeVTK(pathIO.filename().string()+".vtk",1);
  return 0;
}
