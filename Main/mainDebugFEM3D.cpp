#include <Deformable/FEMMesh.h>
#include <Deformable/FEMCell.h>
#include <Deformable/FEMSystem.h>
#include <Deformable/FEMOctreeMesh.h>

USE_PRJ_NAMESPACE

template <typename T>
void genMesh3DVoxel(const std::string& path,const std::string& pathEmbed,T res,bool oct=true,bool doubling=true)
{
  std::shared_ptr<FEMSystem<T>> systemFEM;
  if(exists(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solver.dat").string())) {
    systemFEM.reset(new FEMSystem<T>);
    systemFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solver.dat").string());
  } else {
    std::shared_ptr<FEMMesh<T>> meshFEM(new FEMMesh<T>(3));
    if(exists(std::experimental::filesystem::v1::path(path).filename().replace_extension(".mesh.dat").string())) {
      meshFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).filename().replace_extension(".mesh.dat").string());
    } else {
      tinyxml2::XMLDocument doc;
      doc.InsertEndChild(doc.NewElement("root"));
      put<std::string>(doc,"pathEmbed",pathEmbed);
      put<std::string>(doc,"path",path);
      put<std::string>(doc,"mesher","Voxel");
      put<T>(doc,"Voxel.res",res);
      put<bool>(doc,"forceRemesh",true);

      meshFEM.reset(new FEMMesh<T>(3));
      meshFEM->reset(doc);
      meshFEM->computeMatDist();
      meshFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).filename().replace_extension(".mesh.dat").string());
    }

    if(oct) {
      FEMOctreeMesh<T> oct(meshFEM->getB(0),0);
      meshFEM=oct.getMesh();
    }

    if(doubling) {
      FEMMesh<T> meshFEM2=*meshFEM;
      typename FEMMesh<T>::Mat4T R;
      R.setIdentity();
      R(2,3)=5.0f;
      meshFEM2.applyTrans(NULL,R,-1);
      *meshFEM+=meshFEM2;
    }

    systemFEM.reset(new FEMSystem<T>(meshFEM));
    systemFEM->addGravity(typename FEMMesh<T>::Vec3T(0,-9.81f,0));
    systemFEM->addNonHookeanElastic(1000.0f,10000.0f);
    systemFEM->fixVertices([&](const FEMVertex<T>& v) {
      return v.pos0().y()<-0.79f;
    });
    systemFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solver.dat").string());
  }

  ObjMesh obj;
  systemFEM->mesh().writeObj(NULL,obj);
  obj.writeVTK(std::experimental::filesystem::v1::path(path).filename().replace_extension(".obj.vtk").string(),true);

  systemFEM->debugSystem(10,0.0001f);
  systemFEM->mesh().writeVTK(NULL,std::experimental::filesystem::v1::path(path).filename().replace_extension(".vtk").string());
  systemFEM->mesh().writeEmbeddedMeshVTK(NULL,std::experimental::filesystem::v1::path(path).filename().replace_extension(".embed.vtk").string());
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  //genMesh3D<double>("data/FEM/deformedDino.abq","data/FEM/deformedDino.obj",0.05f);
  genMesh3DVoxel<double>("data/FEM/deformedDino.obj","data/FEM/deformedDino.obj",0.05f);
  return 0;
}
