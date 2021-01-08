#include <Deformable/FEMMesh.h>
#include <Deformable/FEMCell.h>
#include <Deformable/FEMEnergy.h>
#include <Deformable/FEMReducedSystem.h>
#include <Deformable/FEMOctreeMesh.h>

USE_PRJ_NAMESPACE

template <typename T>
void genMesh3DVoxel(const std::string& path,const std::string& pathEmbed,T res,bool free,bool oct=true)
{
  std::string freeStr=free?"Free":"";
  std::shared_ptr<FEMReducedSystem<T>> systemFEM;
  if(exists(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solverReduced"+freeStr+".dat").string())) {
    systemFEM.reset(new FEMReducedSystem<T>);
    systemFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solverReduced"+freeStr+".dat").string());
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

    systemFEM.reset(new FEMReducedSystem<T>(meshFEM));
    systemFEM->addGravity(typename FEMMesh<T>::Vec3T(0,-9.81f,0));
    systemFEM->addStVKElastic(1000.0f,10000.0f);
    if(!free)
      systemFEM->fixVertices([&](const FEMVertex<T>& v) {
      return v.pos0().y()<-0.79f;
    });
    systemFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solverReduced"+freeStr+".dat").string());
  }
  if(exists(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solverReducedBasis"+freeStr+".dat").string()))
    systemFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solverReducedBasis"+freeStr+".dat").string());
  else {
    systemFEM->buildBasis(6,true,true,true,1e-3f);
    systemFEM->optimizeCubature();
    systemFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solverReducedBasis"+freeStr+".dat").string());
  }
  systemFEM->debugFST();
  systemFEM->debugSystem(10,0.0001f);
  systemFEM->writeBasisVTK("basis",0.3f);
  systemFEM->writeCubatureVTK("basis/cub.vtk");
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  genMesh3DVoxel<double>("data/FEM/deformedDino.obj","data/FEM/deformedDino.obj",0.05f,false);
  genMesh3DVoxel<double>("data/FEM/deformedDino.obj","data/FEM/deformedDino.obj",0.05f,true);
  return 0;
}
