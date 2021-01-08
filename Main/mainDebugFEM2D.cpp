#include <Deformable/FEMMesh.h>
#include <Deformable/FEMCell.h>
#include <Deformable/FEMReducedSystem.h>

USE_PRJ_NAMESPACE

template <typename T>
void genMesh2D(const std::string& path)
{
  std::shared_ptr<FEMSystem<T>> systemFEM;
  if(exists(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solver.dat").string())) {
    systemFEM.reset(new FEMSystem<T>);
    systemFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solver.dat").string());
  } else {
    std::shared_ptr<FEMMesh<T>> meshFEM(new FEMMesh<T>(2));
    if(exists(std::experimental::filesystem::v1::path(path).filename().replace_extension(".mesh.dat").string())) {
      meshFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).filename().replace_extension(".mesh.dat").string());
    } else {
      tinyxml2::XMLDocument doc;
      doc.InsertEndChild(doc.NewElement("root"));
      put<std::string>(doc,"path",path);
      put<bool>(doc,"SVG.normalize",true);
      put<scalar>(doc,"SVG.sz",0.001f);
      put<bool>(doc,"forceRemesh",true);

      meshFEM.reset(new FEMMesh<T>(2));
      meshFEM->reset(doc);
      meshFEM->computeMatDist();
      meshFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).filename().replace_extension(".mesh.dat").string());
    }

    FEMMesh<T> meshFEM2=*meshFEM;
    typename FEMMesh<T>::Mat4T R;
    R.setIdentity();
    R(0,3)=5.0f;
    meshFEM2.applyTrans(NULL,R,-1);
    *meshFEM+=meshFEM2;
    *meshFEM-=1;

    systemFEM.reset(new FEMSystem<T>(meshFEM));
    systemFEM->addGravity(typename FEMMesh<T>::Vec3T(0,-9.81f,0));
    systemFEM->addNonHookeanElastic(1000.0f,10000.0f);
    systemFEM->fixVertices([&](const FEMVertex<T>& v) {
      return v.pos0().y()<-0.45f;
    });
    systemFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).filename().replace_extension(".solver.dat").string());
  }

  systemFEM->debugSystem(10,0.0001f);
  systemFEM->mesh().writeVTK(NULL,std::experimental::filesystem::v1::path(path).filename().replace_extension(".vtk").string());
  systemFEM->mesh().writeObjExtrudeVTK(NULL,std::experimental::filesystem::v1::path(path).filename().replace_extension(".extrude.vtk").string(),0.01f);
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  genMesh2D<double>("data/FEM/Cross.svg");
  return 0;
}
