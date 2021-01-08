#include <Articulated/ArticulatedLoader.h>
#include <Articulated/ArticulatedUtils.h>
#include <Deformable/FEMReducedSystem.h>
#include <Deformable/FEMOctreeMesh.h>
#include <Deformable/SoftFinger.h>
#include <Deformable/Hexahedron.h>
#include <Deformable/FEMEnergy.h>
#include <Deformable/FEMMesh.h>
#include <Deformable/FEMCell.h>
#include <Utils/RotationUtil.h>
#include <CommonFile/MakeMesh.h>
#include <Environment/Simulator.h>

USE_PRJ_NAMESPACE

//#define USE_POLY
template <typename T>
void genMesh3D(sizeType nrFinger,sizeType fingerTrans,sizeType baseThick,Vec3i objSz,Vec3i objPos,scalar res,
               sizeType nBasis=6,sizeType nCub=50,bool oct=true,bool forceRecreate=false)
{
  typename FEMMesh<T>::Vec3T g(0,-9.81f,0);

  //create finger
  std::string path="fingerSoft/finger.abq";
  if(forceRecreate || !exists("fingerSoft"))
    recreate("fingerSoft");
  std::shared_ptr<FEMReducedSystem<T>> systemFEM;
  if(exists(std::experimental::filesystem::v1::path(path).replace_extension(".solver.dat").string())) {
    systemFEM.reset(new FEMReducedSystem<T>);
    systemFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).replace_extension(".solver.dat").string());
  } else {
    std::shared_ptr<FEMMesh<T>> meshFEM(new FEMMesh<T>(3));
    if(exists(std::experimental::filesystem::v1::path(path).replace_extension(".mesh.dat").string())) {
      meshFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).replace_extension(".mesh.dat").string());
    } else {
      std::unordered_set<Vec3i,Hash> css=createSoftFinger(10,Vec4i(8,6,1,2),Vec2i(8,2),Vec2i(2,1),Vec2i(2,0));
      writeABQVoxel(css,"fingerSoft/finger.abq",res,true);
      writeVTKVoxel(css,"fingerSoft/finger.vtk",res);

      tinyxml2::XMLDocument doc;
      doc.InsertEndChild(doc.NewElement("root"));
      put<std::string>(doc,"path",path);

      meshFEM.reset(new FEMMesh<T>(3));
      meshFEM->reset(doc);
      meshFEM->computeMatDist();
      meshFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).replace_extension(".mesh.dat").string());
    }

    if(oct) {
      if(exists(std::experimental::filesystem::v1::path(path).replace_extension(".meshOCT.dat").string())) {
        meshFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).replace_extension(".meshOCT.dat").string());
      } else {
        FEMOctreeMesh<T> oct(meshFEM->getB(0),0,8,2);
        meshFEM=oct.getMesh();
        meshFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).replace_extension(".meshOCT.dat").string());
      }
    }
    meshFEM->writeVTK(NULL,"fingerSoft/fingerOctree.vtk");

    systemFEM.reset(new FEMReducedSystem<T>(meshFEM));
    systemFEM->addPressureEmbedded(0,0,1,typename FEMMesh<T>::Vec3T(1,0,0),0.9f,0.1f);
    systemFEM->addGravity(typename FEMMesh<T>::Vec3T(0,-9.81f,0));
    systemFEM->addStVKElastic(10000.0f,30000.0f);
    systemFEM->fixVertices([&](const FEMVertex<T>& v) {
      return v.pos0().x()<res*0.5f;
    });
    systemFEM->addGravity(g);
    systemFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).replace_extension(".solver.dat").string());
  }
  systemFEM->mesh().writeVTK(NULL,std::experimental::filesystem::v1::path(path).replace_extension(".vtk").string());
  systemFEM->mesh().writeEmbeddedMeshVTK(NULL,std::experimental::filesystem::v1::path(path).replace_extension(".embed.vtk").string());
  T fingerHeight=systemFEM->mesh().getBB(NULL).getExtent()[2];

  //build reduced model
  if(exists(std::experimental::filesystem::v1::path(path).replace_extension(".solverReduced.dat").string())) {
    systemFEM.reset(new FEMReducedSystem<T>);
    systemFEM->SerializableBase::read(std::experimental::filesystem::v1::path(path).replace_extension(".solverReduced.dat").string());
  } else {
    if(!systemFEM->buildBasis(nBasis,true,true))
      return;
    systemFEM->writeBasisVTK("fingerSoft/basis",0.1f);
    systemFEM->template buildPoly<FEMEmbeddedPressureEnergy<T>>(0);
#ifdef USE_POLY
    systemFEM->template buildPoly<FEMStVKElasticEnergy<T>>();
#else
    if(!systemFEM->optimizeCubature(nCub))
      return;
#endif
    systemFEM->SerializableBase::write(std::experimental::filesystem::v1::path(path).replace_extension(".solverReduced.dat").string());
  }

  //create base
  std::shared_ptr<ArticulatedBody> base;
  if(exists(std::experimental::filesystem::v1::path(path).replace_extension(".base.dat").string())) {
    base.reset(new ArticulatedBody);
    base->SerializableBase::read(std::experimental::filesystem::v1::path(path).replace_extension(".base.dat").string());
  } else {
    ObjMesh m;
    std::ofstream os(path);
    MakeMesh::makeCylinder3D(m,fingerTrans*res+fingerHeight,baseThick*res,32);
    m.getT()=expWGradV<scalar,Vec3>(Vec3::UnitZ()*M_PI/2);
    m.getPos()=-Vec3::UnitX()*baseThick*res;
    m.applyTrans(Vec3::Zero());
    ArticulatedBody b=ArticulatedLoader::createMesh(m);
    base.reset(new ArticulatedBody(b));
    base->SerializableBase::write(std::experimental::filesystem::v1::path(path).replace_extension(".base.dat").string());
  }
  {
    PBDArticulatedGradientInfo<T> info(*base,FEMMesh<T>::Vec::Zero(base->nrDOF()));
    base->writeVTK(info._TM.unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    }),"fingerSoft/base.vtk",Joint::MESH);
  }

  //create object
  std::shared_ptr<ArticulatedBody> object;
  if(exists(std::experimental::filesystem::v1::path(path).replace_extension(".object.dat").string())) {
    object.reset(new ArticulatedBody);
    object->SerializableBase::read(std::experimental::filesystem::v1::path(path).replace_extension(".object.dat").string());
  } else {
    ObjMesh m;
    std::ofstream os(path);
    MakeMesh::makeBox3D(m,objSz.template cast<scalar>()*res);
    ArticulatedBody b=ArticulatedLoader::createMesh(m);
    object.reset(new ArticulatedBody(b));
    ArticulatedUtils(*object).addBase(3,Vec3d::Zero());
    ArticulatedUtils(*object).simplify(10);
    object->SerializableBase::write(std::experimental::filesystem::v1::path(path).replace_extension(".object.dat").string());
  }
  {
    PBDArticulatedGradientInfo<T> info(*object,FEMMesh<T>::Vec::Zero(object->nrDOF()));
    object->writeVTK(info._TM.unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    }),"fingerSoft/object.vtk",Joint::MESH);
  }

  //create simulator
  std::shared_ptr<Simulator<T>> sim(new Simulator<T>);
  sim->addPBDObject(base,g,FEMMesh<T>::Vec::Zero(0));
  for(sizeType fid=0; fid<nrFinger; fid++) {
    typename Simulator<T>::Vec xInit;
    xInit.setZero(systemFEM->nrDOF());
    typename FEMMesh<T>::Mat3T R=expWGradV<T,typename FEMMesh<T>::Vec3T>(FEMMesh<T>::Vec3T::UnitX()*M_PI*2*fid/nrFinger);
    typename FEMMesh<T>::Vec3T t=R*FEMMesh<T>::Vec3T::UnitZ()*fingerTrans*res;
    xInit.template segment<3>(xInit.size()-6)=t;
    xInit.template segment<3>(xInit.size()-3)=invExpW<T>(R);
    sim->addFEMObject(systemFEM,xInit);
  }
  sim->attach(0,1);
  sim->attach(0,2);
  sim->attach(0,3);
  typename FEMMesh<T>::Vec6T pose;
  pose << objPos[0]*res,objPos[1]*res,objPos[2]*res,0,0,0;
  sim->addPBDObject(object,g,pose);
  //systemFEM->debugFST();
  //sim->debug(10,0.f);

  //simulate
  recreate("fingerSoft/sim");
  for(sizeType id=0; id<300; id++) {
    typename FEMMesh<T>::Mat3X4T t;
    t.setIdentity();
    t(2,3)=std::sin((T)id*0.3f)*0.1f;
    sim->setTrans(0,0,t);
    typename FEMMesh<T>::Vec ctrl=FEMMesh<T>::Vec::Zero(sim->nrControl());
    ctrl << 1,1,1,0,0,0,0,0,0;
    sim->solve(ctrl*4000);
    sim->advance();
    sim->writeVTK(false,"fingerSoft/sim/frm"+std::to_string(id)+".vtk");
  }
  sim->writeVTK(false,"fingerSoft/initSimulator.vtk");
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  genMesh3D<double>(3,28,10,Vec3i(10,10,10),Vec3i(40,0,0),0.025f,10,100);
  return 0;
}
