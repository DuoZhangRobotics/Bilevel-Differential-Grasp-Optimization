#include "PythonInterface.h"
#include "CommonFile/Timing.h"
#include "Utils/Utils.h"
#include <fstream>

extern "C"
{
  //distance
  std::vector<std::string> objPaths;
  std::vector<COMMON::ObjMeshGeomCellExact> objs;
  void distance(const char** paths,int nrShape,const double* pss,double* dss,double* nss,double* hss,int nrPt)
  {
    using namespace COMMON;
    disableTiming();

    //get shape
    TBEG("ReadObjMeshGeomCell");
    objPaths.resize(nrShape,"");
    objs.resize(nrShape);
    //OMP_PARALLEL_FOR_
    for(sizeType i=0; i<nrShape; i++)
      if(objPaths[i]!=paths[i]) {
        ASSERT_MSGV(COMMON::exists(paths[i]),"Cannot find: %s!",paths[i])
        ASSERT_MSGV(COMMON::endsWith(paths[i],"obj") || COMMON::endsWith(paths[i],"dat"),"Path %s is not ending with obj or dat!",paths[i])
        if(endsWith(paths[i],"obj")) {
          ObjMesh mesh;
          std::ifstream is(paths[i]);
          mesh.read(is,false,false);
          mesh.makeUniform();
          if(mesh.getVolume()<0)
            mesh.insideOut();
          COMMON::ObjMeshGeomCell cell(Mat4::Identity(),mesh,0,true);
          objs[i]=COMMON::ObjMeshGeomCellExact(cell);
        } else {
          objs[i].SerializableBase::read(paths[i]);
        }
        //std::cout << "loaded mesh with " << objs[i].nrV() << " vertices " << objs[i].nrI() << " faces!" << std::endl;
        objPaths[i]=paths[i];
      }
    TEND();

    //get distance
    TBEG("GetDistance");
    Eigen::Map<const Eigen::Matrix<double,-1,-1,Eigen::RowMajor>> pssMap(pss,3*nrShape,nrPt);
    Eigen::Map<Eigen::Matrix<double,-1,1>> dssMap(dss,nrPt*nrShape);
    Eigen::Map<Eigen::Matrix<double,3,-1>> nssMap(nss,3,nrPt*nrShape);
    Eigen::Map<Eigen::Matrix<double,3,-1>> hssMap(hss,3,nrPt*3*nrShape);
    Vec3 n,normal;
    Mat3 hessian;
    Vec2i feat;
//#define DEBUG_MODE
#ifndef DEBUG_MODE
    OMP_PARALLEL_FOR_I(OMP_PRI(n,normal,hessian,feat))
#endif
    for(sizeType i=0; i<nrShape*nrPt; i++) {
      sizeType shapeId=i/nrPt;
      sizeType ptId=i%nrPt;
      Vec3 pt=pssMap.block<3,1>(shapeId*3,ptId).cast<scalar>();
      dssMap[i]=objs[shapeId].closest(pt,n,normal,hessian,feat);
#ifdef DEBUG_MODE
      std::cout << "Shape" << shapeId << " Point" << ptId << std::endl;
      std::cout << "Point:" << pt.transpose() << std::endl;
      std::cout << "Normal:" << normal.transpose() << std::endl;
      std::cout << "Hessian:\n" << hessian << std::endl;
#endif
      nssMap.col(i)=normal.cast<double>();
      hssMap.block<3,3>(0,i*3)=hessian.cast<double>();
    }
    TEND();
  }
}
