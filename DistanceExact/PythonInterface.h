#ifndef PYTHON_INTERFACE_H
#define PYTHON_INTERFACE_H

#include "ObjMeshGeomCellExact.h"

extern "C"
{
  //distance
  void distance(const char** paths,int nrShape,const double* pss,double* dss,double* nss,double* hss,int nrPt);
}

#endif
