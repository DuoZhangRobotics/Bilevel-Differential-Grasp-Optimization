#include "ImplicitFuncMeshExact.h"
#include "MPQZIO.h"

USE_PRJ_NAMESPACE

ImplicitFuncMeshExact::ImplicitFuncMeshExact(const ObjMeshGeomCellExact& cell):_cell(cell) {}
ImplicitFuncMeshExact::~ImplicitFuncMeshExact() {}
scalar ImplicitFuncMeshExact::operator()(const PT& pos) const
{
  Vec2i feat;
  Mat3 hessian;
  Vec3 n,normal;
  return _cell.closest(pos,n,normal,hessian,feat);
}
BBox<scalar> ImplicitFuncMeshExact::getBB() const
{
  BBox<scalar> ret;
  const BBoxExact& bb=_cell.getBB();
  ret._minC=castRational<Vec3,ObjMeshGeomCellExact::PT>(bb._minC);
  ret._maxC=castRational<Vec3,ObjMeshGeomCellExact::PT>(bb._maxC);
  return ret;
}
