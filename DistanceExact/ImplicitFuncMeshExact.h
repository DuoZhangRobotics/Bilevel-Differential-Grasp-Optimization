#ifndef IMPLICIT_FUNC_MESH_EXACT_H
#define IMPLICIT_FUNC_MESH_EXACT_H

#include "ObjMeshGeomCellExact.h"
#include "CommonFile/ImplicitFunc.h"

PRJ_BEGIN

class ImplicitFuncMeshExact : public ImplicitFunc<scalar>
{
public:
  ImplicitFuncMeshExact(const ObjMeshGeomCellExact& cell);
  virtual ~ImplicitFuncMeshExact();
  virtual scalar operator()(const PT& pos) const override;
  virtual BBox<ST> getBB() const override;
protected:
  const ObjMeshGeomCellExact& _cell;
};

PRJ_END

#endif
