#ifndef OBJ_MESH_H
#define OBJ_MESH_H

#include "MathBasic.h"
#include "IOFwd.h"
#include "CollisionDetection.h"
#include <vector>
#include <map>

PRJ_BEGIN

template <typename T>
struct ObjMeshTraits;

#define SCALAR_NAME scalarF
#define OBJMESH ObjMeshF
#include "ObjMeshHeader.h"
#undef OBJMESH
#undef SCALAR_NAME

#define SCALAR_NAME scalarD
#define OBJMESH ObjMeshD
#include "ObjMeshHeader.h"
#undef OBJMESH
#undef SCALAR_NAME

#ifdef DOUBLE_PRECISION
typedef ObjMeshTraits<scalarD>::Type ObjMesh;
#else
typedef ObjMeshTraits<scalarF>::Type ObjMesh;
#endif

PRJ_END

#endif
