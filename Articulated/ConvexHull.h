#ifndef CONVEX_HULL_H
#define CONVEX_HULL_H

#include <CommonFile/ObjMesh.h>

PRJ_BEGIN

extern ObjMeshF makeConvex(const ObjMeshF& in);
extern ObjMeshD makeConvex(const ObjMeshD& in);

PRJ_END

#endif
