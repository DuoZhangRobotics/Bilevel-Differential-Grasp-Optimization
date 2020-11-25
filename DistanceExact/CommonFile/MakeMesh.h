#ifndef MAKE_MESH_H
#define MAKE_MESH_H

#include "ObjMesh.h"

PRJ_BEGIN

class MakeMesh
{
public:
  //basic
  static void makeTet3D(ObjMesh& m);
  static void makeBox3D(ObjMesh& m,const Vec3& ext);
  static void makeBox3D(ObjMesh& m,const Vec3& ext,scalar thick);
  static void makeDiscreteBox3D(ObjMesh& m,const Vec3& ext);
  static void makeBox2D(ObjMesh& m,const Vec3& ext);
  static void makeBox2D(ObjMesh& m,const Vec3& ext,scalar thick);
  static void makeSphere3D(ObjMesh& m,const scalar rad,const sizeType slice);
  static void makeSphere2D(ObjMesh& m,const scalar rad,const sizeType slice);
  static void makeSphere3D(ObjMesh& m,const scalar rad,const sizeType slice,scalar thick);
  static void makeSphere2D(ObjMesh& m,const scalar rad,const sizeType slice,scalar thick);
  //capsule3D
  static void makeCapsule3D(ObjMesh& m,const scalar rad,const Vec3& offD,const sizeType slice);
  static void makeCapsule3D(ObjMesh& m,const scalar rad,const scalar x,const scalar y,const sizeType slice);
  static void makeCapsule3D(ObjMesh& m,const scalar rad,const scalar y,const sizeType slice);
  static void makeCapsule3D(ObjMesh& m,const scalar rad,const scalar x,const scalar y,const scalar z,const sizeType slice);
  //capsule2D
  static void makeCapsule2D(ObjMesh& m,const scalar rad,const Vec2& offD,const sizeType slice);
  static void makeCapsule2D(ObjMesh& m,const scalar rad,const scalar x,const scalar y,const sizeType slice);
  static void makeCapsule2D(ObjMesh& m,const scalar rad,const scalar y,const sizeType slice);
  //other
  static void makeCylinder3D(ObjMesh& m,const scalar rad,const scalar y,const sizeType slice,const sizeType sliceY=16,bool cap=true);
  static void makeTorus3D(ObjMesh& m,const scalar rad1,const scalar rad2,const sizeType slice1,const sizeType slice2);
  static void makeRing3D(ObjMesh& m,const scalar rad1,const scalar rad2,const scalar rad3,const sizeType slice);
  static void makeGridWithHole(ObjMesh& m,Vec4i slice);
  static void makeGrid(ObjMesh& m,const Vec2i& slice);
  //two sphere mesh
  static bool makeTwoSphereMesh3D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar lenY,const sizeType slice);
  static bool makeTwoSphereMesh2D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar lenY,const sizeType slice);
  static bool makeTwoSphereMesh(int dim,ObjMesh& m,const scalar rad1,const scalar rad2,scalar lenY,const sizeType slice);
  //three sphere mesh
  static bool makeThreeSphereMesh3D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar rad3,const Vec2& ctr1,const Vec2& ctr2,const Vec2& ctr3,const sizeType slice);
  static bool makeThreeSphereMesh2D(ObjMesh& m,const scalar rad1,const scalar rad2,scalar rad3,const Vec2& ctr1,const Vec2& ctr2,const Vec2& ctr3,const sizeType slice);
  static bool makeThreeSphereMesh(int dim,ObjMesh& m,scalar rad1,scalar rad2,scalar rad3,Vec2 ctr1,Vec2 ctr2,Vec2 ctr3,const sizeType slice);
protected:
  static Mat2 rot2D(scalar theta);
  static scalar getTheta(scalar rad1,scalar rad2,scalar lenY);
  static void hook(int dim,ObjMesh& m,std::vector<sizeType>& id0,std::vector<sizeType>& id1,bool close);
  static std::vector<sizeType> addSpherePatch(int dim,ObjMesh& m,const Vec3& n,scalar theta,sizeType slice,Vec3& x,Vec3& y,const Vec3& off);
  static std::pair<std::vector<sizeType>,std::vector<sizeType> > addSpherePatch(int dim,ObjMesh& mm,const Vec2& n1,const Vec2& n2,const Vec3& nn1,const Vec3& nn2,sizeType slice,const Vec2& off);
  static sizeType GI(sizeType i,sizeType j,sizeType si,sizeType sj,const std::map<sizeType,sizeType>& vMap);
  static sizeType GI(sizeType i,sizeType j,sizeType si,sizeType sj);
};

PRJ_END

#endif
