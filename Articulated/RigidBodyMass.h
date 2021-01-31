#ifndef RIGID_BODY_MASS_H
#define RIGID_BODY_MASS_H

#include <CommonFile/ObjMesh.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct RigidBodyMass
{
  DECL_MAP_TYPES_T
  RigidBodyMass(const ObjMesh& mesh);
  Vec3T getCtr() const;
  Mat6T getMass() const;
  Mat6T getMassCOM() const;
  T getM() const;
  Vec3T getMC() const;
  Mat3T getMCCT() const;
private:
  Mat6T _mat,_matCOM;
  Vec3T _ctr;
  Mat3T _MCCT;
};

PRJ_END

#endif
