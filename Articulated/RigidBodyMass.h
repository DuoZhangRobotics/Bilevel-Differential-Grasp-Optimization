#ifndef RIGID_BODY_MASS_H
#define RIGID_BODY_MASS_H

#include <CommonFile/ObjMesh.h>

PRJ_BEGIN

struct RigidBodyMass {
  RigidBodyMass(const ObjMesh& mesh);
  Vec3 getCtr() const;
  Mat6 getMass() const;
  Mat6 getMassCOM() const;
  scalar getM() const;
  Vec3 getMC() const;
  Mat3 getMCCT() const;
private:
  Mat6 _mat,_matCOM;
  Vec3 _ctr;
  Mat3 _MCCT;
};

PRJ_END

#endif
