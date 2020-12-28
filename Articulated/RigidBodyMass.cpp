#include "RigidBodyMass.h"
#include <Utils/RotationUtil.h>
#include <CommonFile/MakeMesh.h>
#include <CommonFile/geom/StaticGeom.h>

USE_PRJ_NAMESPACE

#define X 0
#define Y 1
#define Z 2
#define SQR(A) ((A)*(A))
#define CUBE(A) ((A)*(A)*(A))

struct RigidBodyMass2D {
  typedef std::vector<Vec3,Eigen::aligned_allocator<Vec3> > POS;
  typedef std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > IDS;
public:
  //just use Christopher Batty's code
  RigidBodyMass2D(const POS& V,const IDS& I):_V(V),_I(I) {
    _m = 0;
    _Cx = 0;
    _Cy = 0;
    _xx = 0;
    _yy = 0;
    _xy = 0;

    for(sizeType i=0; i<(sizeType)I.size(); i++) {
      Vec3 v1=V[I[i][0]];
      Vec3 v2=V[I[i][1]];

      // Signed area of this triangle (where 3rd vertex is the origin).
      scalar v = v1[0]*v2[1] - v2[0]*v1[1];
      // Contribution to the mass
      _m += v;
      // Contribution to the centroid
      scalar x4 = v1[0] + v2[0];
      _Cx += (v * x4);
      scalar y4 = v1[1] + v2[1];
      _Cy += (v * y4);
      // Contribution to moment of inertia monomials
      _xx += v * (v1[0]*v1[0] + v2[0]*v2[0] + x4*x4);
      _yy += v * (v1[1]*v1[1] + v2[1]*v2[1] + y4*y4);
      _xy += v * (v1[0]*v1[1] + v2[0]*v2[1] + x4*y4);
    }

    _mat.setZero();
    // Centroid.
    // The case _m = 0 needs to be addressed here.
    scalar r = 1.0f / (3.0f * _m);
    scalar Cx = _Cx * r;
    scalar Cy = _Cy * r;
    _ctr=Vec3(Cx,Cy,0);
    // Mass
    scalar m=_m / 2.0f;
    // Moment of inertia about the centroid.
    r = 1.0f / 24.0f;
    Mat3 J=Mat3::Zero();
    _xx = _xx * r;
    _yy = _yy * r;
    _xy = _xy * r;
    J(2,2)=_xx+_yy;
    // translate inertia tensor to center of mass
    translateCOMInv(J,Vec3(Cx,Cy,0.0f),m);
    _mat(0,0)=_mat(1,1)=m;
    _mat.block<3,3>(3,3) = J;
    // translate again
    _matCOM=_mat;
    translateCOM(_mat,Vec3(Cx,Cy,0.0f),m);
  }
  //getter
  virtual Vec3 getCtr() const {
    return _ctr;
  }
  virtual Mat6 getMass() const {
    return _mat;
  }
  virtual Mat6 getMassCOM() const {
    return _matCOM;
  }
  virtual Mat3 getMCCT() const {
    Mat3 ret=Mat3::Zero();
    ret(0,0)=_xx;
    ret(1,1)=_yy;
    ret(0,1)=ret(1,0)=_xy;
    return ret;
  }
protected:
  static void translateCOM(Mat6& J,const Vec3& r,const scalar& mass) {
    Mat3 cr=cross<scalar>(r);
    J.block<3,3>(3,3)+=cr*cr.transpose()*mass;
    J.block<3,3>(0,3)=cr.transpose()*mass;
    J.block<3,3>(3,0)=cr*mass;
  }
  static void translateCOMInv(Mat3& J,const Vec3& r,const scalar& mass) {
    J(0,0) -= mass * (r[Y]*r[Y] + r[Z]*r[Z]);
    J(1,1) -= mass * (r[Z]*r[Z] + r[X]*r[X]);
    J(2,2) -= mass * (r[X]*r[X] + r[Y]*r[Y]);
    J(1,0) = J(0,1) += mass * r[X] * r[Y];
    J(2,1) = J(1,2) += mass * r[Y] * r[Z];
    J(2,0) = J(0,2) += mass * r[Z] * r[X];
  }
  Mat6 _mat,_matCOM;
  Vec3 _ctr;
  POS _V;
  IDS _I;
private:
  scalar _m;
  scalar _Cx, _Cy;
  scalar _xx, _yy, _xy;
};
struct RigidBodyMass3D : public RigidBodyMass2D {
public:
  //just use Brian Mirtich's code
  RigidBodyMass3D(const POS& V,const IDS& I):RigidBodyMass2D(V,I) {
    _TN.resize(_I.size());
    for(sizeType i=0; i<(sizeType)_I.size(); i++) {
      Vec3i T=_I[i];
      _TN[i]=(_V[T[1]]-_V[T[0]]).cross(_V[T[2]]-_V[T[0]]);
      if(_TN[i].norm()>ScalarUtil<scalar>::scalar_eps())
        _TN[i].normalize();
      else {
        //remove zero-volume triangles
        _TN[i]=_TN.back();
        _TN.pop_back();
        _I[i]=_I.back();
        _I.pop_back();
        i--;
      }
    }
    compVolumeIntegrals();
    _mat.setIdentity();
    // compute center of mass
    Vec3 r(_T1[X] / _T0, _T1[Y] / _T0, _T1[Z] / _T0);
    _ctr=r;
    // compute inertia tensor
    Mat3 J=Mat3::Zero();
    J(0,0) = (_T2[Y] + _T2[Z]);
    J(1,1) = (_T2[Z] + _T2[X]);
    J(2,2) = (_T2[X] + _T2[Y]);
    J(1,0) = J(0,1) = -_TP[X];
    J(2,1) = J(1,2) = -_TP[Y];
    J(2,0) = J(0,2) = -_TP[Z];
    // translate inertia tensor to center of mass
    translateCOMInv(J,r,_T0);
    //set to matrix
    _mat.block<3,3>(0,0)=Mat3::Identity()*_T0;
    _mat.block<3,3>(3,3)=J;
    // translate again
    _matCOM=_mat;
    translateCOM(_mat,r,_T0);
  }
  void compProjectionIntegrals(const Vec3i& f) {
    scalar a0, a1, da;
    scalar b0, b1, db;
    scalar a0_2, a0_3, a0_4, b0_2, b0_3, b0_4;
    scalar a1_2, a1_3, b1_2, b1_3;
    scalar C1, Ca, Caa, Caaa, Cb, Cbb, Cbbb;
    scalar Cab, Kab, Caab, Kaab, Cabb, Kabb;
    sizeType i;

    _P1 = _Pa = _Pb = _Paa = _Pab = _Pbb = _Paaa = _Paab = _Pabb = _Pbbb = 0.0f;

    for (i = 0; i < 3; i++) {
      a0 = _V[f[i]][_A];
      b0 = _V[f[i]][_B];
      a1 = _V[f[(i+1) % 3]][_A];
      b1 = _V[f[(i+1) % 3]][_B];
      da = a1 - a0;
      db = b1 - b0;
      a0_2 = a0 * a0;
      a0_3 = a0_2 * a0;
      a0_4 = a0_3 * a0;
      b0_2 = b0 * b0;
      b0_3 = b0_2 * b0;
      b0_4 = b0_3 * b0;
      a1_2 = a1 * a1;
      a1_3 = a1_2 * a1;
      b1_2 = b1 * b1;
      b1_3 = b1_2 * b1;

      C1 = a1 + a0;
      Ca = a1*C1 + a0_2;
      Caa = a1*Ca + a0_3;
      Caaa = a1*Caa + a0_4;
      Cb = b1*(b1 + b0) + b0_2;
      Cbb = b1*Cb + b0_3;
      Cbbb = b1*Cbb + b0_4;
      Cab = 3*a1_2 + 2*a1*a0 + a0_2;
      Kab = a1_2 + 2*a1*a0 + 3*a0_2;
      Caab = a0*Cab + 4*a1_3;
      Kaab = a1*Kab + 4*a0_3;
      Cabb = 4*b1_3 + 3*b1_2*b0 + 2*b1*b0_2 + b0_3;
      Kabb = b1_3 + 2*b1_2*b0 + 3*b1*b0_2 + 4*b0_3;

      _P1 += db*C1;
      _Pa += db*Ca;
      _Paa += db*Caa;
      _Paaa += db*Caaa;
      _Pb += da*Cb;
      _Pbb += da*Cbb;
      _Pbbb += da*Cbbb;
      _Pab += db*(b1*Cab + b0*Kab);
      _Paab += db*(b1*Caab + b0*Kaab);
      _Pabb += da*(a1*Cabb + a0*Kabb);
    }

    _P1 /= 2.0;
    _Pa /= 6.0;
    _Paa /= 12.0;
    _Paaa /= 20.0;
    _Pb /= -6.0;
    _Pbb /= -12.0;
    _Pbbb /= -20.0;
    _Pab /= 24.0;
    _Paab /= 60.0;
    _Pabb /= -60.0;
  }
  void compFaceIntegrals(const Vec3i& f,const Vec3& n,const scalar& w) {
    scalar k1, k2, k3, k4;

    compProjectionIntegrals(f);

    k1 = 1 / n[_C];
    k2 = k1 * k1;
    k3 = k2 * k1;
    k4 = k3 * k1;

    _Fa = k1 * _Pa;
    _Fb = k1 * _Pb;
    _Fc = -k2 * (n[_A]*_Pa + n[_B]*_Pb + w*_P1);

    _Faa = k1 * _Paa;
    _Fbb = k1 * _Pbb;
    _Fcc = k3 * (SQR(n[_A])*_Paa + 2*n[_A]*n[_B]*_Pab + SQR(n[_B])*_Pbb
                 + w*(2*(n[_A]*_Pa + n[_B]*_Pb) + w*_P1));

    _Faaa = k1 * _Paaa;
    _Fbbb = k1 * _Pbbb;
    _Fccc = -k4 * (CUBE(n[_A])*_Paaa + 3*SQR(n[_A])*n[_B]*_Paab
                   + 3*n[_A]*SQR(n[_B])*_Pabb + CUBE(n[_B])*_Pbbb
                   + 3*w*(SQR(n[_A])*_Paa + 2*n[_A]*n[_B]*_Pab + SQR(n[_B])*_Pbb)
                   + w*w*(3*(n[_A]*_Pa + n[_B]*_Pb) + w*_P1));

    _Faab = k1 * _Paab;
    _Fbbc = -k2 * (n[_A]*_Pabb + n[_B]*_Pbbb + w*_Pbb);
    _Fcca = k3 * (SQR(n[_A])*_Paaa + 2*n[_A]*n[_B]*_Paab + SQR(n[_B])*_Pabb
                  + w*(2*(n[_A]*_Paa + n[_B]*_Pab) + w*_Pa));
  }
  void compVolumeIntegrals() {
    Vec3i f;
    Vec3 n;
    sizeType i;

    _T0 = _T1[X] = _T1[Y] = _T1[Z]
                            = _T2[X] = _T2[Y] = _T2[Z]
                                       = _TP[X] = _TP[Y] = _TP[Z] = 0;

    for (i = 0; i < (sizeType)_I.size(); i++) {
      f = _I[i];
      n=_TN[i];
      scalar nx=std::abs(n[0]);
      scalar ny=std::abs(n[1]);
      scalar nz=std::abs(n[2]);
      if (nx > ny && nx > nz) _C = X;
      else _C = (ny > nz) ? Y : Z;
      _A = (_C + 1) % 3;
      _B = (_A + 1) % 3;

      compFaceIntegrals(f,n,-n.dot(_V[f[0]]));

      _T0 += n[X] * ((_A == X) ? _Fa : ((_B == X) ? _Fb : _Fc));

      _T1[_A] += n[_A] * _Faa;
      _T1[_B] += n[_B] * _Fbb;
      _T1[_C] += n[_C] * _Fcc;
      _T2[_A] += n[_A] * _Faaa;
      _T2[_B] += n[_B] * _Fbbb;
      _T2[_C] += n[_C] * _Fccc;
      _TP[_A] += n[_A] * _Faab;
      _TP[_B] += n[_B] * _Fbbc;
      _TP[_C] += n[_C] * _Fcca;
    }

    _T1[X] /= 2;
    _T1[Y] /= 2;
    _T1[Z] /= 2;
    _T2[X] /= 3;
    _T2[Y] /= 3;
    _T2[Z] /= 3;
    _TP[X] /= 2;
    _TP[Y] /= 2;
    _TP[Z] /= 2;
  }
  Mat3 getMCCT() const {
    Mat3 ret=Mat3::Identity();
    ret(0,0)=_T2[X];
    ret(1,1)=_T2[Y];
    ret(2,2)=_T2[Z];
    ret(1,0)=ret(0,1)=_TP[X];
    ret(2,1)=ret(1,2)=_TP[Y];
    ret(2,0)=ret(0,2)=_TP[Z];
    return ret;
  }
private:
  sizeType _A;
  sizeType _B;
  sizeType _C;
  scalar _P1, _Pa, _Pb, _Paa, _Pab, _Pbb, _Paaa, _Paab, _Pabb, _Pbbb;
  scalar _Fa, _Fb, _Fc, _Faa, _Fbb, _Fcc, _Faaa, _Fbbb, _Fccc, _Faab, _Fbbc, _Fcca;
  scalar _T0, _T1[3], _T2[3], _TP[3];
  POS _TN;
};
RigidBodyMass::RigidBodyMass(const ObjMesh& mesh)
{
  if(mesh.getDim() == 2) {
    RigidBodyMass2D mass(mesh.getV(),mesh.getI());
    _ctr=mass.getCtr();
    _mat=mass.getMass();
    _matCOM=mass.getMassCOM();
    _MCCT=mass.getMCCT();
  } else {
    RigidBodyMass3D mass(mesh.getV(),mesh.getI());
    _ctr=mass.getCtr();
    _mat=mass.getMass();
    _matCOM=mass.getMassCOM();
    _MCCT=mass.getMCCT();
  }
}
Vec3 RigidBodyMass::getCtr() const
{
  return _ctr;
}
Mat6 RigidBodyMass::getMass() const
{
  return _mat;
}
Mat6 RigidBodyMass::getMassCOM() const
{
  return _matCOM;
}
scalar RigidBodyMass::getM() const
{
  return _mat(0,0);
}
Vec3 RigidBodyMass::getMC() const
{
  return invCross<scalar>(_mat.block<3,3>(3,0));
}
Mat3 RigidBodyMass::getMCCT() const
{
  return _MCCT;
}
