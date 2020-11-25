#ifndef IMPLICIT_FUNC_H
#define IMPLICIT_FUNC_H

#include "GridBasic.h"
#include "ObjMesh.h"
#include "ImplicitFuncInterface.h"
#include "CollisionDetection.h"

PRJ_BEGIN

class ImplicitFuncCSG : public ImplicitFunc<scalar>
{
public:
  enum OP_TYPE {
    INTERSECT,
    UNION,
    SUBTRACT,
  };
  ImplicitFuncCSG(OP_TYPE op);
  ImplicitFuncCSG(const BBox<scalar,3>& bb,OP_TYPE op);
  ImplicitFuncCSG(const BBox<scalar,2>& bb,OP_TYPE op);
  virtual scalar operator()(const Vec3& pos) const;
  void setAlpha(const scalar& alpha);
  std::shared_ptr<ImplicitFunc<scalar> > _a;
  std::shared_ptr<ImplicitFunc<scalar> > _b;
  scalar _alpha;
  OP_TYPE _op;
};
class ImplicitFuncPlane : public ImplicitFunc<scalar>
{
public:
  ImplicitFuncPlane();
  ImplicitFuncPlane(const Vec3& x0,const Vec3& n);
  ImplicitFuncPlane(const Vec3& a,const Vec3& b,const Vec3& c);
  virtual scalar operator()(const Vec3& pos) const;
  PlaneTpl<scalar> _p;
};
class ImplicitFuncGridRef : public ImplicitFunc<scalar>
{
public:
  ImplicitFuncGridRef(const Grid<scalar,scalar>& ls):_lsRef(ls) {}
  virtual BBox<scalar> getBB() const {
    return _lsRef.getBB();
  }
  virtual scalar operator()(const Vec3& pos) const;
  const Grid<scalar,scalar>& _lsRef;
};
class ImplicitFuncGrid : public ImplicitFuncGridRef
{
public:
  ImplicitFuncGrid():ImplicitFuncGridRef(_ls) {}
  virtual BBox<scalar> getBB() const {
    return _ls.getBB();
  }
  Grid<scalar,scalar> _ls;
};
class ImplicitFuncReinit : public ImplicitFuncGrid
{
public:
  ImplicitFuncReinit(scalar cellSz,const ImplicitFunc<scalar>& inner);
  ImplicitFuncReinit(scalar cellSz,ImplicitFunc<scalar>& inner,bool verbose);
  ImplicitFuncReinit(const Grid<scalar,scalar> &tpl,const ImplicitFunc<scalar>& inner);
  virtual scalar operator()(const Vec3& pos) const;
};
struct ObjMeshGeomCell;
class ImplicitFuncMeshRef : public ImplicitFunc<scalar>
{
public:
  ImplicitFuncMeshRef(const ObjMeshGeomCell& mesh,scalar eps=0.0f);
  virtual void beginSampleSet(const ScalarField& toBeSampled);
  virtual void endSampleSet(const ScalarField& toBeSampled);
  virtual scalar operator()(const Vec3& pos) const;
  virtual BBox<scalar> getBB() const;
protected:
  void search(Vec3i base,sizeType dir);
  std::shared_ptr<TagField> _tag;
  const ObjMeshGeomCell& _mesh;
  const scalar _eps;
  const Vec3 _ext;
};
class ImplicitFuncOffset : public ImplicitFunc<scalar>
{
public:
  ImplicitFuncOffset():_off(0) {}
  virtual void beginSampleSet(const ScalarField& toBeSampled) {
    _inner->beginSampleSet(toBeSampled);
  }
  virtual void endSampleSet(const ScalarField& toBeSampled) {
    _inner->endSampleSet(toBeSampled);
  }
  scalar operator()(const Vec3& pos) const {
    return _inner->operator()(pos)+_off;
  }
  virtual BBox<scalar> getBB() const {
    BBox<scalar> bb=_inner->getBB();
    bb.enlarged(std::abs(_off),bb.getExtent()[2] == 0 ? 2 : 3);
    return bb;
  }
  //data
  scalar _off;
  std::shared_ptr<ImplicitFunc<scalar> > _inner;
};
class ImplicitFuncNegate : public ImplicitFunc<scalar>
{
public:
  virtual void beginSampleSet(const ScalarField& toBeSampled) {
    _inner->beginSampleSet(toBeSampled);
  }
  virtual void endSampleSet(const ScalarField& toBeSampled) {
    _inner->endSampleSet(toBeSampled);
  }
  scalar operator()(const Vec3& pos) const {
    return -_inner->operator()(pos);
  }
  virtual BBox<scalar> getBB() const {
    return _inner->getBB();
  }
  //data
  std::shared_ptr<ImplicitFunc<scalar> > _inner;
};
class ImplicitFuncRosy : public ImplicitFunc<scalar>
{
public:
  ImplicitFuncRosy(const Vec3& origin,const Vec3& X,const scalar& step,const scalar& coef);
  virtual scalar operator()(const Vec3& pos) const;
  virtual scalar dist(const Vec2& p) const;
  virtual scalar y(const scalar& x) const =0;
  //axis
  scalar _step;
  scalar _coef;
  Vec3 _origin;
  Vec3 _X;
};

PRJ_END

#endif
