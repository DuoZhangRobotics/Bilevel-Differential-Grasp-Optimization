#ifndef COLLISION_FUNC_H
#define COLLISION_FUNC_H

#include "ParticleSet.h"
#include <set>

PRJ_BEGIN

template<typename P_TYPE>
class CollisionFunc
{
public:
  virtual ~CollisionFunc() {}
  virtual void operator()(const P_TYPE& p,const sizeType& id) =0;
};
template<typename P_TYPE>
class CollisionFuncEarlyStop
{
public:
  virtual ~CollisionFuncEarlyStop() {}
  virtual bool operator()(const P_TYPE& p,const sizeType& id) =0;
};
template<typename P_TYPE>
class CollisionFuncWrapper : public CollisionFuncEarlyStop<P_TYPE>
{
public:
  CollisionFuncWrapper(CollisionFunc<P_TYPE>& wrapped):_wrapped(wrapped) {}
  virtual bool operator()(const P_TYPE& p,const sizeType& id) {
    _wrapped(p,id);
    return false;
  }
  CollisionFunc<P_TYPE>& _wrapped;
};
template<typename P_TYPE>
class CollisionFuncRecord : public CollisionFunc<P_TYPE>
{
public:
  CollisionFuncRecord() {}
  virtual void operator()(const P_TYPE& p,const sizeType& id) {
    _neigh.insert(id);
  }
  std::set<sizeType> _neigh;
};
template<typename P_TYPE>
class CollisionFuncAlwaysTrue : public CollisionFuncEarlyStop<P_TYPE>
{
public:
  virtual bool operator()(const P_TYPE& p,const sizeType& id) {
    return true;
  }
};
template<typename VEC3>
struct ExtractPosDirect {
  static FORCE_INLINE const VEC3& extract(const VEC3& p) {
    return p;
  }
  template <typename V>
  static FORCE_INLINE VEC3 assign(VEC3 p,const V& v) {
    p=v;
    return p;
  }
  template <typename V>
  static FORCE_INLINE VEC3 add(VEC3 p,const V& v) {
    p+=v;
    return p;
  }
};
template<typename P_TYPE>
struct ExtractPosParticle {
  typedef typename ScalarUtil<typename P_TYPE::scalarType>::ScalarVec3 Vec3Type;
  static FORCE_INLINE const Vec3Type& extract(const P_TYPE& p) {
    return p._pos;
  }
  template <typename V>
  static FORCE_INLINE P_TYPE assign(P_TYPE p,const V& v) {
    p._pos=v;
    return p;
  }
  template <typename V>
  static FORCE_INLINE P_TYPE add(P_TYPE p,const V& v) {
    p._pos+=v;
    return p;
  }
};

PRJ_END

#endif
