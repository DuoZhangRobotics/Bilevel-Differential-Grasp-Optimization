#ifdef OPTIMIZER_SUPPORT
#ifdef ENVIRONMENT_SUPPORT
#ifdef DEFORMABLE_SUPPORT
#ifndef COLLISION_CONSTRAINT_H
#define COLLISION_CONSTRAINT_H

#include "Environment.h"
#include <Articulated/Joint.h>
#include <CommonFile/geom/BVHNode.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct FEMBody;
template <typename T>
struct FEMGradientInfo;
template <typename T>
struct PBDArticulatedGradientInfo;
template <typename T>
struct SimulatorObject;
template <typename T>
class CollisionNode;
template <typename T>
struct CollisionConstraint
{
  DECL_MAP_TYPES_T
  void buildFrame(sizeType n);
  T assembleSystem(const SimulatorObject<T>& objL,const SimulatorObject<T>& objR,sizeType row,STrips* J) const;
  std::shared_ptr<CollisionNode<T>> _L,_R;
  Mat3XT _frame;
  Vec3d _pL,_pR;
  Vec4d _bary;
  Vec4i _vid;
  bool _isP2T;
};
template <typename T>
class CollisionNode
{
public:
  DECL_MAP_TYPES_T
  typedef KDOP18<scalarD> BBOX;
  CollisionNode(const SimulatorObject<T>& obj,sizeType objId,sizeType jointId);
  void updateBVH(const SimulatorObject<T>& obj,sizeType id);
  //check collision
  static bool checkCollisionE2T(std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R);
  static bool checkCollision(std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R);
  //add constraint
  static void addConstraintP2T(std::vector<CollisionConstraint<T>>& cons,std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R,T dist);
  static void addConstraintE2E(std::vector<CollisionConstraint<T>>& cons,std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R,T dist);
  static void addConstraint(std::vector<CollisionConstraint<T>>& cons,std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R,T dist);
  T assembleSystem(const SimulatorObject<T>& obj,sizeType row,STrips* J,sizeType id,const Vec3T& coef) const;
  BBox<scalarD> getBB() const;
  sizeType objId() const;
  sizeType jointId() const;
  sizeType nrV() const;
protected:
  void updateBVH();
  void buildBVH(const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>>& iss,
                const std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>>& ess);
  std::vector<Vec3d,Eigen::aligned_allocator<Vec3d>> _vlss,_vgss;
  std::vector<Node<sizeType,BBOX>> _bvhV;
  std::vector<Node<Vec2i,BBOX>> _bvhE;
  std::vector<Node<Vec3i,BBOX>> _bvhT;
  sizeType _objId,_jointId;
};

PRJ_END

#endif
#endif
#endif
#endif
