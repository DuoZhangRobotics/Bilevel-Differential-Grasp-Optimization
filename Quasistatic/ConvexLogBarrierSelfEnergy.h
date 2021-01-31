#ifndef CONVEX_LOG_BARRIER_SELF_ENERGY_H
#define CONVEX_LOG_BARRIER_SELF_ENERGY_H

#include "ArticulatedObjective.h"
#include <CommonFile/Hash.h>
#include <unordered_set>

PRJ_BEGIN

template <typename T>
class ConvexLogBarrierSelfEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::updateBVH;
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  struct SeparatingPlane
  {
    Mat3XT _pss[2];
    Vec4T _plane;
  };
  ConvexLogBarrierSelfEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T d0,T mu,bool allPairs=false);
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,STrips* fhess) override;
  void updatePlanes();
  virtual void setUpdateCache(const Vec& x,bool update) override;
protected:
  bool initializePlane(sizeType idL,sizeType idR);
  Vec4T updatePlane(const Vec2i& linkId,const SeparatingPlane& sp,bool callback) const;
  void addTerm(bool& valid,const std::tuple<Vec2i,sizeType,sizeType>& termId,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  std::unordered_map<Vec2i,SeparatingPlane,Hash> _plane;
  std::unordered_set<Vec2i,Hash> _exclude;
  bool _updateCache;
  bool _allPairs;
  T _d0,_mu;
};

PRJ_END

#endif
