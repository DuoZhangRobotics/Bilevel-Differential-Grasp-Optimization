#ifndef LOG_BARRIER_SELF_ENERGY_H
#define LOG_BARRIER_SELF_ENERGY_H

#include "GraspPlanner.h"
#include <Utils/Hash.h>
#include <unordered_set>

PRJ_BEGIN

template <typename T>
class LogBarrierSelfEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::updateBVH;
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  struct SeparatingPlane
  {
    Mat3XT _pss[2];
    Vec4T _plane;
  };
  LogBarrierSelfEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T d0,T mu,bool allPairs=false);
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h) override;
  void updatePlanes(const PBDArticulatedGradientInfo<T>& info);
protected:
  bool initializePlane(const PBDArticulatedGradientInfo<T>& info,sizeType idL,sizeType idR);
  std::unordered_map<Vec2i,SeparatingPlane,Hash> _plane;
  std::unordered_set<Vec2i,Hash> _exclude;
  bool _allPairs;
  T _d0,_mu;
};

PRJ_END

#endif
