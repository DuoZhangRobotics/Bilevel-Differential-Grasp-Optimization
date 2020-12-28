#ifndef CENTROID_CLOSEDNESS_ENERGY_H
#define CENTROID_CLOSEDNESS_ENERGY_H

#include "GraspPlanner.h"

PRJ_BEGIN

template <typename T>
class CentroidClosednessEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  CentroidClosednessEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T coef);
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h) override;
protected:
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> _centroid;
  T _coef;
};

PRJ_END

#endif
