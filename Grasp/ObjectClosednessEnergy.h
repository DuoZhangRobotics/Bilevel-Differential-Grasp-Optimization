#ifndef OBJECT_CLOSEDNESS_ENERGY_H
#define OBJECT_CLOSEDNESS_ENERGY_H

#include "GraspPlanner.h"

PRJ_BEGIN

template <typename T>
class ObjectClosednessEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  ObjectClosednessEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T coef);
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h) override;
protected:
  T _coef;
};

PRJ_END

#endif
