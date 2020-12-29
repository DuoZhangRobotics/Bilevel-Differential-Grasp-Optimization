#ifndef METRIC_ENERGY_H
#define METRIC_ENERGY_H

#include "GraspPlanner.h"

PRJ_BEGIN

template <typename T>
class MetricEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::updateBVH;
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  MetricEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T alpha,METRIC_TYPE m,T coef);
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h) override;
protected:
  T activation(T param,T* D=NULL,T* DD=NULL) const;
  T _alpha,_coef;
  METRIC_TYPE _type;
};

PRJ_END

#endif
