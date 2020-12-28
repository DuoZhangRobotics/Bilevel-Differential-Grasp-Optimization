#ifndef LOG_BARRIER_OBJ_ENERGY_H
#define LOG_BARRIER_OBJ_ENERGY_H

#include "GraspPlanner.h"

PRJ_BEGIN

template <typename T>
class LogBarrierObjEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::updateBVH;
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  LogBarrierObjEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T d0,T mu);
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h) override;
protected:
  T _d0,_mu;
};

PRJ_END

#endif
