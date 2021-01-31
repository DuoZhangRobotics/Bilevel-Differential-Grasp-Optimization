#ifndef PRIMAL_DUAL_QINF_METRIC_ENERGY_FGT_H
#define PRIMAL_DUAL_QINF_METRIC_ENERGY_FGT_H

#include "PrimalDualQInfMetricEnergy.h"

PRJ_BEGIN

template <typename T>
struct FGTTreeNode;
template <typename T>
class PrimalDualQInfMetricEnergyFGT : public PrimalDualQInfMetricEnergy<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  using MetricEnergy<T>::_alpha;
  using MetricEnergy<T>::_coef;
  using MetricEnergy<T>::_pss;
  using PrimalDualQInfMetricEnergy<T>::values;
  PrimalDualQInfMetricEnergyFGT(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const T& alpha,T coef,T normalExtrude=0,T FGTThres=1e-6f);
  //constraints
  virtual int operator()(const Vec& x,Vec& fvec,STrips* fjac=NULL) override;
protected:
  std::vector<std::shared_ptr<FGTTreeNode<T>>> _gripperFGT;
  std::shared_ptr<FGTTreeNode<T>> _objectFGT;
  T _FGTThres;
};

PRJ_END

#endif
