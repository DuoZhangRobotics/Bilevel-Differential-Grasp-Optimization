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
  using ArticulatedObjective<T>::_obj;
  using MetricEnergy<T>::_alpha;
  using MetricEnergy<T>::_coef;
  using MetricEnergy<T>::_pss;
  using PrimalDualQInfMetricEnergy<T>::nrCons;
  PrimalDualQInfMetricEnergyFGT(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,const T& alpha,T coef,T normalExtrude=0,T FGTThres=1e-6f);
  //constraints
  virtual void cons(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType offr,sizeType offc,Vec& c,MatT* cjac=NULL) override;
  virtual std::string name() const override;
protected:
  std::vector<std::shared_ptr<FGTTreeNode<T>>> _gripperFGT;
  std::shared_ptr<FGTTreeNode<T>> _objectFGT;
  T _FGTThres;
};

PRJ_END

#endif
