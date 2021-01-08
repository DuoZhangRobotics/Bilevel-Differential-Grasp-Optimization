#ifndef PRIMAL_DUAL_QINF_METRIC_ENERGY_H
#define PRIMAL_DUAL_QINF_METRIC_ENERGY_H

#include "MetricEnergy.h"

PRJ_BEGIN

template <typename T>
class PrimalDualQInfMetricEnergy : public MetricEnergy<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  using MetricEnergy<T>::_alpha;
  using MetricEnergy<T>::_coef;
  using MetricEnergy<T>::_pss;
  PrimalDualQInfMetricEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,const T& alpha,T coef,METRIC_ACTIVATION a,T normalExtrude=0);
  //additional objective
  virtual int operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* gFinal,MatT* hFinal) override;
  //constraints
  virtual void cons(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType offr,sizeType offc,Vec& c,MatT* cjac=NULL) override;
  virtual sizeType nrAdditionalDOF() const override;
  virtual sizeType nrCons() const override;
  virtual std::string name() const override;
protected:
  T addTerm(T area,sizeType linkId,sizeType oid,const PBDArticulatedGradientInfo<T>& info,Mat3X4TM cjacG) const;
};

PRJ_END

#endif
