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
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  using MetricEnergy<T>::_alpha;
  using MetricEnergy<T>::_coef;
  using MetricEnergy<T>::_pss;
  PrimalDualQInfMetricEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const T& alpha,T coef,METRIC_ACTIVATION a,T normalExtrude=0);
  //additional objective
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,STrips* fhess) override;
  //constraints
  virtual int operator()(const Vec& x,Vec& fvec,STrips* fjac=NULL) override;
  virtual int values() const override;
protected:
  T addTerm(T area,sizeType linkId,sizeType oid,Mat3X4TM cjacG) const;
};

PRJ_END

#endif
