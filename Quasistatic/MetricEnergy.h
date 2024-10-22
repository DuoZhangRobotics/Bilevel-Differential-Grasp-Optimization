#ifndef METRIC_ENERGY_H
#define METRIC_ENERGY_H

#include "ArticulatedObjective.h"

PRJ_BEGIN

enum METRIC_TYPE
{
  Q_1,
  Q_INF,
  Q_INF_CONSTRAINT,
  Q_INF_CONSTRAINT_FGT,
  Q_INF_BARRIER,
  NO_METRIC,
};
enum METRIC_ACTIVATION
{
  EXP_ACTIVATION,
  SQR_EXP_ACTIVATION,
  INVERSE_ACTIVATION,
  SQR_INVERSE_ACTIVATION,
};
template <typename T>
class MetricEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  MetricEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,
               T d0,const T& alpha,T coef,METRIC_TYPE m,METRIC_ACTIVATION a,T normalExtrude=0);
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,STrips* fhess) override;
  T Quality(const Vec& x);
protected:
  void addTerm(T area,const Vec& gw,sizeType linkId,sizeType linkPId,sizeType oid,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  T activation(T param,T* D=NULL,T* DD=NULL) const;
  T _d0,_coef;
  const T& _alpha;
  METRIC_TYPE _type;
  METRIC_ACTIVATION _activation;
  sizeType _off;
  Mat3XT _pss;
};

PRJ_END

#endif
