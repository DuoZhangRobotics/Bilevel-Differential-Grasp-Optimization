#ifndef GRAVITY_ENERGY_H
#define GRAVITY_ENERGY_H

#include "ArticulatedObjective.h"

PRJ_BEGIN

template <typename T>
class GravityEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  GravityEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const Vec3T& g,T coef);
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,STrips* fhess) override;
protected:
  Vec3T _g;
  T _coef;
};

PRJ_END

#endif
