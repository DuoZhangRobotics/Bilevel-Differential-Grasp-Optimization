#ifndef CENTROID_CLOSEDNESS_ENERGY_H
#define CENTROID_CLOSEDNESS_ENERGY_H

#include "ArticulatedObjective.h"

PRJ_BEGIN

template <typename T>
class CentroidClosednessEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  CentroidClosednessEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T coef);
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,STrips* fhess) override;
protected:
  void addTerm(sizeType i,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> _centroid;
  T _coef;
};

PRJ_END

#endif
