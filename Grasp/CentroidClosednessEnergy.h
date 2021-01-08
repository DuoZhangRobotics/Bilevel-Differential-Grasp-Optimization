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
  using ArticulatedObjective<T>::_obj;
  CentroidClosednessEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T coef);
  virtual int operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* gFinal,MatT* hFinal) override;
  virtual std::string name() const override;
protected:
  void addTerm(sizeType i,const PBDArticulatedGradientInfo<T>& info,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> _centroid;
  T _coef;
};

PRJ_END

#endif
