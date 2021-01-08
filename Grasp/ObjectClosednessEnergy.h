#ifndef OBJECT_CLOSEDNESS_ENERGY_H
#define OBJECT_CLOSEDNESS_ENERGY_H

#include "ArticulatedObjective.h"

PRJ_BEGIN

template <typename T>
class ObjectClosednessEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  ObjectClosednessEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T coef);
  virtual int operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* gFinal,MatT* hFinal) override;
  virtual std::string name() const override;
protected:
  void addTerm(const Vec2i& termId,const PBDArticulatedGradientInfo<T>& info,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  T _coef;
};

PRJ_END

#endif
