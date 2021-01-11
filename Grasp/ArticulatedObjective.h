#ifndef ARTICULATED_OBJECTIVE_H
#define ARTICULATED_OBJECTIVE_H

#include "GraspQualityMetric.h"
#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Utils/ParallelVector.h>

PRJ_BEGIN

template <typename T>
class ArticulatedObjective;
template <typename T>
class GraspPlanner;
template <typename T>
class ArticulatedObjective
{
public:
  DECL_MAP_FUNCS
  DECL_MAP_TYPES_T
  ArticulatedObjective(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj);
  //hand-related objective
  virtual int operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g=NULL,ParallelMatrix<Mat12XT>* h=NULL,Vec* gFinal=NULL,MatT* hFinal=NULL);
  virtual int operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,T& e,Vec* g=NULL,MatT* h=NULL);
  //constraints
  virtual void cons(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType offr,sizeType offc,Vec& c,MatT* cjac=NULL);
  //whether modifying objective function expression is allowed
  virtual void setUpdateCache(bool);
  virtual sizeType nrAdditionalDOF() const;
  virtual sizeType nrCons() const;
  virtual std::string name() const=0;
  void debug(Vec x,T deltaCustom=0);
protected:
  std::vector<KDOP18<scalar>> updateBVH(const PBDArticulatedGradientInfo<T>& info) const;
  const GraspPlanner<T>& _planner;
  const GraspQualityMetric<T>& _obj;
};

PRJ_END

#endif
