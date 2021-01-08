#ifndef LOG_BARRIER_OBJ_ENERGY_H
#define LOG_BARRIER_OBJ_ENERGY_H

#include "ArticulatedObjective.h"
#include <Utils/Hash.h>
#include <unordered_map>

PRJ_BEGIN

template <typename T>
class LogBarrierObjEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::updateBVH;
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  LogBarrierObjEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T d0,T mu,const bool& useGJK);
  virtual int operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* gFinal,MatT* hFinal) override;
  virtual void setUpdateCache(bool update) override;
  virtual std::string name() const override;
protected:
  void addTerm(bool& valid,const Vec2i& termId,Vec2i& feat,const PBDArticulatedGradientInfo<T>& info,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  std::unordered_map<Vec2i,Vec2i,Hash> _cache;
  const bool& _useGJK;
  bool _updateCache;
  T _d0,_d1,_mu;
};

PRJ_END

#endif
