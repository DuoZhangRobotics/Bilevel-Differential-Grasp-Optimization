#ifndef LOG_BARRIER_OBJ_ENERGY_H
#define LOG_BARRIER_OBJ_ENERGY_H

#include "ArticulatedObjective.h"
#include <CommonFile/Hash.h>
#include <unordered_map>

PRJ_BEGIN

template <typename T>
class LogBarrierObjEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::updateBVH;
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  LogBarrierObjEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T d0,T mu,const bool& useGJK);
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,STrips* fhess) override;
  virtual void setUpdateCache(const Vec& x,bool update) override;
protected:
  void addTerm(bool& valid,const Vec2i& termId,Vec2i& feat,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  std::unordered_map<Vec2i,Vec2i,Hash> _cache;
  const bool& _useGJK;
  bool _updateCache;
  T _d0,_d1,_mu;
};

PRJ_END

#endif
