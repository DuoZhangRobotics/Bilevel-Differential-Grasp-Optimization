#ifndef LOG_BARRIER_SELF_ENERGY_H
#define LOG_BARRIER_SELF_ENERGY_H

#include "ArticulatedObjective.h"
#include <Utils/Hash.h>
#include <unordered_set>

PRJ_BEGIN

template <typename T>
class LogBarrierSelfEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::updateBVH;
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_obj;
  struct SeparatingPlane
  {
    Mat3XT _pss[2];
    Vec4T _plane;
  };
  LogBarrierSelfEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T d0,T mu,bool allPairs=false);
  virtual int operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* gFinal,MatT* hFinal) override;
  void updatePlanes(const PBDArticulatedGradientInfo<T>& info);
  virtual void setUpdateCache(bool update) override;
  virtual std::string name() const override;
protected:
  bool initializePlane(const PBDArticulatedGradientInfo<T>& info,sizeType idL,sizeType idR);
  Vec4T updatePlane(const PBDArticulatedGradientInfo<T>& info,const Vec2i& linkId,const SeparatingPlane& sp,bool callback) const;
  void addTerm(bool& valid,const std::tuple<Vec2i,sizeType,sizeType>& termId,const PBDArticulatedGradientInfo<T>& info,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  std::unordered_map<Vec2i,SeparatingPlane,Hash> _plane;
  std::unordered_set<Vec2i,Hash> _exclude;
  bool _updateCache;
  bool _allPairs;
  T _d0,_mu;
};

PRJ_END

#endif
