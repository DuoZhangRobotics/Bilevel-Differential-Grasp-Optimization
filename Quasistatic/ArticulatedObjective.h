#ifndef ARTICULATED_OBJECTIVE_H
#define ARTICULATED_OBJECTIVE_H

#include <Optimizer/DSSQPObjective.h>
#include <Utils/ParallelVector.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct PBDArticulatedGradientInfo;
template <typename T>
class PointCloudObject;
template <typename T>
class GraspPlanner;
template <typename T>
class ArticulatedObjective : public DSSQPObjectiveComponent<T>
{
public:
  DECL_MAP_FUNCS
  DECL_MAP_TYPES_T
  using DSSQPObjectiveComponent<T>::inputs;
  using DSSQPObjectiveComponent<T>::values;
  using DSSQPObjectiveComponent<T>::_name;
  using DSSQPObjectiveComponent<T>::_tmp;
  using DSSQPObjectiveComponent<T>::setUpdateCache;
  ArticulatedObjective(DSSQPObjectiveCompound<T>& obj,const std::string& name,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object);
  //hand-related objective
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g=NULL,ParallelMatrix<Mat12XT>* h=NULL,Vec* fgrad=NULL,STrips* fhess=NULL);
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g=NULL,ParallelMatrix<Mat12XT>* h=NULL,Vec* fgrad=NULL,SMat* fhess=NULL);
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g=NULL,ParallelMatrix<Mat12XT>* h=NULL,Vec* fgrad=NULL,DMat* fhess=NULL);
  virtual T operator()(const Vec& x,Vec* fgrad=NULL,STrips* fhess=NULL) override;
  //whether modifying objective function expression is allowed
  virtual void setUpdateCache(const Vec& x,bool) override;
  const PBDArticulatedGradientInfo<T>& info() const;
  PBDArticulatedGradientInfo<T>& info();
  const PointCloudObject<T>& object() const;
protected:
  std::vector<KDOP18<scalar>> updateBVH() const;
  const GraspPlanner<T>& _planner;
  const PointCloudObject<T>& _object;
  const PBDArticulatedGradientInfo<T>& _info;
};

PRJ_END

#endif
