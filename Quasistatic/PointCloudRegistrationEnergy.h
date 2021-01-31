#ifndef POINT_CLOUD_REGISTRATION_ENERGY_H
#define POINT_CLOUD_REGISTRATION_ENERGY_H

#include "ArticulatedObjective.h"

PRJ_BEGIN

template <typename T>
class PointCloudRegistrationEnergy : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  struct PointCloudMap
  {
    bool operator<(const PointCloudMap& other) const;
    bool operator!=(const PointCloudMap& other) const;
    bool operator==(const PointCloudMap& other) const;
    sizeType _pid,_oid;
    Vec3T _objPos;
  };
  PointCloudRegistrationEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T coef);
  void addPointCloudMap(sizeType pid,sizeType oid,const Vec3T& objPos);
  bool existPointCloudMap(sizeType pid,sizeType oid) const;
  void clearPointCloudMap();
  //objective
  virtual int operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,STrips* fhess) override;
protected:
  void addTerm(const PointCloudMap& pm,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const;
  std::vector<PointCloudMap> _pmss;
  T _coef;
};

PRJ_END

#endif
