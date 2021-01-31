#ifndef CONTACT_INDEX_CONSTRAINT_H
#define CONTACT_INDEX_CONSTRAINT_H

#include "ArticulatedObjective.h"

PRJ_BEGIN

template <typename T>
class ContactIndexConstraint : public ArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using ArticulatedObjective<T>::_planner;
  using ArticulatedObjective<T>::_object;
  using ArticulatedObjective<T>::_info;
  using ArticulatedObjective<T>::_offset;
  struct ContactIndexMap
  {
    bool operator<(const ContactIndexMap& other) const;
    bool operator!=(const ContactIndexMap& other) const;
    bool operator==(const ContactIndexMap& other) const;
    sizeType _oid,_oidOther,_pid;
    Mat3XT _objDir;
    Vec3T _objPos;
    Coli _off;
  };
  ContactIndexConstraint(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const Vec3T& g,const T& tau);
  void addContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid,T mu,sizeType nDir);
  bool existContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid) const;
  void clearContactIndexMap();
  void writeContactVTK(const Vec& x,const PBDArticulatedGradientInfo<T>& info,const std::string& path) const;
  //constraints
  virtual int operator()(const Vec& x,Vec& fvec,STrips* fjac=NULL) override;
  virtual int values() const override;
protected:
  std::vector<ContactIndexMap> _cmss;
  DSSQPObjectiveCompound<T>& _obj;
  const Vec3T _g;
  const T& _tau;
};

PRJ_END

#endif
