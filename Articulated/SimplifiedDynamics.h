#ifndef SIMPLIFIED_DYNAMICS_H
#define SIMPLIFIED_DYNAMICS_H

#include "PBDArticulatedGradientInfo.h"
#include <Optimizer/DSSQPObjective.h>
#include <Utils/Hash.h>

PRJ_BEGIN

struct EndEffectorBounds : public SerializableBase
{
  EndEffectorBounds();
  EndEffectorBounds(sizeType jid);
  EndEffectorBounds(sizeType jid,const Vec3d& localPos,scalarD phi0);
  bool exists(const BBox<sizeType>& bb) const;
  bool isBetterBB(const BBox<sizeType>& bb,const std::vector<sizeType>& sym) const;
  bool sameOrthant(const Vec3i& id,const Vec3d& pos,const std::vector<sizeType>& sym) const;
  static BBox<sizeType> closeFar(const BBox<sizeType>& in,const std::vector<sizeType>& sym);
  sizeType nrDOF(const ArticulatedBody& body) const;
  Cold getDOF(const Cold& pt,const Vec3d& ctr) const;
  Vec3d globalPosAll(const ArticulatedBody& body,const Cold& x) const;
  Vec3d globalPos(const ArticulatedBody& body,const Cold& x) const;
  Vec3d randomPos(const Vec3d& ctr) const;
  sizeType jointId() const;
  bool operator==(const EndEffectorBounds& other) const;
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  //data
  std::unordered_map<Vec3i,Cold,Hash> _indices;
  std::vector<sizeType> _JID;
  BBox<sizeType> _bb;
  Vec3d _localPos;
  scalarD _phi0;
  scalarD _res;
};
#ifdef OPTIMIZER_SUPPORT
template <typename T>
class InverseKinematicsObj : public DSSQPObjectiveComponent<T>
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  using DSSQPObjectiveComponent<T>::_gl;
  using DSSQPObjectiveComponent<T>::_gu;
  using DSSQPObjectiveComponent<T>::_offset;
  InverseKinematicsObj(DSSQPObjectiveCompound<T>& obj,const ArticulatedBody& body,const EndEffectorBounds& ee,sizeType off,bool allDOF=false,const Vec* init=NULL);
  int operator()(const Vec& x,Vec& fvec,STrips* fjac) override;
  static Vec assignDOF(const ArticulatedBody& body,const EndEffectorBounds& ee,const Vec& x,const std::vector<sizeType>* vid);
  static void assignDOF(const ArticulatedBody& body,const EndEffectorBounds& ee,STrips& fjac,const Vec& G,sizeType offset,const std::vector<sizeType>& vid);
  static std::string DOF(sizeType r,sizeType c);
  Vec3T _target;
  Mat3T _H;
protected:
  PBDArticulatedGradientInfo<T> _info;
  const ArticulatedBody& _body;
  const EndEffectorBounds& _ee;
  std::vector<sizeType> _vid;
};
template <typename T>
class VariationalSimilarityObj : public DSSQPObjectiveComponent<T>
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  using DSSQPObjectiveComponent<T>::_gl;
  using DSSQPObjectiveComponent<T>::_gu;
  using DSSQPObjectiveComponent<T>::_offset;
  VariationalSimilarityObj(DSSQPObjectiveCompound<T>& obj,const ArticulatedBody& body,const Vec& init,sizeType off,bool rotBase,bool transBase);
  T operator()(const Vec& x,Vec* fgrad) override;
  PBDArticulatedGradientInfo<T> _info[2];
protected:
  std::vector<sizeType> _vid,_vidBase;
  const ArticulatedBody& _body;
};
#endif
class SimplifiedDynamics : public SerializableBase
{
public:
  SimplifiedDynamics();
  SimplifiedDynamics(std::shared_ptr<ArticulatedBody> body,scalarD res=1,sizeType symDir=-1,const Vec3d& zRange=Vec3d::Zero(),scalarD zMargin=0);
  static void detectEndEffector(const ArticulatedBody& body,sizeType jid,Vec3d& localPos,scalarD& phi0,const Vec3d& zRange=Vec3d::Zero());
  static Cold inverseKinematics(const ArticulatedBody& body,const std::vector<EndEffectorBounds>& ee,
                                const std::vector<Vec3d,Eigen::aligned_allocator<Vec3d>>& pos,
                                const Mat3d& H,const Cold* init,bool rotBase=false,bool transBase=true);
  Cold inverseKinematics(const EndEffectorBounds& ee,const Vec3d& pos,const Cold* init) const;
  PBDArticulatedGradientInfo<scalarD> getDOF(const Vec6d& base,const Mat3Xd& xss) const;
  PBDArticulatedGradientInfo<scalarD> getDOF(const Mat3X4d& base,const Mat3Xd& xss) const;
  void debugDOF(const std::string& path,sizeType nrPose=10) const;
  void optimizeBB(EndEffectorBounds& ee,const Vec3d& zRange,scalarD zMargin) const;
  void initializeBB(EndEffectorBounds& ee) const;
  void detectSymmetry();
  void cullSymmetry(sizeType symDir);
  void cullSymmetry(sizeType i,sizeType symDir);
  void writeVTK(const std::string& path) const;
  std::shared_ptr<ArticulatedBody> getBody() const;
  EndEffectorBounds getForceInfo(sizeType jid) const;
  const std::vector<EndEffectorBounds>& getForceInfos() const;
  Vec3d lb(sizeType eeid) const;
  Vec3d ub(sizeType eeid) const;
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
protected:
  std::shared_ptr<ArticulatedBody> _body;
  std::shared_ptr<ArticulatedBody> _bodySimplified;
  std::vector<EndEffectorBounds> _ee;
  std::vector<sizeType> _sym;
  sizeType _rootId;
  Vec3d _ctr;
};

PRJ_END

#endif
