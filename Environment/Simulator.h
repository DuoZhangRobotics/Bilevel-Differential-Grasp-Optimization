#ifdef OPTIMIZER_SUPPORT
#ifdef ENVIRONMENT_SUPPORT
#ifdef DEFORMABLE_SUPPORT
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <CommonFile/geom/BVHNode.h>
#include <CommonFile/geom/SweepAndPrune.h>
#include <Deformable/FEMGradientInfo.h>
#include <Articulated/PBDArticulatedGradientInfo.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct CollisionConstraint;
template <typename T>
class CollisionNode;
template <typename T>
class PBDSimulator;
template <typename T>
struct SimulatorObject
{
  DECL_MAP_TYPES_T
  sizeType nrControl() const;
  sizeType nrDOF() const;
  std::shared_ptr<FEMSystem<T>> _sysFEM;//soft
  FEMGradientInfo<T> _xFEM[3];
  std::shared_ptr<ArticulatedBody> _body;//rigid
  std::shared_ptr<PBDSimulator<T>> _sysPBD;
  PBDArticulatedGradientInfo<T> _xPBD[3];
  sizeType _offControl,_offDOF;
  //this is a convenient way to attach two objects A and B
  //when you apply a transform to object A, all attached objects will be updated
  std::unordered_map<sizeType,Mat3X4T> _attachedObjects;
};
template <typename T>
struct SimulatorState
{
  FEMGradientInfo<T> _xFEM;
  PBDArticulatedGradientInfo<T> _xPBD;
};
template <typename T>
class Simulator
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef ParallelMatrix<T> PVal;
  typedef ParallelMatrix<Vec> PVec;
  typedef ParallelVector<Eigen::Triplet<T,sizeType> > TRIPS;
  Simulator();
  Vec getState(sizeType id=0) const;
  void setState(const Vec& x,sizeType id=0);
  sizeType addFEMObject(std::shared_ptr<FEMSystem<T>> sys,const Vec& xInit);
  sizeType addPBDObject(std::shared_ptr<ArticulatedBody> body,const Vec3T& g,const Vec& xInit);
  void writeConstraintsVTK(const std::string& path) const;
  void writeVTK(bool embed,const std::string& path) const;
  void writeObj(bool embed,ObjMesh& mesh) const;
  bool step(const Vec& ctrl,const std::string& callbackPrefix);
  void attach(sizeType oidA,sizeType oidB); //attach two objects in current pose
  void setTrans(sizeType oid,sizeType id,const Mat3X4T& t);
  Mat3X4T getTrans(sizeType oid,sizeType id=0) const;
  void debugCollisionConstraint(sizeType N,sizeType nrIter,T scale);
  void debug(sizeType nrIter,T scale);
  //parameter
  sizeType nrControl() const;
  sizeType nrDOF() const;
  const T& dtMin() const;
  T& dtMin();
  const T& dtMax() const;
  T& dtMax();
  const T& dist() const;
  T& dist();
  const T& dThres() const;
  T& dThres();
  const T& rThres() const;
  T& rThres();
  const T& aThres() const;
  T& aThres();
  const sizeType& itMax() const;
  sizeType& itMax();
protected:
  void advance();
  bool stepInnerAdaptive(T dt,const Vec& ctrl,const std::string& callbackPrefix);
  bool stepInner(const Vec& ctrl,const std::string& callbackPrefix);
  std::unordered_set<sizeType> attachedChildren(sizeType root) const;
  void saveState(std::vector<SimulatorState<T>>& s,sizeType id=0) const;
  void loadState(const std::vector<SimulatorState<T>>& s,sizeType id=0);
  void updateState(std::vector<SimulatorState<T>>& s,const Vec& x) const;
  bool buildSystem(const Vec& ctrl,T* E,Vec* G,SMat* H,SMat* diagonalPerturb,bool pad=true) const;
  //collision
  Vec getCollisionConstraintJacobian(SMat* J);
  void getCollisionConstraint(sizeType id=0);
  bool checkCollision(sizeType id=0);
  bool filterCollision(sizeType i,sizeType j) const;
  //data
  std::vector<CollisionConstraint<T>> _constraints;
  std::vector<std::shared_ptr<CollisionNode<T>>> _css;
  SweepAndPruneIncremental<scalarD,3> _sap;
  std::vector<SimulatorObject<T>> _oss;
  T _dt,_lastDt,_dtMin,_dtMax,_dist;
  T _dThres,_rThres,_aThres;
  sizeType _itMax,_fDir;
  Vec _lb,_ub;
};

PRJ_END

#endif
#endif
#endif
#endif
