#ifdef OPTIMIZER_SUPPORT
#ifdef ENVIRONMENT_SUPPORT
#ifdef DEFORMABLE_SUPPORT
#ifndef SIMULATOR_H
#define SIMULATOR_H

#include <Deformable/FEMGradientInfo.h>
#include <Articulated/PBDArticulatedGradientInfo.h>

PRJ_BEGIN

#include <Utils/ArticulatedBodyPragma.h>
template <typename T>
class PBDSimulator;
template <typename T>
class Simulator
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef ParallelMatrix<T> PVal;
  typedef ParallelMatrix<Vec> PVec;
  typedef ParallelVector<Eigen::Triplet<T,sizeType> > TRIPS;
  struct SimulatorObject
  {
    std::shared_ptr<FEMSystem<T>> _sysFEM;//soft
    FEMGradientInfo<T> _xFEM[3];
    std::shared_ptr<ArticulatedBody> _body;//rigid
    std::shared_ptr<PBDSimulator<T>> _sysPBD;
    PBDArticulatedGradientInfo<T> _xPBD[3];
    //this is a convenient way to attach two objects A and B
    //when you apply a transform to object A, all attached objects will be updated
    std::unordered_map<sizeType,Mat3X4T> _attachedObjects;
  };
  struct SimulatorState
  {
    FEMGradientInfo<T> _xFEM;
    PBDArticulatedGradientInfo<T> _xPBD;
  };
  Simulator();
  sizeType addFEMObject(std::shared_ptr<FEMSystem<T>> sys,const Vec& xInit);
  sizeType addPBDObject(std::shared_ptr<ArticulatedBody> body,const Vec3T& g,const Vec& xInit);
  void writeVTK(bool embed,const std::string& path) const;
  void writeObj(bool embed,ObjMesh& mesh) const;
  bool solve(const Vec& ctrl,T GThres=1e-4f,T alphaThres=1e-4f,sizeType itMax=100,bool direct=true,bool callback=true);
  void attach(sizeType oidA,sizeType oidB); //attach two objects in current pose
  void setTrans(sizeType oid,sizeType id,const Mat3X4T& t);
  Mat3X4T getTrans(sizeType oid,sizeType id=0) const;
  void advance();
  void debug(sizeType nrIter,T scale);
  sizeType nrControl() const;
  sizeType nrDOF() const;
  const T& lastDt() const;
  T& lastDt();
  const T& dt() const;
  T& dt();
protected:
  std::unordered_set<sizeType> attachedChildren(sizeType root) const;
  void saveState(std::vector<SimulatorState>& s,sizeType id=0) const;
  void loadState(const std::vector<SimulatorState>& s,sizeType id=0);
  void updateState(std::vector<SimulatorState>& s,const Vec& x) const;
  bool buildSystem(const Vec& ctrl,T* E,Vec* G,Eigen::SparseMatrix<T,0,sizeType>* H,bool pad=true) const;
  std::vector<SimulatorObject> _oss;
  T _dt,_lastDt;
};

PRJ_END

#endif
#endif
#endif
#endif
