#ifndef PHYSICS_REGISTRATION_H
#define PHYSICS_REGISTRATION_H

#include "GraspPlanner.h"
#include "ContactIndexConstraint.h"
#include "PointCloudRegistrationEnergy.h"

PRJ_BEGIN

struct PhysicsRegistrationParameter
{
  PhysicsRegistrationParameter(Options& ops);
  void reset(Options& ops);
  static void initOptions(PhysicsRegistrationParameter& sol);
  //data
  scalarD _coefPotential;
  scalarD _coefRegister;
  scalarD _tau;
  Vec3d _g;
  //solver
  scalarD _sigma0;
  scalarD _cThres;
  scalarD _alphaThres;
  scalarD _newtonThres;
  scalarD _initPenalty;
  scalarD _maxPenalty;
  scalarD _lineSearchThres;
  bool _useAugLag;
  bool _callback;
  bool _sparse;
  sizeType _maxIterNewton;
  sizeType _maxIterAugLag;
};
template <typename T>
class PhysicsRegistration : public GraspPlanner<T>
{
public:
  DECL_MAP_TYPES_T
  using GraspPlanner<T>::env;
  using GraspPlanner<T>::rad;
  using GraspPlanner<T>::solveDenseQP;
  using GraspPlanner<T>::solveSparseQP;
  using GraspPlanner<T>::_pnss;
  using GraspPlanner<T>::_body;
  using GraspPlanner<T>::_info;
  using GraspPlanner<T>::_objs;
  using GraspPlanner<T>::_A;
  using GraspPlanner<T>::_b;
  using GraspPlanner<T>::_l;
  using GraspPlanner<T>::_u;
  using GraspPlanner<T>::_gl;
  using GraspPlanner<T>::_gu;
  struct Penetration
  {
    bool operator<(const Penetration& other) const;
    bool operator!=(const Penetration& other) const;
    bool operator==(const Penetration& other) const;
    void writeVTK(const std::string& path) const;
    std::vector<std::tuple<sizeType,Vec3T,T>> _penetratedPoints;
    std::tuple<sizeType,Vec3T,T> _deepestPenetration;
    sizeType _oid,_oidOther;
  };
  PhysicsRegistration();
  void reset(const std::vector<ObjMesh>& objs,T rad,bool convex=true,T SDFRes=0,T SDFExtension=0,bool SDFRational=false);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  //index set helper
  void setIndexModifier(std::function<void(sizeType,const Vec&)> func);
  bool existContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid) const;
  bool existPointCloudMap(sizeType pid,sizeType oid) const;
  void addContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid,T mu,sizeType nDir);
  void addPointCloudMap(sizeType pid,sizeType oid,const Vec3T& objPos);
  void clearContactIndexMap();
  void clearPointCloudMap();
  void writeContactVTK(const Vec& x,const std::string& path) const;
  void getDeepestPenetration(std::set<Penetration>& pss) const;
  void getDeepestPenetration(Penetration& p) const;
  void debugPenetration(sizeType iter,T scale);
  sizeType nrObject() const;
  //optimize
  Vec optimize(bool debug,const Vec& init,const PointCloudObject<T>& env,const PointCloudObject<T>* obj,PhysicsRegistrationParameter& ops);
  bool assembleAugLag(Vec x,const std::vector<T>& lambda,T penalty,bool update,T& e,Vec* g=NULL,MatT* h=NULL);
  bool assembleAugLag(Vec x,const std::vector<T>& lambda,T penalty,bool update,T& e,Vec* g=NULL,SMat* h=NULL);
  bool lineSearchThresViolated(const Vec& x,const Vec& xNew,const PhysicsRegistrationParameter& ops) const;
  Vec optimizeNewton(Vec x,std::vector<T>& lambda,T penalty,PhysicsRegistrationParameter& ops);
  Vec optimizeAugLag(Vec x,PhysicsRegistrationParameter& ops);
  Vec optimizeSQP(Vec x,PhysicsRegistrationParameter& ops);
  T computePhi(T e,const Vec& c,T sigma) const;
  template <typename MAT>
  T computePhiPred(const Vec& d,T e,const MAT& h,const Vec& g,Vec c,const MAT& cjac,T sigma) const;
  template <typename MAT>
  T computeLinearizedCInf(const Vec& d,const Vec& c,const MAT& cjac) const;
  void debugSystemAugLag(const Vec& x);
private:
  std::function<void(sizeType,const Vec&)> _indexModifier;
  std::vector<std::vector<Node<sizeType,BBox<scalarD>>>> _bvhss;
  std::vector<scalarD> _pLMax;
};

PRJ_END

#endif
