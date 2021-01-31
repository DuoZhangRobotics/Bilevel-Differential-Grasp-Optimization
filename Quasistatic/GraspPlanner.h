#ifndef GRASP_PLANNER_H
#define GRASP_PLANNER_H

#include "MetricEnergy.h"
#include "PointCloudObject.h"
#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Optimizer/QCQPSolverMosek.h>
#include <Optimizer/QCQPSolverGurobi.h>
#include <Utils/ParallelVector.h>
#include <Utils/Options.h>

PRJ_BEGIN

struct ConvexHullExact;
struct ObjMeshGeomCellExact;
struct GraspPlannerParameter
{
  GraspPlannerParameter(Options& ops);
  void reset(Options& ops);
  static void initOptions(GraspPlannerParameter& sol);
  //data
  scalarD _d0;
  scalarD _alpha;
  sizeType _metric,_activation;
  scalarD _normalExtrude;
  scalarD _FGTThres;
  scalarD _coefM;
  scalarD _coefOC;
  scalarD _coefCC;
  scalarD _coefO;
  scalarD _coefS;
  bool _useGJK;
  //solver
  scalarD _rho0;
  scalarD _thres;
  scalarD _alphaThres;
  bool _callback;
  bool _sparse;
  sizeType _maxIter;
};
template <typename T>
struct PBDArticulatedGradientInfo;
template <typename T>
class ArticulatedObjective;
template <typename T>
class Environment;
template <typename T>
class GraspPlanner : public SerializableBase
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef std::vector<std::pair<Mat3XT,Mat3XT>> PNSS;
  GraspPlanner();
  void reset(T rad,bool convex=true,T SDFRes=0,T SDFExtension=0,bool SDFRational=false,bool checkValid=true);
  void reset(const std::string& path,T rad,bool convex=true,T SDFRes=0,T SDFExtension=0,bool SDFRational=false);
  void fliterSample(std::function<bool(sizeType lid,const Vec3T& p,const Vec3T& n)> f);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  ArticulatedBody& body();
  const ArticulatedBody& body() const;
  const Environment<T>& env(sizeType jid) const;
  const ObjMeshGeomCellExact& dist(sizeType jid) const;
  const std::vector<std::pair<Mat3XT,Mat3XT>>& pnss() const;
  void writeVTK(const Vec& x,const std::string& path,T len) const;
  void writeLocalVTK(const std::string& path,T len) const;
  void writeLimitsVTK(const std::string& path) const;
  //optimize
  Vec optimize(bool debug,const Vec& init,PointCloudObject<T>& object,GraspPlannerParameter& ops);
  bool solveDenseQP(Vec& d,const Vec& x,const Vec& g,MatT& h,const Vec* c,const MatT* cjac,T TR,T rho);
  bool solveSparseQP(Vec& d,const Vec& x,const Vec& g,SMat& h,const Vec* c,const SMat* cjac,T TR,T rho,T& reg);
  bool assemble(Vec x,bool update,T& e,Vec* g=NULL,MatT* h=NULL,Vec* c=NULL,MatT* cjac=NULL);
  bool assemble(Vec x,bool update,T& e,Vec* g=NULL,SMat* h=NULL,Vec* c=NULL,SMat* cjac=NULL);
  Vec optimizeSQP(Vec x,GraspPlannerParameter& ops);
  void debugSystem(const Vec& x);
  const SMat& A() const;
  const Vec& b() const;
  T area() const;
  T rad() const;
  bool validSample(sizeType l,const PBDArticulatedGradientInfo<T>& info,const Vec3T& p) const;
protected:
  std::vector<std::shared_ptr<Environment<T>>> _env;
  PBDArticulatedGradientInfo<T> _info;
  DSSQPObjectiveCompound<T> _objs;
  ArticulatedBody _body;
  //mimic
  SMat _A;
  QCQPSolverGurobi<T> _sol;
  Vec _b,_l,_u,_gl,_gu;
  //sampled points
  PNSS _pnss;
  T _rad;
};

PRJ_END

#endif
