#ifndef GRASP_PLANNER_H
#define GRASP_PLANNER_H

#include "MetricEnergy.h"
#include "GraspQualityMetric.h"
#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Utils/ParallelVector.h>
#include <Utils/Options.h>

PRJ_BEGIN

struct ALIGN_16 ConvexHullExact;
struct ALIGN_16 ObjMeshGeomCellExact;
template <typename T>
class ArticulatedObjective;
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
  sizeType _maxIter;
};
template <typename T>
class GraspPlanner : public SerializableBase
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef std::vector<std::pair<Mat3XT,Mat3XT>> PNSS;
  GraspPlanner();
  void reset(const std::string& path,T rad,bool convex=true);
  void fliterSample(std::function<bool(sizeType lid,const Vec3T& p,const Vec3T& n)> f);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  ArticulatedBody& body();
  const ArticulatedBody& body() const;
  const ObjMeshGeomCellExact& dist(sizeType jid) const;
  const std::vector<std::pair<Mat3XT,Mat3XT>>& pnss() const;
  void writeVTK(const Vec& x,const std::string& path,T len) const;
  void writeLocalVTK(const std::string& path,T len) const;
  void writeLimitsVTK(const std::string& path) const;
  Vec optimize(bool debug,const Vec& init,GraspQualityMetric<T>& obj,GraspPlannerParameter& ops);
  Vec optimizeNewton(Vec x,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs,GraspPlannerParameter& ops) const;
  bool assemble(Vec x,PBDArticulatedGradientInfo<T>& info,bool update,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs,T& e,Vec* g=NULL,MatT* h=NULL,Vec* c=NULL,MatT* cjac=NULL) const;
  void debugSystem(const Vec& x,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs) const;
  sizeType nrAdditionalDOF(std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs) const;
  sizeType nrCons(std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs) const;
  const MatT& A() const;
  const Vec& b() const;
  T area() const;
  T rad() const;
protected:
  bool validSample(sizeType l,const PBDArticulatedGradientInfo<T>& info,const Vec3T& p) const;
  std::vector<std::shared_ptr<ObjMeshGeomCellExact>> _distExact;
  ArticulatedBody _body;
  //mimic
  MatT _A;
  Vec _b,_l,_u;
  //sampled points
  PNSS _pnss;
  T _rad;
};

PRJ_END

#endif
