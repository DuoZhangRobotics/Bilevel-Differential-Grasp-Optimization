#ifndef GRASP_PLANNER_H
#define GRASP_PLANNER_H

#include "GraspQualityMetric.h"
#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>

PRJ_BEGIN

enum METRIC_TYPE
{
  Q_INF_MEAN,
  Q_INF,
  Q_1,
  NO_METRIC,
};
struct ALIGN_16 ObjMeshGeomCellExact;
template <typename T>
class ArticulatedObjective;
template <typename T>
class GraspPlanner : public SerializableBase
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef std::vector<std::pair<Mat3XT,Mat3XT>> PNSS;
  GraspPlanner();
  void reset(const std::string& path,T rad);
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
  Vec optimize(const Vec& init,GraspQualityMetric<T>& obj,T d0=1,T alpha=1,METRIC_TYPE m=NO_METRIC,T coefM=1000,T coefOC=1,T coefCC=1,T coefO=1e-3f,T coefS=1e-3f) const;
  Vec optimizeNewton(Vec x,std::vector<std::shared_ptr<ArticulatedObjective<T>>>& objs,sizeType maxIter=1000,T gThres=1e-7f,T alphaThres=1e-6f,bool callback=true) const;
  T area() const;
  T rad() const;
protected:
  std::vector<std::shared_ptr<ObjMeshGeomCellExact>> _distExact;
  ArticulatedBody _body;
  //mimic
  MatT _A;
  Vec _b,_l,_u;
  //sampled points
  PNSS _pnss;
  T _rad;
};
template <typename T>
class ArticulatedObjective
{
public:
  DECL_MAP_FUNCS
  DECL_MAP_TYPES_T
  ArticulatedObjective(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj);
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g=NULL,MatT* h=NULL)=0;
  virtual int operator()(const Vec& x,T& e,Vec* g=NULL,MatT* h=NULL);
  void debug(const Vec& x);
protected:
  std::vector<KDOP18<scalar>> updateBVH(const PBDArticulatedGradientInfo<T>& info) const;
  const GraspPlanner<T>& _planner;
  const GraspQualityMetric<T>& _obj;
};

PRJ_END

#endif
