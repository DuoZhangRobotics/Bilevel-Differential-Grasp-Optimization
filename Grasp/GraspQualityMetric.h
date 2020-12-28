#ifndef GRASP_QUALITY_METRIC_H
#define GRASP_QUALITY_METRIC_H

#include <CommonFile/geom/BVHNode.h>
#include <CommonFile/ObjMesh.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

#include <Utils/ArticulatedBodyPragma.h>
struct ALIGN_16 ObjMeshGeomCellExact;
template <typename T>
class GraspQualityMetric : public SerializableBase
{
public:
  DECL_MAP_TYPES_T
  GraspQualityMetric();
  void reset(ObjMesh& obj,T rad,sizeType dRes=4,const Mat6T& M=Mat6T::Identity(),T mu=0.7f);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  void writeVTK(const std::string& path,T len) const;
  const std::vector<Node<sizeType,BBox<scalar>>>& getBVH() const;
  T computeQInfMean(const Vec& w,Vec* g=NULL) const;
  T computeQInf(const Vec& w,Vec* g=NULL) const;
  T computeQ1(const Vec& w,Vec* g=NULL) const;
  const ObjMeshGeomCellExact& dist() const;
  const Mat3XT& pss() const;
  const Mat3XT& nss() const;
  const MatT& gij() const;
  void debug(sizeType iter);
protected:
  static T computeGij(const Vec3T& p,const Vec3T& n,const Vec6T& d,const Mat6T& M,T mu);
  std::vector<Node<sizeType,BBox<scalar>>> _bvh;
  std::shared_ptr<ObjMeshGeomCellExact> _distExact;
  Mat3XT _pss,_nss;
  ObjMesh _m;
  MatT _gij;
  T _rad;
};

PRJ_END

#endif
