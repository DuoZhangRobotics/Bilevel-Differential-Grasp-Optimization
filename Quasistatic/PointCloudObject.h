#ifndef POINT_CLOUD_OBJECT_H
#define POINT_CLOUD_OBJECT_H

#include <CommonFile/geom/BVHNode.h>
#include <CommonFile/ObjMesh.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
struct ObjMeshGeomCellExact;
struct ArticulatedBody;
template <typename T>
class PointCloudObject : public SerializableBase
{
public:
  DECL_MAP_TYPES_T
  PointCloudObject();
  void reset(ObjMesh& obj,T rad);
  void resetPointCloud(ArticulatedBody& body,const Vec& x,const Mat4& m,const Mat4& prj,const Vec2i& res);
  void resetGraspable(ObjMesh& obj,T rad,sizeType dRes=4,const Mat6T& M=Mat6T::Identity(),T mu=0.1f,bool torque=false);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  const std::vector<Node<sizeType,BBox<scalarD>>>& getBVH() const;
  void writeVTK(const std::string& path,T len,T normalExtrude=0) const;
  T computeQInfBarrier(const Vec& w,T r,T d0,Vec* g=NULL) const;
  T computeQInf(const Vec& w,Vec* g=NULL) const;
  T computeQ1(const Vec& w,Vec* g=NULL) const;
  const ObjMeshGeomCellExact& dist() const;
  Mat3XT pss(T normalExtrude) const;
  const Mat3XT& pss() const;
  const Mat3XT& nss() const;
  const Coli& idss() const;
  const MatT& gij() const;
  void debug(sizeType iter);
protected:
  static T computeGij(const Vec3T& p,const Vec3T& n,const Vec6T& d,const Mat6T& M,T mu);
  void buildGij(sizeType dRes,const Mat6T& M,T mu,bool torque);
  void samplePoints(T rad);
  void buildBVH();
  //data
  std::vector<Node<sizeType,BBox<scalarD>>> _bvh;
  std::shared_ptr<ObjMeshGeomCellExact> _distExact;
  Mat3XT _pss,_nss;
  Coli _idss;
  ObjMesh _m;
  MatT _gij;
  T _rad;
};

PRJ_END

#endif
