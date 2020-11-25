#ifndef STATIC_GEOM_H
#define STATIC_GEOM_H

#include "../MathBasic.h"
#include "../ObjMesh.h"
#include "BVHNode.h"

PRJ_BEGIN

template <typename T,int dim> class OBBTpl;
template <typename T> class PlaneTpl;

struct ALIGN_16 StaticGeomCell : public Serializable {
  friend struct VertexCallback;
  friend class StaticGeom;
public:
  using Serializable::read;
  using Serializable::write;
  EIGEN_DEVICE_FUNC StaticGeomCell(const std::string& type);
  EIGEN_DEVICE_FUNC StaticGeomCell(const Mat4& T,sizeType dim,const std::string& type);
  EIGEN_DEVICE_FUNC virtual ~StaticGeomCell() {}
  template <typename T> void setUserData(T* data) {
    _userData=(void*)(data);
  }
  template <typename T> T* getUserData() {
    return (T*)(_userData);
  }
  void generateUV(scalar scale);  //used for rendering, but the collision detection ob ObjMeshGeomCell will be affected!
  void subdivideMesh(sizeType nrSubd);  //used for rendering, but the collision detection ob ObjMeshGeomCell will be affected!
  virtual void getMesh(ObjMesh& mesh,bool ref=false,bool render=false) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBB(bool ref=false) const;
  DEVICE_ONLY_FUNC virtual bool dist(const Vec3& pt,Vec3& n) const;
  DEVICE_ONLY_FUNC virtual bool closest(const Vec3& pt,Vec3& n,Vec3* normal=NULL) const;
  DEVICE_ONLY_FUNC virtual scalar rayQuery(Vec3 x0,Vec3 dir) const;
  virtual bool read(std::istream& is,IOData* dat);
  virtual bool write(std::ostream& os,IOData* dat) const;
  const std::vector<Node<sizeType>>& bvh() const;
  const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss() const;
  const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss() const;
  void debugDistQuery(bool closestTest,scalar nScale);
  void debugDistQuery(scalar nScale);
  const Mat4& getInvT() const;
  const Mat4& getT() const;
  void setT(const Mat4& T);
  sizeType getRes() const;
  virtual void setRes(sizeType res);
  sizeType dim() const;
  ALIGN_16 sizeType _index;
protected:
  virtual void getMeshInner(ObjMesh& mesh) const;
  DEVICE_ONLY_FUNC virtual BBox<scalar> getBBInner() const;
  DEVICE_ONLY_FUNC virtual bool distInner(const Vec3& pt,Vec3& n) const;
  DEVICE_ONLY_FUNC virtual bool closestInner(const Vec3& pt,Vec3& n,Vec3* normal=NULL) const;
  DEVICE_ONLY_FUNC virtual scalar rayQueryInner(const Vec3& x0,const Vec3& dir) const;
  void buildBVHBottomUp();
  virtual void build(bool buildBVHAsWell,bool bottomUp=false);
  virtual void generateUVInner(ObjMesh& mesh,scalar scale) const;
  //vertex
  ALIGN_16 std::vector<Vec2,Eigen::aligned_allocator<Vec2> > _tss;
  ALIGN_16 std::vector<Vec3,Eigen::aligned_allocator<Vec3> > _vss;
  ALIGN_16 std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _iss;
  ALIGN_16 std::vector<Node<sizeType> > _bvh;
  //data
  ALIGN_16 Mat4 _T,_invT;
  ALIGN_16 sizeType _dim;
  ALIGN_16 void* _userData;
  ALIGN_16 sizeType _res;
};
class StaticGeomCallback
{
public:
  virtual ~StaticGeomCallback() {}
  virtual void onCollideVertex(const Vec3& x,const Vec3& n,std::shared_ptr<StaticGeomCell> c) =0;
};
class StaticGeom : public Serializable
{
public:
  using Serializable::read;
  using Serializable::write;
  //method
  EIGEN_DEVICE_FUNC StaticGeom();
  EIGEN_DEVICE_FUNC StaticGeom(sizeType dim);
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar> > >& getBVH() const;
  std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar> > >& getBVH();
  const StaticGeomCell& getG(sizeType i) const;
  std::shared_ptr<StaticGeomCell> getGPtr(sizeType i) const;
  StaticGeomCell& getG(sizeType i);
  sizeType nrG() const;
  sizeType depth() const;
  void clear();
  void assemble();
  void parityCheck();
  bool update(scalar expand=1E-3f,bool dynamic=true);
  bool dist(const Vec3& pt,Vec3& n,std::shared_ptr<StaticGeomCell>& cell) const;
  void closest(const Vec3& pt,Vec3& n,Vec3* normal,std::shared_ptr<StaticGeomCell>& cell) const;
  bool rayQuery(const Vec3& pt0,Vec3& dir,std::shared_ptr<StaticGeomCell>& cell,Vec3& r) const;
  void collideVertex(const StaticGeom& other,StaticGeomCallback& cb) const;
  void forceBuild() const;
  virtual void setRes(sizeType res);
  //remove geometry
  void removeGeomCell(std::shared_ptr<StaticGeomCell> c);
  //add geometry
  void addGeomCell(std::shared_ptr<StaticGeomCell> c);
  void addGeomBox(const Mat4& trans,const BBox<scalar>& bb,scalar depth=0.0f);
  void addGeomBox(const Mat4& trans,const Vec3& ext,scalar depth=0.0f);
  void addGeomBox(const OBBTpl<scalar,2>& obb,scalar depth=0.0f);
  void addGeomBox(const OBBTpl<scalar,3>& obb,scalar depth=0.0f);
  void addGeomSphericalBox(const Mat4& trans,const Vec4& ext);
  void addGeomCylinder(const Mat4& trans,scalar rad,scalar y,bool capsule=false);
  void addGeomPlane(const Mat4& trans,const Vec4& plane,scalar ext=100.0f);
  void addGeomPlane(const Mat4& trans,const PlaneTpl<scalar>& plane,scalar ext=100.0f);
  void addGeomSphere(const Vec3& ctr,scalar rad,scalar depth=0.0f);
  void addGeomMesh(const Mat4& trans,const ObjMesh& mesh,scalar depth=0.0f);
  void addGeomMesh(const Mat4& trans,const std::string& path,scalar depth=0.0f);
  void addGeomHeightField(sizeType dimH,scalar h0,scalar hr,scalar sz,scalar cellSz);
  void addGeomSolidBox(const Mat4& trans,const Vec3& ext,scalar thick);
  void addGeomSolidSphere(const Mat4& trans,scalar rad,scalar thick);
  //IO
  void writeVTK(const std::string& path) const;
  void writeBVH() const;
  bool read(std::istream& is,IOData* dat);
  bool write(std::ostream& os,IOData* dat) const;
  std::shared_ptr<SerializableBase> copy() const;
  //tools
  static bool write(const std::shared_ptr<StaticGeomCell>& cell,std::ostream& os);
  static bool read(std::shared_ptr<StaticGeomCell>& cell,std::istream& is);
  static void registerType(IOData* dat);
  static void debugRayQuery(const std::string& path,std::shared_ptr<StaticGeomCell> cell,sizeType nr=100);
  static void debugVertexQuery(const std::string& path,std::shared_ptr<StaticGeomCell> cell,sizeType nr=100);
  static void debugVertexDistQuery(const std::string& path,std::shared_ptr<StaticGeomCell> cell,sizeType nr=100);
  static void debugRayVertexQuery(bool box,bool sphere,bool cylinder,bool sphericalBox,bool capsule,bool objMesh,bool height,bool twoSphere,bool threeSphere);
private:
  //mesh data
  std::vector<std::shared_ptr<StaticGeomCell> > _css;
  std::shared_ptr<std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar> > > > _bvh;
  sizeType _dim;
};

PRJ_END

#endif
