#ifndef ENVIRONMENT_H
#define ENVIRONMENT_H

#include <CommonFile/geom/ObjMeshGeomCell.h>

PRJ_BEGIN

struct ArticulatedBody;
class ObjMeshGeomCellExact;
class ConvexHullExact;
template <typename T>
class Environment
{
public:
  typedef Eigen::Matrix<T,3,1> Vec3T;
  typedef Eigen::Matrix<T,2,2> Mat2T;
  typedef Eigen::Matrix<T,3,3> Mat3T;
  virtual ~Environment() {}
  virtual T phi(const Vec3T& x,Vec3T* g=NULL) const=0;
  virtual Vec3T phiGrad(const Vec3T& x,Mat3T* h=NULL) const=0;
  virtual BBox<scalarD> getBB() const=0;
  virtual ObjMesh getMesh() const=0;
  virtual bool empty() const;
  virtual ObjMesh getMeshProj2D() const;
  virtual void debug(sizeType nrIter=100);
};
template <typename T>
class EnvironmentExact : public Environment<T>, public SerializableBase
{
public:
  using typename Environment<T>::Vec3T;
  using typename Environment<T>::Mat2T;
  using typename Environment<T>::Mat3T;
  EnvironmentExact();
  EnvironmentExact(const ObjMeshGeomCellExact& obj);
  EnvironmentExact(const ConvexHullExact& obj);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual T phi(const Vec3T& x,Vec3T* g=NULL) const override;
  virtual Vec3T phiGrad(const Vec3T& x,Mat3T* h=NULL) const override;
  virtual BBox<scalarD> getBB() const override;
  virtual ObjMesh getMesh() const override;
  virtual bool empty() const override;
  virtual void debugVTK(const std::string& path,sizeType res=100) const;
  virtual ObjMeshGeomCell createStair(scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n);
  virtual ObjMeshGeomCell createHills(scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res);
  virtual ObjMeshGeomCell createZigZag(scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z=1);
  virtual ObjMeshGeomCell createFloor(scalarD x,scalarD y,scalarD z);
  virtual ObjMeshGeomCell createFloor(const Vec4d& plane);
  const ObjMeshGeomCellExact& getObj() const;
private:
  std::shared_ptr<ObjMeshGeomCellExact> _obj;
};
template <typename T>
class EnvironmentCubic : public Environment<T>, public SerializableBase
{
public:
  using typename Environment<T>::Vec3T;
  using typename Environment<T>::Mat2T;
  using typename Environment<T>::Mat3T;
  EnvironmentCubic();
  EnvironmentCubic(scalarD dx);
  EnvironmentCubic(const ObjMeshGeomCell& obj,scalarD dx,scalarD enlarge=0);
  EnvironmentCubic(const ObjMeshGeomCellExact& obj,scalarD dx,scalarD enlarge=0);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  virtual T phi(const Vec3T& x,Vec3T* g=NULL) const override;
  virtual Vec3T phiGrad(const Vec3T& x,Mat3T* h=NULL) const override;
  virtual BBox<scalarD> getBB() const override;
  virtual ObjMesh getMesh() const override;
  virtual bool empty() const override;
  virtual void debugVTK(const std::string& path) const;
  virtual void writeDistVTK(const std::string& path) const;
  virtual void createStair(scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n);
  virtual void createHills(scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res);
  virtual void createZigZag(scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z=1);
  virtual void createFloor(scalarD x,scalarD y,scalarD z);
  virtual void createFloor(const Vec4d& plane);
protected:
  void buildDist(const ObjMeshGeomCell& env);
  void buildDist(const EnvironmentExact<T>& env);
  static T interp1D(T t,std::function<T(sizeType)> f);
  static T interp1DDiff(T t,std::function<T(sizeType)> f);
  static T interp1DDDiff(T t,std::function<T(sizeType)> f);
  ScalarFieldD _dist;
  ObjMesh _mesh;
  scalarD _dx;
  scalarD _enlarge;
};
template <typename T>
class EnvironmentHeight : public EnvironmentCubic<T>
{
public:
  using typename Environment<T>::Vec3T;
  using typename Environment<T>::Mat2T;
  using typename Environment<T>::Mat3T;
  using EnvironmentCubic<T>::_dist;
  using EnvironmentCubic<T>::_mesh;
  using EnvironmentCubic<T>::_dx;
  using EnvironmentCubic<T>::interp1D;
  using EnvironmentCubic<T>::interp1DDiff;
  using EnvironmentCubic<T>::interp1DDDiff;
  EnvironmentHeight();
  EnvironmentHeight(scalarD dx);
  EnvironmentHeight(const std::string& path,bool is2D=false,scalar dxMul=1.0);
  EnvironmentHeight(const ObjMeshGeomCell& obj,scalarD dx);
  virtual T phi(const Vec3T& x,Vec3T* g=NULL) const override;
  virtual Vec3T phiGrad(const Vec3T& x,Mat3T* h=NULL) const override;
  virtual BBox<scalarD> getBB() const override;
  virtual ObjMesh getMesh() const override;
  virtual void createStair(scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n) override;
  virtual void createHills(scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res) override;
  virtual void createZigZag(scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z=1) override;
  virtual void createFloor(scalarD x,scalarD y,scalarD z) override;
  virtual void createFloor(const Vec4d& plane) override;
private:
  void buildDist(const ObjMeshGeomCell& cell);
};

PRJ_END

#endif
