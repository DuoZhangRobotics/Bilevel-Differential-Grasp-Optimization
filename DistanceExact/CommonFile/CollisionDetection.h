#ifndef COLLISION_DETECTION_H
#define COLLISION_DETECTION_H

#include "MathBasic.h"
#include "IOFwd.h"

PRJ_BEGIN

template <typename T,int DIM=3>
struct ALIGN_16 BBox {
  static const int dim=DIM;
  typedef Eigen::Matrix<T,DIM,1> PT;
  typedef Eigen::Matrix<T,2,1> PT2;
  EIGEN_DEVICE_FUNC BBox();
  EIGEN_DEVICE_FUNC BBox(const PT& p);
  EIGEN_DEVICE_FUNC BBox(const PT& minC,const PT& maxC);
  template <typename T2>
  EIGEN_DEVICE_FUNC BBox(const BBox<T2, DIM>& other) {
    copy(other);
  }
  EIGEN_DEVICE_FUNC virtual ~BBox();
  template <typename T2>
  EIGEN_DEVICE_FUNC BBox& operator=(const BBox<T2,DIM>& other);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC static BBox createMM(const PT& minC,const PT& maxC);
  EIGEN_DEVICE_FUNC static BBox createME(const PT& minC,const PT& extent);
  EIGEN_DEVICE_FUNC static BBox createCE(const PT& center,const PT& extent);
  EIGEN_DEVICE_FUNC BBox getIntersect(const BBox& other) const;
  EIGEN_DEVICE_FUNC BBox getUnion(const BBox& other) const;
  EIGEN_DEVICE_FUNC BBox getUnion(const PT& point) const;
  EIGEN_DEVICE_FUNC BBox getUnion(const PT& ctr,const T& rad) const;
  EIGEN_DEVICE_FUNC void setIntersect(const BBox& other);
  EIGEN_DEVICE_FUNC void setUnion(const BBox& other);
  EIGEN_DEVICE_FUNC void setUnion(const PT& point);
  EIGEN_DEVICE_FUNC void setUnion(const PT& ctr,const T& rad);
  EIGEN_DEVICE_FUNC void setPoints(const PT& a,const PT& b,const PT& c);
  EIGEN_DEVICE_FUNC PT minCorner() const;
  EIGEN_DEVICE_FUNC PT maxCorner() const;
  EIGEN_DEVICE_FUNC void enlargedEps(T eps);
  EIGEN_DEVICE_FUNC BBox enlargeEps(T eps) const;
  EIGEN_DEVICE_FUNC void enlarged(T len,const sizeType d=DIM);
  EIGEN_DEVICE_FUNC BBox enlarge(T len,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC PT lerp(const PT& frac) const;
  EIGEN_DEVICE_FUNC bool empty() const;
  template <int DIM2>
  EIGEN_DEVICE_FUNC bool containDim(const PT& point) const;
  EIGEN_DEVICE_FUNC bool contain(const BBox& other,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC bool contain(const PT& point,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC bool contain(const PT& point,const T& rad,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC void reset();
  EIGEN_DEVICE_FUNC PT getExtent() const;
  EIGEN_DEVICE_FUNC T distTo(const BBox& other,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC T distTo(const PT& pt,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC T distToSqr(const PT& pt,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC PT closestTo(const PT& pt,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC bool intersect(const PT& p,const PT& q,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC bool intersect(const PT& p,const PT& q,T& s,T& t,const sizeType d=DIM) const;
  EIGEN_DEVICE_FUNC bool intersect(const BBox& other,const sizeType& d=DIM) const;
  EIGEN_DEVICE_FUNC PT2 project(const PT& a,const sizeType d=DIM) const;
  template<typename T2>
  EIGEN_DEVICE_FUNC BBox& copy(const BBox<T2,DIM>& other);
  EIGEN_DEVICE_FUNC T perimeter(const sizeType d=DIM) const;
  ALIGN_16 PT _minC;
  ALIGN_16 PT _maxC;
};
template <typename T>
class ALIGN_16 LineSegTpl
{
public:
  typedef typename Eigen::Matrix<T,3,1> PT;
  typedef typename Eigen::Matrix<T,2,1> PT2;
  typedef typename Eigen::Matrix<T,2,2> MAT2;
  EIGEN_DEVICE_FUNC LineSegTpl();
  EIGEN_DEVICE_FUNC LineSegTpl(const PT& x,const PT& y);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC T length() const;
  EIGEN_DEVICE_FUNC PT circumcenter() const;
  EIGEN_DEVICE_FUNC PT masscenter() const;
  EIGEN_DEVICE_FUNC PT normal() const;
  EIGEN_DEVICE_FUNC PT gradDir(sizeType id) const;
  EIGEN_DEVICE_FUNC T signedArea() const;
  EIGEN_DEVICE_FUNC bool intersect(const LineSegTpl<T>& l,T& t,bool infinite=false) const;
  EIGEN_DEVICE_FUNC void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT& b) const;
  EIGEN_DEVICE_FUNC void calcLineDist(const LineSegTpl<T>& l,T& sqrDistance,T& a,T& b) const;
  EIGEN_DEVICE_FUNC void writeVTK(VTKWriter<T>& os) const;
public:
  EIGEN_DEVICE_FUNC T getClampedRoot(T slope,T h0,T h1) const;
  EIGEN_DEVICE_FUNC void computeIntersection(T mF00,T mF10,T mB,const T sValue[2],const sizeType classify[2],sizeType edge[2],T end[2][2]) const;
  EIGEN_DEVICE_FUNC void computeMinimumParameters(T mB,T mC,T mE,T mG00,T mG01,T mG10,T mG11,const sizeType edge[2],const T end[2][2],T parameter[2]) const;
  //data
  ALIGN_16 PT _x;
  ALIGN_16 PT _y;
};
template <typename T>
class ALIGN_16 PlaneTpl
{
public:
  typedef typename Eigen::Matrix<T,3,1> PT;
  typedef typename Eigen::Matrix<T,2,1> PT2;
  EIGEN_DEVICE_FUNC PlaneTpl();
  EIGEN_DEVICE_FUNC PlaneTpl(const PT& x0,const PT& n);
  EIGEN_DEVICE_FUNC PlaneTpl(const PT& a,const PT& b,const PT& c);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC T side(const PT& p) const;
  EIGEN_DEVICE_FUNC bool intersect(const BBox<T>& bb) const;
  EIGEN_DEVICE_FUNC PT2 project(const PT& d) const;
  EIGEN_DEVICE_FUNC void writeVTK(VTKWriter<T>& os) const;
public:
  //data
  ALIGN_16 PT _x0;
  ALIGN_16 PT _n;
};
template <typename T>
class ALIGN_16 TriangleTpl
{
public:
  typedef typename Eigen::Matrix<T,3,1> PT;
  typedef typename Eigen::Matrix<T,2,1> PT2;
  typedef typename Eigen::Matrix<T,2,2> MAT2;
  typedef typename Eigen::Matrix<T,3,3> MAT3;
  EIGEN_DEVICE_FUNC TriangleTpl();
  EIGEN_DEVICE_FUNC TriangleTpl(const PT& a,const PT& b,const PT& c);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC PT circumcenter() const;
  EIGEN_DEVICE_FUNC PT masscenter() const;
  EIGEN_DEVICE_FUNC PT bary(const PT& pt) const;
  EIGEN_DEVICE_FUNC T area() const;
  EIGEN_DEVICE_FUNC T signedVolume() const;
  EIGEN_DEVICE_FUNC PT normal() const;
  EIGEN_DEVICE_FUNC PT gradDir(sizeType id) const;
  EIGEN_DEVICE_FUNC PT height(const PT& a,const PT& b,const PT& c) const;
  //intersect
  EIGEN_DEVICE_FUNC bool isInside(const PT& pt) const;
  EIGEN_DEVICE_FUNC bool isPlanePointInside(const PT& pt) const;
  EIGEN_DEVICE_FUNC bool updateIntr(T& s,T& t,T num,T denom) const;
  EIGEN_DEVICE_FUNC bool intersect(const LineSegTpl<T>& l,T& s,T& t) const;
  EIGEN_DEVICE_FUNC bool intersect(const LineSegTpl<T>& l,T& t,bool infinite=false,PT* abtOut=NULL) const;
  EIGEN_DEVICE_FUNC bool intersect(const TriangleTpl<T>& t) const;
  EIGEN_DEVICE_FUNC bool intersect(const BBox<T>& bb) const;
  EIGEN_DEVICE_FUNC bool calcLineDist(const LineSegTpl<T>& l,T& sqrDistance,PT& bt,PT2& bl) const;
  template <typename PT_BARY>
  EIGEN_DEVICE_FUNC void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT_BARY& b) const;
  EIGEN_DEVICE_FUNC bool calcTriangleDist(const TriangleTpl<T>& t2,T& sqrDistance,PT& bt,PT& bt2) const;
  EIGEN_DEVICE_FUNC PT2 project(const PT& d) const;
  EIGEN_DEVICE_FUNC void writeVTK(VTKWriter<T>& os) const;
  EIGEN_DEVICE_FUNC const PT& operator[](sizeType d) const;
public:
  //data
  ALIGN_16 PT _a;
  ALIGN_16 PT _b;
  ALIGN_16 PT _c;
};
template <typename T>
class ALIGN_16 TetrahedronTpl
{
public:
  typedef typename Eigen::Matrix<T,2,1> PT2;
  typedef typename Eigen::Matrix<T,3,1> PT;
  typedef typename Eigen::Matrix<T,4,1> PT4;
  typedef typename Eigen::Matrix<T,6,1> PT6;
  typedef typename Eigen::Matrix<T,3,3> MAT3;
  EIGEN_DEVICE_FUNC TetrahedronTpl();
  EIGEN_DEVICE_FUNC TetrahedronTpl(const PT& a,const PT& b,const PT& c,const PT& d);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC PT circumcenter() const;
  EIGEN_DEVICE_FUNC PT masscenter() const;
  EIGEN_DEVICE_FUNC PT4 bary(const PT& pt) const;
  EIGEN_DEVICE_FUNC bool isInside(const PT& pt) const;
  EIGEN_DEVICE_FUNC T volume() const;
  EIGEN_DEVICE_FUNC bool dualCellVolume(PT4& dv);
  EIGEN_DEVICE_FUNC static bool dualFaceVolume(T& va,T& vb,T& vc,const PT& a,const PT& b,const PT& c,const PT& cc);
  EIGEN_DEVICE_FUNC static T dihedralAngle(const PT& a,const PT& b,const PT& c,const PT& d);
  EIGEN_DEVICE_FUNC PT6 dihedralAngleTet();
  //intersection
  EIGEN_DEVICE_FUNC bool calcLineDist(const LineSegTpl<T>& l,PT4& bt,PT2& bl) const;
  EIGEN_DEVICE_FUNC void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT4& bc) const;
  EIGEN_DEVICE_FUNC void getPlane(const sizeType& i,PlaneTpl<T>& p,Vec3i& ind) const;
  EIGEN_DEVICE_FUNC const PT& getNode(const sizeType& i) const;
  EIGEN_DEVICE_FUNC void writeVTK(VTKWriter<T>& os) const;
public:
  //data
  ALIGN_16 PT _a;
  ALIGN_16 PT _b;
  ALIGN_16 PT _c;
  ALIGN_16 PT _d;
  ALIGN_16 bool _swap;
};
template <typename T,int DIM=3>
class OBBTpl;
template <typename T>
class ALIGN_16 OBBTpl<T,2>
{
public:
  typedef typename Eigen::Matrix<T,2,1> PT;
  typedef typename Eigen::Matrix<T,3,1> PT3;
  typedef typename Eigen::Matrix<T,2,2> MAT;
  typedef typename Eigen::Matrix<T,3,3> MAT3;
  EIGEN_DEVICE_FUNC OBBTpl();
  EIGEN_DEVICE_FUNC OBBTpl(const BBox<T,2>& bb);
  EIGEN_DEVICE_FUNC OBBTpl(const BBox<T,3>& bb);
  EIGEN_DEVICE_FUNC OBBTpl(const MAT& rot,const PT& trans,const BBox<T,2>& bb);
  EIGEN_DEVICE_FUNC OBBTpl(const MAT3& rot,const PT3& trans,const BBox<T,3>& bb);
  EIGEN_DEVICE_FUNC OBBTpl(const MAT& rot,const PT& trans,const PT& ext);
  EIGEN_DEVICE_FUNC OBBTpl(const MAT3& rot,const PT3& trans,const PT3& ext);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC bool closest(PT pt,PT& n,PT* normal) const;
  EIGEN_DEVICE_FUNC bool closestInner(const PT& pt,PT& n,PT* normal) const;
  EIGEN_DEVICE_FUNC bool intersect(const OBBTpl<T,2>& other) const;
  EIGEN_DEVICE_FUNC bool intersect(const BBox<T,2>& other) const;
  EIGEN_DEVICE_FUNC void writeVTK(VTKWriter<T>& os) const;
  ALIGN_16 MAT _rot;
  ALIGN_16 PT _trans;
  ALIGN_16 PT _ext;
};
template <typename T>
class ALIGN_16 OBBTpl<T,3>
{
public:
  typedef typename Eigen::Matrix<T,3,1> PT;
  typedef typename Eigen::Matrix<T,2,1> PT2;
  typedef typename Eigen::Matrix<T,3,3> MAT;
  EIGEN_DEVICE_FUNC OBBTpl();
  EIGEN_DEVICE_FUNC OBBTpl(const BBox<T,3>& bb);
  EIGEN_DEVICE_FUNC OBBTpl(const MAT& rot,const PT& trans,const BBox<T,3>& bb);
  EIGEN_DEVICE_FUNC OBBTpl(const MAT& rot,const PT& trans,const PT& ext);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC bool closest(PT pt,PT& n,PT* normal) const;
  EIGEN_DEVICE_FUNC bool closestInner(const PT& pt,PT& n,PT* normal) const;
  EIGEN_DEVICE_FUNC bool intersect(const OBBTpl<T,3>& other) const;
  EIGEN_DEVICE_FUNC bool intersect(const BBox<T,3>& other) const;
  EIGEN_DEVICE_FUNC void writeVTK(VTKWriter<T>& os) const;
  ALIGN_16 MAT _rot;
  ALIGN_16 PT _trans;
  ALIGN_16 PT _ext;
};
template <typename T>
class ALIGN_16 KDOP18
{
public:
  static const int dim=3;
  typedef typename Eigen::Matrix<T,3,1> PT;
  EIGEN_DEVICE_FUNC KDOP18();
  EIGEN_DEVICE_FUNC KDOP18(const PT& v);
  EIGEN_DEVICE_FUNC KDOP18(const PT& a,const PT& b);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC void reset();
  EIGEN_DEVICE_FUNC void empty();
  EIGEN_DEVICE_FUNC void enlarged(T len);
  EIGEN_DEVICE_FUNC void enlarged(T len,sizeType dim);
  EIGEN_DEVICE_FUNC void setPoints(const PT& a,const PT& b,const PT& c);
  EIGEN_DEVICE_FUNC void setUnion(const PT& p);
  EIGEN_DEVICE_FUNC void setUnion(const KDOP18& b);
  EIGEN_DEVICE_FUNC KDOP18<T> getUnion(const KDOP18& b) const;
  EIGEN_DEVICE_FUNC PT minCorner() const;
  EIGEN_DEVICE_FUNC PT maxCorner() const;
  EIGEN_DEVICE_FUNC bool intersect(const KDOP18& b,sizeType DIM=3) const;
  EIGEN_DEVICE_FUNC bool intersect(const BBox<T>& b,sizeType DIM=3) const;
  EIGEN_DEVICE_FUNC bool contain(const PT& p) const;
protected:
  EIGEN_DEVICE_FUNC static void getDistances(const PT& p,T &d3, T &d4, T &d5, T &d6, T &d7, T &d8);
  EIGEN_DEVICE_FUNC static void getDistances(const PT& p, T d[]);
  EIGEN_DEVICE_FUNC static T getDistances(const PT &p, int i);
  ALIGN_16 T _dist[18];
};
template <typename T>
class ALIGN_16 Sphere
{
public:
  typedef typename Eigen::Matrix<T,3,1> PT;
  EIGEN_DEVICE_FUNC Sphere();
  EIGEN_DEVICE_FUNC Sphere(const PT& ctr,const T& rad);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  EIGEN_DEVICE_FUNC bool closest(const PT& pt,PT& n,PT* normal) const;
  EIGEN_DEVICE_FUNC bool intersect(const PT& a,const PT& b) const;
  //data
  ALIGN_16 PT _ctr;
  ALIGN_16 T _rad;
};

typedef BBox<scalar,2> BOX2D;
typedef BBox<scalar,3> BOX3D;
typedef LineSegTpl<scalar> LineSeg;
typedef PlaneTpl<scalar> Plane;
typedef TriangleTpl<scalar> Triangle;
typedef TetrahedronTpl<scalar> Tetrahedron;
typedef OBBTpl<scalar,2> OBB2D;
typedef OBBTpl<scalar,3> OBB3D;

PRJ_END

#endif
