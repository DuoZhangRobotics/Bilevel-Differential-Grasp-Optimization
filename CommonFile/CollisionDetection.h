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
  BBox();
  BBox(const PT& p);
  BBox(const PT& minC,const PT& maxC);
  template <typename T2>
  BBox(const BBox<T2, DIM>& other) {
    copy(other);
  }
  virtual ~BBox();
  template <typename T2>
  BBox& operator=(const BBox<T2,DIM>& other);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  static BBox createMM(const PT& minC,const PT& maxC);
  static BBox createME(const PT& minC,const PT& extent);
  static BBox createCE(const PT& center,const PT& extent);
  BBox getIntersect(const BBox& other) const;
  BBox getUnion(const BBox& other) const;
  BBox getUnion(const PT& point) const;
  BBox getUnion(const PT& ctr,const T& rad) const;
  void setIntersect(const BBox& other);
  void setUnion(const BBox& other);
  void setUnion(const PT& point);
  void setUnion(const PT& ctr,const T& rad);
  void setPoints(const PT& a,const PT& b,const PT& c);
  PT minCorner() const;
  PT maxCorner() const;
  void enlargedEps(T eps);
  BBox enlargeEps(T eps) const;
  void enlarged(T len,const sizeType d=DIM);
  BBox enlarge(T len,const sizeType d=DIM) const;
  PT lerp(const PT& frac) const;
  bool empty() const;
  template <int DIM2>
  bool containDim(const PT& point) const;
  bool contain(const BBox& other,const sizeType d=DIM) const;
  bool contain(const PT& point,const sizeType d=DIM) const;
  bool contain(const PT& point,const T& rad,const sizeType d=DIM) const;
  void reset();
  PT getExtent() const;
  T distTo(const BBox& other,const sizeType d=DIM) const;
  T distTo(const PT& pt,const sizeType d=DIM) const;
  T distToSqr(const PT& pt,const sizeType d=DIM) const;
  PT closestTo(const PT& pt,const sizeType d=DIM) const;
  bool intersect(const PT& p,const PT& q,const sizeType d=DIM) const;
  bool intersect(const PT& p,const PT& q,T& s,T& t,const sizeType d=DIM) const;
  bool intersect(const BBox& other,const sizeType& d=DIM) const;
  PT2 project(const PT& a,const sizeType d=DIM) const;
  template<typename T2>
  BBox& copy(const BBox<T2,DIM>& other);
  T perimeter(const sizeType d=DIM) const;
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
  LineSegTpl();
  LineSegTpl(const PT& x,const PT& y);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  T length() const;
  PT circumcenter() const;
  PT masscenter() const;
  PT normal() const;
  PT gradDir(sizeType id) const;
  T signedArea() const;
  bool intersect(const LineSegTpl<T>& l,T& t,bool infinite=false) const;
  void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT& b) const;
  void calcLineDist(const LineSegTpl<T>& l,T& sqrDistance,T& a,T& b) const;
  void writeVTK(VTKWriter<T>& os) const;
public:
  T getClampedRoot(T slope,T h0,T h1) const;
  void computeIntersection(T mF00,T mF10,T mB,const T sValue[2],const sizeType classify[2],sizeType edge[2],T end[2][2]) const;
  void computeMinimumParameters(T mB,T mC,T mE,T mG00,T mG01,T mG10,T mG11,const sizeType edge[2],const T end[2][2],T parameter[2]) const;
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
  PlaneTpl();
  PlaneTpl(const PT& x0,const PT& n);
  PlaneTpl(const PT& a,const PT& b,const PT& c);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  T side(const PT& p) const;
  bool intersect(const BBox<T>& bb) const;
  PT2 project(const PT& d) const;
  void writeVTK(VTKWriter<T>& os) const;
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
  TriangleTpl();
  TriangleTpl(const PT& a,const PT& b,const PT& c);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  PT circumcenter() const;
  PT masscenter() const;
  PT bary(const PT& pt) const;
  T area() const;
  T signedVolume() const;
  PT normal() const;
  PT gradDir(sizeType id) const;
  PT height(const PT& a,const PT& b,const PT& c) const;
  //intersect
  bool isInside(const PT& pt) const;
  bool isPlanePointInside(const PT& pt) const;
  bool updateIntr(T& s,T& t,T num,T denom) const;
  bool intersect(const LineSegTpl<T>& l,T& s,T& t) const;
  bool intersect(const LineSegTpl<T>& l,T& t,bool infinite=false,PT* abtOut=NULL) const;
  bool intersect(const TriangleTpl<T>& t) const;
  bool intersect(const BBox<T>& bb) const;
  bool calcLineDist(const LineSegTpl<T>& l,T& sqrDistance,PT& bt,PT2& bl) const;
  template <typename PT_BARY>
  void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT_BARY& b) const;
  bool calcTriangleDist(const TriangleTpl<T>& t2,T& sqrDistance,PT& bt,PT& bt2) const;
  PT2 project(const PT& d) const;
  void writeVTK(VTKWriter<T>& os) const;
  const PT& operator[](sizeType d) const;
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
  TetrahedronTpl();
  TetrahedronTpl(const PT& a,const PT& b,const PT& c,const PT& d);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  PT circumcenter() const;
  PT masscenter() const;
  PT4 bary(const PT& pt) const;
  bool isInside(const PT& pt) const;
  T volume() const;
  bool dualCellVolume(PT4& dv);
  static bool dualFaceVolume(T& va,T& vb,T& vc,const PT& a,const PT& b,const PT& c,const PT& cc);
  static T dihedralAngle(const PT& a,const PT& b,const PT& c,const PT& d);
  PT6 dihedralAngleTet();
  //intersection
  bool calcLineDist(const LineSegTpl<T>& l,PT4& bt,PT2& bl) const;
  void calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT4& bc) const;
  void getPlane(const sizeType& i,PlaneTpl<T>& p,Vec3i& ind) const;
  const PT& getNode(const sizeType& i) const;
  void writeVTK(VTKWriter<T>& os) const;
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
  OBBTpl();
  OBBTpl(const BBox<T,2>& bb);
  OBBTpl(const BBox<T,3>& bb);
  OBBTpl(const MAT& rot,const PT& trans,const BBox<T,2>& bb);
  OBBTpl(const MAT3& rot,const PT3& trans,const BBox<T,3>& bb);
  OBBTpl(const MAT& rot,const PT& trans,const PT& ext);
  OBBTpl(const MAT3& rot,const PT3& trans,const PT3& ext);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  bool closest(PT pt,PT& n,PT* normal) const;
  bool closestInner(const PT& pt,PT& n,PT* normal) const;
  bool intersect(const OBBTpl<T,2>& other) const;
  bool intersect(const BBox<T,2>& other) const;
  void writeVTK(VTKWriter<T>& os) const;
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
  OBBTpl();
  OBBTpl(const BBox<T,3>& bb);
  OBBTpl(const MAT& rot,const PT& trans,const BBox<T,3>& bb);
  OBBTpl(const MAT& rot,const PT& trans,const PT& ext);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  bool closest(PT pt,PT& n,PT* normal) const;
  bool closestInner(const PT& pt,PT& n,PT* normal) const;
  bool intersect(const OBBTpl<T,3>& other) const;
  bool intersect(const BBox<T,3>& other) const;
  void writeVTK(VTKWriter<T>& os) const;
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
  KDOP18();
  KDOP18(const PT& v);
  KDOP18(const PT& a,const PT& b);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  void reset();
  void empty();
  void enlarged(T len);
  void enlarged(T len,sizeType dim);
  void setPoints(const PT& a,const PT& b,const PT& c);
  void setUnion(const PT& p);
  void setUnion(const KDOP18& b);
  KDOP18<T> getUnion(const KDOP18& b) const;
  PT minCorner() const;
  PT maxCorner() const;
  bool intersect(const KDOP18& b,sizeType DIM=3) const;
  bool intersect(const BBox<T>& b,sizeType DIM=3) const;
  bool contain(const PT& p) const;
protected:
  static void getDistances(const PT& p,T &d3, T &d4, T &d5, T &d6, T &d7, T &d8);
  static void getDistances(const PT& p, T d[]);
  static T getDistances(const PT &p, int i);
  ALIGN_16 T _dist[18];
};
template <typename T>
class ALIGN_16 Sphere
{
public:
  typedef typename Eigen::Matrix<T,3,1> PT;
  Sphere();
  Sphere(const PT& ctr,const T& rad);
  bool write(std::ostream& os) const;
  bool read(std::istream& is);
  bool closest(const PT& pt,PT& n,PT* normal) const;
  bool intersect(const PT& a,const PT& b) const;
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
