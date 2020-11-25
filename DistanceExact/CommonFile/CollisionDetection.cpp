#include "CollisionDetection.h"
#include "IO.h"
using std::max;
using std::min;

PRJ_BEGIN

template <typename MAT,bool isInt>
struct InvTraits
{
  EIGEN_DEVICE_FUNC static MAT inverse(const MAT& m) {
    return m.inverse();
  }
};
template <typename MAT>
struct InvTraits<MAT,true>
{
  EIGEN_DEVICE_FUNC static MAT inverse(const MAT& m) {
    FUNCTION_NOT_IMPLEMENTED
    return m;
  }
};
template <typename T>
class SAT
{
public:
  typedef typename Eigen::Matrix<T,3,1> PT;
  typedef typename Eigen::Matrix<T,2,1> PT2;
public:
  static bool testSegment(const PT2& a,const PT2& b) {
    return a.x() > b.y() || b.x() > a.y();
  }
  template <typename T1,typename T2>
  static bool sep(const PT& d,const T1& A,const T2& B) {
    return testSegment(A.project(d),B.project(d));
  }
};
//BBox
template <typename T,int DIM>
EIGEN_DEVICE_FUNC BBox<T,DIM>::BBox() {
  reset();
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC BBox<T,DIM>::BBox(const PT& p):_minC(p),_maxC(p) {}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC BBox<T,DIM>::BBox(const PT& minC,const PT& maxC):_minC(minC),_maxC(maxC) {}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC BBox<T,DIM>::~BBox() {}
template <typename T,int DIM>
template <typename T2>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox& BBox<T,DIM>::operator=(const BBox<T2,DIM>& other) {
  copy(other);
  return *this;
}
template <typename T,int DIM>
bool BBox<T,DIM>::write(std::ostream& os) const
{
  writeBinaryData(_minC,os);
  writeBinaryData(_maxC,os);
  return os.good();
}
template <typename T,int DIM>
bool BBox<T,DIM>::read(std::istream& is)
{
  readBinaryData(_minC,is);
  readBinaryData(_maxC,is);
  return is.good();
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::createMM(const PT& minC,const PT& maxC) {
  return BBox(minC,maxC);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::createME(const PT& minC,const PT& extent) {
  return BBox(minC,minC+extent);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::createCE(const PT& center,const PT& extent) {
  return BBox(center-extent,center+extent);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::getIntersect(const BBox& other) const {
  return createMM(compMax(_minC,other._minC),compMin(_maxC,other._maxC));
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::getUnion(const BBox& other) const {
  return createMM(compMin(_minC,other._minC),compMax(_maxC,other._maxC));
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::getUnion(const PT& point) const {
  return createMM(compMin(_minC,point),compMax(_maxC,point));
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::getUnion(const PT& ctr,const T& rad) const {
  return createMM(compMin(_minC,ctr-PT::Constant(rad)),compMax(_maxC,ctr+PT::Constant(rad)));
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::setIntersect(const BBox& other) {
  _minC=compMax(_minC,other._minC);
  _maxC=compMin(_maxC,other._maxC);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::setUnion(const BBox& other) {
  _minC=compMin<PT>(_minC,other._minC);
  _maxC=compMax<PT>(_maxC,other._maxC);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::setUnion(const PT& point) {
  _minC=compMin(_minC,point);
  _maxC=compMax(_maxC,point);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::setUnion(const PT& ctr,const T& rad) {
  _minC=compMin(_minC,ctr-PT::Constant(rad));
  _maxC=compMax(_maxC,ctr+PT::Constant(rad));
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::setPoints(const PT& a,const PT& b,const PT& c) {
  _minC=compMin(compMin(a,b),c);
  _maxC=compMax(compMax(a,b),c);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::PT BBox<T,DIM>::minCorner() const {
  return _minC;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::PT BBox<T,DIM>::maxCorner() const {
  return _maxC;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::enlargedEps(T eps) {
  PT d=(_maxC-_minC)*(typename EigenTraits<PT>::ScalarType)(eps*0.5f);
  _minC-=d;
  _maxC+=d;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::enlargeEps(T eps) const {
  PT d=(_maxC-_minC)*(typename EigenTraits<PT>::ScalarType)(eps*0.5f);
  return createMM(_minC-d,_maxC+d);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::enlarged(T len,const sizeType d) {
  for(sizeType i=0; i<d; i++) {
    _minC[i]-=len;
    _maxC[i]+=len;
  }
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox BBox<T,DIM>::enlarge(T len,const sizeType d) const {
  BBox ret=createMM(_minC,_maxC);
  ret.enlarged(len,d);
  return ret;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::PT BBox<T,DIM>::lerp(const PT& frac) const {
  return (_maxC.array()*frac.array()-_minC.array()*(frac.array()-(typename EigenTraits<PT>::ScalarType)1.0f)).matrix();
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::empty() const {
  return !compL(_minC,_maxC);
}
template <typename T,int DIM>
template <int DIM2>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::containDim(const PT& point) const {
  for(int i=0; i<DIM2; i++)
    if(_minC[i] > point[i] || _maxC[i] < point[i])
      return false;
  return true;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::contain(const BBox& other,const sizeType d) const {
  for(int i=0; i<d; i++)
    if(_minC[i] > other._minC[i] || _maxC[i] < other._maxC[i])
      return false;
  return true;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::contain(const PT& point,const sizeType d) const {
  for(int i=0; i<d; i++)
    if(_minC[i] > point[i] || _maxC[i] < point[i])
      return false;
  return true;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::contain(const PT& point,const T& rad,const sizeType d) const {
  for(int i=0; i<d; i++)
    if(_minC[i]+rad > point[i] || _maxC[i]-rad < point[i])
      return false;
  return true;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC void BBox<T,DIM>::reset() {
  _minC=PT::Constant( ScalarUtil<T>::scalar_max());
  _maxC=PT::Constant(-ScalarUtil<T>::scalar_max());
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::PT BBox<T,DIM>::getExtent() const {
  return _maxC-_minC;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC T BBox<T,DIM>::distTo(const BBox& other,const sizeType d) const {
  PT dist=PT::Zero();
  for(sizeType i=0; i<d; i++) {
    if (other._maxC[i] < _minC[i])
      dist[i] = other._maxC[i] - _minC[i];
    else if (other._minC[i] > _maxC[i])
      dist[i] = other._minC[i] - _maxC[i];
  }
  return dist.norm();
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC T BBox<T,DIM>::distTo(const PT& pt,const sizeType d) const {
  return std::sqrt(distToSqr(pt,d));
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC T BBox<T,DIM>::distToSqr(const PT& pt,const sizeType d) const {
  PT dist=PT::Zero();
  for(sizeType i=0; i<d; i++) {
    if (pt[i] < _minC[i])
      dist[i] = pt[i] - _minC[i];
    else if (pt[i] > _maxC[i])
      dist[i] = pt[i] - _maxC[i];
  }
  return dist.squaredNorm();
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::PT BBox<T,DIM>::closestTo(const PT& pt,const sizeType d) const {
  PT dist(pt);
  for(sizeType i=0; i<d; i++) {
    if (pt[i] < _minC[i])
      dist[i] = _minC[i];
    else if (pt[i] > _maxC[i])
      dist[i] = _maxC[i];
  }
  return dist;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::intersect(const PT& p,const PT& q,const sizeType d) const {
  T s=0, t=1;
  return intersect(p,q,s,t,d);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::intersect(const PT& p,const PT& q,T& s,T& t,const sizeType d) const {
  const T lo=1-5*ScalarUtil<T>::scalar_eps();
  const T hi=1+5*ScalarUtil<T>::scalar_eps();

  s=0;
  t=1;
  for(sizeType i=0; i<d; ++i) {
    T D=q[i]-p[i];
    if(p[i]<q[i]) {
      T s0=lo*(_minC[i]-p[i])/D, t0=hi*(_maxC[i]-p[i])/D;
      if(s0>s) s=s0;
      if(t0<t) t=t0;
    } else if(p[i]>q[i]) {
      T s0=lo*(_maxC[i]-p[i])/D, t0=hi*(_minC[i]-p[i])/D;
      if(s0>s) s=s0;
      if(t0<t) t=t0;
    } else {
      if(p[i]<_minC[i] || p[i]>_maxC[i])
        return false;
    }

    if(s>t)
      return false;
  }
  return true;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC bool BBox<T,DIM>::intersect(const BBox& other,const sizeType& d) const {
  for(sizeType i=0; i<d; i++) {
    if(_maxC[i] < other._minC[i] || other._maxC[i] < _minC[i])
      return false;
  }
  return true;
  //return compLE(_minC,other._maxC) && compLE(other._minC,_maxC);
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::PT2 BBox<T,DIM>::project(const PT& a,const sizeType d) const {
  PT ctr=(_minC+_maxC)*(typename EigenTraits<PT>::ScalarType)0.5f;
  T ctrD=a.dot(ctr);
  T delta=0.0f;
  ctr=_maxC-ctr;
  for(sizeType i=0; i<d; i++)
    delta+=std::abs(ctr[i]*a[i]);
  return PT2(ctrD-delta,ctrD+delta);
}
template <typename T,int DIM>
template<typename T2>
EIGEN_DEVICE_FUNC typename BBox<T,DIM>::BBox& BBox<T,DIM>::copy(const BBox<T2,DIM>& other) {
  for(sizeType i=0; i<DIM; i++) {
    _minC[i]=std::convert<T>()(other._minC[i]);
    _maxC[i]=std::convert<T>()(other._maxC[i]);
  }
  return *this;
}
template <typename T,int DIM>
EIGEN_DEVICE_FUNC T BBox<T,DIM>::perimeter(const sizeType d) const {
  PT ext=getExtent();
  if(d <= 2)
    return ext.sum()*2.0f;
  else {
    ASSERT(d == 3);
    return (ext[0]*ext[1]+ext[1]*ext[2]+ext[0]*ext[2])*2.0f;
  }
}
//LineSegTpl
template <typename T>
EIGEN_DEVICE_FUNC LineSegTpl<T>::LineSegTpl() {}
template <typename T>
EIGEN_DEVICE_FUNC LineSegTpl<T>::LineSegTpl(const PT& x,const PT& y):_x(x),_y(y) {}
template <typename T>
bool LineSegTpl<T>::write(std::ostream& os) const
{
  writeBinaryData(_x,os);
  writeBinaryData(_y,os);
  return os.good();
}
template <typename T>
bool LineSegTpl<T>::read(std::istream& is)
{
  readBinaryData(_x,is);
  readBinaryData(_y,is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC T LineSegTpl<T>::length() const
{
  return (_y-_x).norm();
}
template <typename T>
EIGEN_DEVICE_FUNC typename LineSegTpl<T>::PT LineSegTpl<T>::circumcenter() const
{
  return (_x+_y)*(T)0.5f;
}
template <typename T>
EIGEN_DEVICE_FUNC typename LineSegTpl<T>::PT LineSegTpl<T>::masscenter() const
{
  return (_x+_y)*(T)0.5f;
}
template <typename T>
EIGEN_DEVICE_FUNC typename LineSegTpl<T>::PT LineSegTpl<T>::normal() const
{
  PT ret=_y-_x;
  return PT(ret.y(),-ret.x(),0.0f)/
         max(ret.norm(),ScalarUtil<T>::scalar_eps());
}
template <typename T>
EIGEN_DEVICE_FUNC typename LineSegTpl<T>::PT LineSegTpl<T>::gradDir(sizeType id) const
{
  PT dir=_x-_y;
  dir/=dir.squaredNorm();
  if(id == 1)
    dir*=(T)-1.0f;
  return dir;
}
template <typename T>
EIGEN_DEVICE_FUNC T LineSegTpl<T>::signedArea() const
{
  return _x.cross(_y)[2]/2;
}
template <typename T>
EIGEN_DEVICE_FUNC bool LineSegTpl<T>::intersect(const LineSegTpl<T>& l,T& t,bool infinite) const
{
  PT2 abt;
  MAT2 mat;
  mat.col(0)=(_y-_x).block(0,0,2,1);
  mat.col(1)=(l._x-l._y).block(0,0,2,1);
  if(std::abs(mat.determinant()) < ScalarUtil<T>::scalar_eps())
    return false;

  abt=InvTraits<MAT2,Eigen::NumTraits<T>::IsInteger>::inverse(mat)*(l._x-_x).block(0,0,2,1);
  t=abt.y();
  if(infinite)
    abt.y()=min(abt.y(),(T)0.5f);
  return (abt.x() >= 0.0f && abt.x() <= 1.0f) &&	//in current LineSeg
         (abt.y() >= 0.0f && abt.y() <= 1.0f);	//in given LineSeg
}
template <typename T>
EIGEN_DEVICE_FUNC void LineSegTpl<T>::calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT& b) const
{
  b[1]=(pt-_x).dot(_y-_x)/max((_y-_x).squaredNorm(),ScalarUtil<T>::scalar_eps());
  b[1]=max((T)0.0f,min((T)1.0f,b[1]));
  b[0]=1.0f-b[1];
  cp=_x*b[0]+_y*b[1];
  sqrDistance=(cp-pt).squaredNorm();
}
template <typename T>
EIGEN_DEVICE_FUNC void LineSegTpl<T>::calcLineDist(const LineSegTpl<T>& l,T& sqrDistance,T& a,T& b) const
{
#define P0 _x
#define P1 _y
#define Q0 l._x
#define Q1 l._y
  // The code allows degenerate line segments; that is, P0 and P1 can be
  // the same point or Q0 and Q1 can be the same point.  The quadratic
  // function for squared distance between the segment is
  //   R(s,t) = a*s^2 - 2*b*s*t + c*t^2 + 2*d*s - 2*e*t + f
  // for (s,t) in [0,1]^2 where
  //   a = Dot(P1-P0,P1-P0), b = Dot(P1-P0,Q1-Q0), c = Dot(Q1-Q0,Q1-Q0),
  //   d = Dot(P1-P0,P0-Q0), e = Dot(Q1-Q0,P0-Q0), f = Dot(P0-Q0,P0-Q0)
  PT P1mP0 = P1 - P0;
  PT Q1mQ0 = Q1 - Q0;
  PT P0mQ0 = P0 - Q0;
  T mA = P1mP0.dot(P1mP0);
  T mB = P1mP0.dot(Q1mQ0);
  T mC = Q1mQ0.dot(Q1mQ0);
  T mD = P1mP0.dot(P0mQ0);
  T mE = Q1mQ0.dot(P0mQ0);

  T mF00 = mD;
  T mF10 = mF00 + mA;
  T mF01 = mF00 - mB;
  T mF11 = mF10 - mB;

  T mG00 = -mE;
  T mG10 = mG00 - mB;
  T mG01 = mG00 + mC;
  T mG11 = mG10 + mC;

  if (mA > (T)0 && mC > (T)0)
  {
    // Compute the solutions to dR/ds(s0,0) = 0 and dR/ds(s1,1) = 0.  The
    // location of sI on the s-axis is stored in classifyI (I = 0 or 1).  If
    // sI <= 0, classifyI is -1.  If sI >= 1, classifyI is 1.  If 0 < sI < 1,
    // classifyI is 0.  This information helps determine where to search for
    // the minimum point (s,t).  The fij values are dR/ds(i,j) for i and j in
    // {0,1}.

    T sValue[2];
    sValue[0] = getClampedRoot(mA, mF00, mF10);
    sValue[1] = getClampedRoot(mA, mF01, mF11);

    sizeType classify[2];
    for (sizeType i = 0; i < 2; ++i)
    {
      if (sValue[i] <= (T)0)
      {
        classify[i] = -1;
      }
      else if (sValue[i] >= (T)1)
      {
        classify[i] = +1;
      }
      else
      {
        classify[i] = 0;
      }
    }

    if (classify[0] == -1 && classify[1] == -1)
    {
      // The minimum must occur on s = 0 for 0 <= t <= 1.
      a = (T)0;
      b = getClampedRoot(mC, mG00, mG01);
    }
    else if (classify[0] == +1 && classify[1] == +1)
    {
      // The minimum must occur on s = 1 for 0 <= t <= 1.
      a = (T)1;
      b = getClampedRoot(mC, mG10, mG11);
    }
    else
    {
      // The line dR/ds = 0 intersects the domain [0,1]^2 in a
      // nondegenerate segment.  Compute the endpoints of that segment,
      // end[0] and end[1].  The edge[i] flag tells you on which domain
      // edge end[i] lives: 0 (s=0), 1 (s=1), 2 (t=0), 3 (t=1).
      sizeType edge[2];
      T end[2][2];
      computeIntersection(mF00, mF10, mB, sValue, classify, edge, end);

      // The directional derivative of R along the segment of
      // intersection is
      //   H(z) = (end[1][1]-end[1][0])*dR/dt((1-z)*end[0] + z*end[1])
      // for z in [0,1].  The formula uses the fact that dR/ds = 0 on
      // the segment.  Compute the minimum of H on [0,1].
      T parameter[2];
      computeMinimumParameters(mB, mC, mE, mG00, mG01, mG10, mG11, edge, end, parameter);
      a=parameter[0];
      b=parameter[1];
    }
  }
  else
  {
    if (mA > (T)0)
    {
      // The Q-segment is degenerate (Q0 and Q1 are the same point) and
      // the quadratic is R(s,0) = a*s^2 + 2*d*s + f and has (half)
      // first derivative F(t) = a*s + d.  The closest P-point is
      // interior to the P-segment when F(0) < 0 and F(1) > 0.
      a = getClampedRoot(mA, mF00, mF10);
      b = (T)0;
    }
    else if (mC > (T)0)
    {
      // The P-segment is degenerate (P0 and P1 are the same point) and
      // the quadratic is R(0,t) = c*t^2 - 2*e*t + f and has (half)
      // first derivative G(t) = c*t - e.  The closest Q-point is
      // interior to the Q-segment when G(0) < 0 and G(1) > 0.
      a = (T)0;
      b = getClampedRoot(mC, mG00, mG01);
    }
    else
    {
      // P-segment and Q-segment are degenerate.
      a = (T)0;
      b = (T)0;
    }
  }

  PT diff=_x*(1-a)+_y*a;
  diff-=l._x*(1-b)+l._y*b;
  sqrDistance=diff.squaredNorm();
#undef P0
#undef P1
#undef Q0
#undef Q1
}
template <typename T>
EIGEN_DEVICE_FUNC void LineSegTpl<T>::writeVTK(VTKWriter<T>& os) const
{
  os.setRelativeIndex();
  std::vector<PT,Eigen::aligned_allocator<PT> > vss;
  vss.push_back(_x);
  vss.push_back(_y);
  os.appendPoints(vss.begin(),vss.end());

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
  iss.push_back(Vec3i(0,1,0));
  os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
}
template <typename T>
EIGEN_DEVICE_FUNC T LineSegTpl<T>::getClampedRoot(T slope,T h0,T h1) const
{
  // Theoretically, r is in (0,1).  However, when the slope is nearly zero,
  // then so are h0 and h1.  Significant numerical rounding problems can
  // occur when using floating-point arithmetic.  If the rounding causes r
  // to be outside the interval, clamp it.  It is possible that r is in
  // (0,1) and has rounding errors, but because h0 and h1 are both nearly
  // zero, the quadratic is nearly constant on (0,1).  Any choice of p
  // should not cause undesirable accuracy problems for the final distance
  // computation.
  //
  // NOTE:  You can use bisection to recompute the root or even use
  // bisection to compute the root and skip the division.  This is generally
  // slower, which might be a problem for high-performance applications.

  T r;
  if (h0 < (T)0)
  {
    if (h1 > (T)0)
    {
      r = -h0 / slope;
      if (r > (T)1)
      {
        r = (T)0.5;
      }
      // The slope is positive and -h0 is positive, so there is no
      // need to test for a negative value and clamp it.
    }
    else
    {
      r = (T)1;
    }
  }
  else
  {
    r = (T)0;
  }
  return r;
}
template <typename T>
EIGEN_DEVICE_FUNC void LineSegTpl<T>::computeIntersection(T mF00,T mF10,T mB,const T sValue[2],const sizeType classify[2],sizeType edge[2],T end[2][2]) const
{
  // The divisions are theoretically numbers in [0,1].  Numerical rounding
  // errors might cause the result to be outside the interval.  When this
  // happens, it must be that both numerator and denominator are nearly
  // zero.  The denominator is nearly zero when the segments are nearly
  // perpendicular.  The numerator is nearly zero when the P-segment is
  // nearly degenerate (mF00 = a is small).  The choice of 0.5 should not
  // cause significant accuracy problems.
  //
  // NOTE:  You can use bisection to recompute the root or even use
  // bisection to compute the root and skip the division.  This is generally
  // slower, which might be a problem for high-performance applications.

  if (classify[0] < 0)
  {
    edge[0] = 0;
    end[0][0] = (T)0;
    end[0][1] = mF00 / mB;
    if (end[0][1] < (T)0 || end[0][1] > (T)1)
    {
      end[0][1] = (T)0.5;
    }

    if (classify[1] == 0)
    {
      edge[1] = 3;
      end[1][0] = sValue[1];
      end[1][1] = (T)1;
    }
    else  // classify[1] > 0
    {
      edge[1] = 1;
      end[1][0] = (T)1;
      end[1][1] = mF10 / mB;
      if (end[1][1] < (T)0 || end[1][1] > (T)1)
      {
        end[1][1] = (T)0.5;
      }
    }
  }
  else if (classify[0] == 0)
  {
    edge[0] = 2;
    end[0][0] = sValue[0];
    end[0][1] = (T)0;

    if (classify[1] < 0)
    {
      edge[1] = 0;
      end[1][0] = (T)0;
      end[1][1] = mF00 / mB;
      if (end[1][1] < (T)0 || end[1][1] > (T)1)
      {
        end[1][1] = (T)0.5;
      }
    }
    else if (classify[1] == 0)
    {
      edge[1] = 3;
      end[1][0] = sValue[1];
      end[1][1] = (T)1;
    }
    else
    {
      edge[1] = 1;
      end[1][0] = (T)1;
      end[1][1] = mF10 / mB;
      if (end[1][1] < (T)0 || end[1][1] > (T)1)
      {
        end[1][1] = (T)0.5;
      }
    }
  }
  else  // classify[0] > 0
  {
    edge[0] = 1;
    end[0][0] = (T)1;
    end[0][1] = mF10 / mB;
    if (end[0][1] < (T)0 || end[0][1] > (T)1)
    {
      end[0][1] = (T)0.5;
    }

    if (classify[1] == 0)
    {
      edge[1] = 3;
      end[1][0] = sValue[1];
      end[1][1] = (T)1;
    }
    else
    {
      edge[1] = 0;
      end[1][0] = (T)0;
      end[1][1] = mF00 / mB;
      if (end[1][1] < (T)0 || end[1][1] > (T)1)
      {
        end[1][1] = (T)0.5;
      }
    }
  }
}
template <typename T>
EIGEN_DEVICE_FUNC void LineSegTpl<T>::computeMinimumParameters(T mB,T mC,T mE,T mG00,T mG01,T mG10,T mG11,const sizeType edge[2],const T end[2][2],T parameter[2]) const
{
  T delta = end[1][1] - end[0][1];
  T h0 = delta * (-mB * end[0][0] + mC * end[0][1] - mE);
  if (h0 >= (T)0)
  {
    if (edge[0] == 0)
    {
      parameter[0] = (T)0;
      parameter[1] = getClampedRoot(mC, mG00, mG01);
    }
    else if (edge[0] == 1)
    {
      parameter[0] = (T)1;
      parameter[1] = getClampedRoot(mC, mG10, mG11);
    }
    else
    {
      parameter[0] = end[0][0];
      parameter[1] = end[0][1];
    }
  }
  else
  {
    T h1 = delta * (-mB * end[1][0] + mC * end[1][1] - mE);
    if (h1 <= (T)0)
    {
      if (edge[1] == 0)
      {
        parameter[0] = (T)0;
        parameter[1] = getClampedRoot(mC, mG00, mG01);
      }
      else if (edge[1] == 1)
      {
        parameter[0] = (T)1;
        parameter[1] = getClampedRoot(mC, mG10, mG11);
      }
      else
      {
        parameter[0] = end[1][0];
        parameter[1] = end[1][1];
      }
    }
    else  // h0 < 0 and h1 > 0
    {
      T z = std::min<T>(std::max<T>(h0 / (h0 - h1), (T)0), (T)1);
      T omz = (T)1 - z;
      parameter[0] = omz * end[0][0] + z * end[1][0];
      parameter[1] = omz * end[0][1] + z * end[1][1];
    }
  }
}
//PlaneTpl
template <typename T>
EIGEN_DEVICE_FUNC PlaneTpl<T>::PlaneTpl() {}
template <typename T>
EIGEN_DEVICE_FUNC PlaneTpl<T>::PlaneTpl(const PT& x0,const PT& n):_x0(x0),_n(n) {}
template <typename T>
EIGEN_DEVICE_FUNC PlaneTpl<T>::PlaneTpl(const PT& a,const PT& b,const PT& c):_x0(a),_n((b-a).cross(c-a)) {}
template <typename T>
bool PlaneTpl<T>::write(std::ostream& os) const
{
  writeBinaryData(_x0,os);
  writeBinaryData(_n,os);
  return os.good();
}
template <typename T>
bool PlaneTpl<T>::read(std::istream& is)
{
  readBinaryData(_x0,is);
  readBinaryData(_n,is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC T PlaneTpl<T>::side(const PT& p) const
{
  return (p-_x0).dot(_n);
}
template <typename T>
EIGEN_DEVICE_FUNC bool PlaneTpl<T>::intersect(const BBox<T>& bb) const
{
  return SAT<T>::testSegment(PT2::Constant(_x0.dot(_n)),bb.project(_n));
}
template <typename T>
EIGEN_DEVICE_FUNC typename PlaneTpl<T>::PT2 PlaneTpl<T>::project(const PT& d) const
{
  return (d==_n) ? PT2::Constant(_x0.dot(_n)) :
         PT2(-ScalarUtil<T>::scalar_max(),ScalarUtil<T>::scalar_max());
}
template <typename T>
EIGEN_DEVICE_FUNC void PlaneTpl<T>::writeVTK(VTKWriter<T>& os) const
{
#define SEG 128
  {
    sizeType mid;
    _n.maxCoeff(&mid);
    PT nx=(mid == 0 || mid == 1) ? PT(-_n[1],_n[0],0.0f) : PT(-_n[2],0.0f,_n[0]);
    nx.normalize();
    PT ny=_n.cross(nx).normalized();
    os.setRelativeIndex();
    std::vector<PT,Eigen::aligned_allocator<PT> > vss;
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
    for(sizeType i=0; i<SEG; i++) {
      scalar ANG=2.0f*M_PI*(scalar)i/(scalar)SEG;
      vss.push_back(nx*(T)cos(ANG)+ny*(T)sin(ANG));
      iss.push_back(Vec3i(i,(i+1)%SEG,0));
    }
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
  }
#undef SEG

  {
    os.setRelativeIndex();
    std::vector<PT,Eigen::aligned_allocator<PT> > vss;
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
    vss.push_back(_x0);
    vss.push_back(_n+_x0);
    iss.push_back(Vec3i(0,1,0));
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
  }
}
//TriangleTpl
template <typename T>
EIGEN_DEVICE_FUNC TriangleTpl<T>::TriangleTpl() {}
template <typename T>
EIGEN_DEVICE_FUNC TriangleTpl<T>::TriangleTpl(const PT& a,const PT& b,const PT& c):_a(a),_b(b),_c(c) {}
template <typename T>
bool TriangleTpl<T>::write(std::ostream& os) const
{
  writeBinaryData(_a,os);
  writeBinaryData(_b,os);
  writeBinaryData(_c,os);
  return os.good();
}
template <typename T>
bool TriangleTpl<T>::read(std::istream& is)
{
  readBinaryData(_a,is);
  readBinaryData(_b,is);
  readBinaryData(_c,is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC typename TriangleTpl<T>::PT TriangleTpl<T>::circumcenter() const
{
  PT alpha=_b-_a;
  PT beta=_c-_a;

  PT2 ls(0.5f,0.0f);
  PT2 le(-alpha.dot(beta)/alpha.dot(alpha),1.0f);

  PT2 rs(0.0f,0.5f);
  PT2 re(1.0f,-alpha.dot(beta)/beta.dot(beta));

  MAT2 A;
  A.col(0)=le;
  A.col(1)=-re;
  T t=(InvTraits<MAT2,Eigen::NumTraits<T>::IsInteger>::inverse(A)*(rs-ls))[0];
  ls+=le*t;

  return _a+alpha*ls.x()+beta*ls.y();
}
template <typename T>
EIGEN_DEVICE_FUNC typename TriangleTpl<T>::PT TriangleTpl<T>::masscenter() const
{
  return (_a+_b+_c)/(T)3.0f;
}
template <typename T>
EIGEN_DEVICE_FUNC typename TriangleTpl<T>::PT TriangleTpl<T>::bary(const PT& pt) const
{
  Eigen::Matrix<T,3,2> A;
  A.col(0)=_a-_c;
  A.col(1)=_b-_c;
  MAT2 ATA=A.transpose()*A;
  PT2 ab=InvTraits<MAT2,Eigen::NumTraits<T>::IsInteger>::inverse(ATA)*(A.transpose()*(pt-_c));
  return PT(ab.x(),ab.y(),1.0f-ab.sum());
}
template <typename T>
EIGEN_DEVICE_FUNC T TriangleTpl<T>::area() const
{
  return (_b-_a).cross(_c-_a).norm()*0.5f;
}
template <typename T>
EIGEN_DEVICE_FUNC T TriangleTpl<T>::signedVolume() const
{
  T v321 = _c.x() * _b.y() * _a.z();
  T v231 = _b.x() * _c.y() * _a.z();
  T v312 = _c.x() * _a.y() * _b.z();
  T v132 = _a.x() * _c.y() * _b.z();
  T v213 = _b.x() * _a.y() * _c.z();
  T v123 = _a.x() * _b.y() * _c.z();
  return (1.0f/6.0f)*(-v321+v231+v312-v132-v213+v123);
}
template <typename T>
EIGEN_DEVICE_FUNC typename TriangleTpl<T>::PT TriangleTpl<T>::normal() const
{
  PT dir=(_b-_a).cross(_c-_a);
  return dir/max(dir.norm(),ScalarUtil<T>::scalar_eps());
}
template <typename T>
EIGEN_DEVICE_FUNC typename TriangleTpl<T>::PT TriangleTpl<T>::gradDir(sizeType id) const
{
  switch(id) {
  case 0:
    return height(_a,_b,_c);
  case 1:
    return height(_b,_a,_c);
  default:
    return height(_c,_a,_b);
  }
}
template <typename T>
EIGEN_DEVICE_FUNC typename TriangleTpl<T>::PT TriangleTpl<T>::height(const PT& a,const PT& b,const PT& c) const
{
  PT ret=a-b;
  PT bc=(c-b).normalized();
  ret=ret-ret.dot(bc)*bc;
  return ret/ret.squaredNorm();
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::isInside(const PT& pt) const
{
  PT bc=bary(pt);
  PT ptPlane=bc.x()*_a+bc.y()*_b+bc.z()*_c;
  return (ptPlane-pt).norm() < ScalarUtil<T>::scalar_eps() &&
         bc.x() >= 0.0f && bc.y() >= 0.0f && bc.z() >= 0.0f;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::isPlanePointInside(const PT& pt) const
{
  PT bc=bary(pt);
  return bc.x() >= 0.0f && bc.y() >= 0.0f && bc.z() >= 0.0f;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::updateIntr(T& s,T& t,T num,T denom) const
{
  if(std::abs(denom) < ScalarUtil<T>::scalar_eps())
    return num <= 0.0f;
  else if(denom > 0.0f)
    s=max(s,T(num/denom));
  else
    t=min(t,T(num/denom));
  return s<t;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::intersect(const LineSegTpl<T>& l,T& s,T& t) const
{
  //compute A,B
  MAT2 m;
  m(0,0)=(_a-_c).x();
  m(0,1)=(_b-_c).x();
  m(1,0)=(_a-_c).y();
  m(1,1)=(_b-_c).y();
  PT2 A=InvTraits<MAT2,Eigen::NumTraits<T>::IsInteger>::inverse(m)*PT2(l._x.x()-_c.x(),l._x.y()-_c.y());
  PT2 B=InvTraits<MAT2,Eigen::NumTraits<T>::IsInteger>::inverse(m)*PT2(l._y.x()-l._x.x(),l._y.y()-l._x.y());

  //get the expected interval
  s=0.0f,t=1.0f;
  return updateIntr(s,t,-A.x(),B.x()) &&
         updateIntr(s,t,-A.y(),B.y()) &&
         updateIntr(s,t,A.x()+A.y()-1.0f,-(B.x()+B.y()));
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::intersect(const LineSegTpl<T>& l,T& t,bool infinite,PT* abtOut) const
{
  PT abt;
  MAT3 mat;
  mat.col(0)=_a-_c;
  mat.col(1)=_b-_c;
  mat.col(2)=l._x-l._y;
  if(std::abs(mat.determinant()) < ScalarUtil<T>::scalar_eps())
    return false;

  abt=InvTraits<MAT3,Eigen::NumTraits<T>::IsInteger>::inverse(mat)*(l._x-_c);
  t=abt.z();
  if(abtOut)
    *abtOut=abt;
  if(infinite)
    abt.z()=min(abt.z(),(T)0.5f);
  return (abt.x() >= 0.0f && abt.y() >= 0.0f && (abt.x()+abt.y()) <= 1.0f) &&	//in triangle
         (abt.z() >= 0.0f && abt.z() <= 1.0f);	//in segment
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::intersect(const TriangleTpl<T>& t) const
{
  PT n=normal();
  PT nt=t.normal();
  if(n.norm() < 0.9 || nt.norm() < 0.9)
    return false;
  else if(n.cross(nt).norm() < 1E-6f)
    return false;

  T b;
  if(intersect(LineSegTpl<T>(t._a,t._b),b,false) ||
      intersect(LineSegTpl<T>(t._b,t._c),b,false) ||
      intersect(LineSegTpl<T>(t._c,t._a),b,false))
    return true;
  if(t.intersect(LineSegTpl<T>(_a,_b),b,false) ||
      t.intersect(LineSegTpl<T>(_b,_c),b,false) ||
      t.intersect(LineSegTpl<T>(_c,_a),b,false))
    return true;
  return false;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::intersect(const BBox<T>& bb) const
{
  PT d=(_b-_a).cross(_c-_a);
  if(SAT<T>::testSegment(bb.project(d),PT2::Constant(_a.dot(d))))
    return false;
  for(sizeType i=0; i<3; i++)
    if(SAT<T>::sep(PT::Unit(i),*this,bb))
      return false;
  for(sizeType i=0; i<3; i++)
    if(SAT<T>::sep(PT::Unit(i).cross(d),*this,bb))
      return false;
  return true;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::calcLineDist(const LineSegTpl<T>& l,T& sqrDistance,PT& bt,PT2& bl) const
{
  T t,t2;
  PT abt,cp;
  if(intersect(l,t,false,&abt)) {
    sqrDistance=0;
    bt=PT(abt[0],abt[1],1-abt[0]-abt[1]);
    bl=PT2(1-t,t);
    return true;
  } else {
    T sqrDist;
    sqrDistance=ScalarUtil<T>::scalar_max();
    //LL
    for(sizeType d=0; d<3; d++) {
      LineSegTpl<T> l2((*this)[d],(*this)[(d+1)%3]);
      l2.calcLineDist(l,sqrDist,t,t2);
      if(sqrDist<sqrDistance) {
        sqrDistance=sqrDist;
        bt.setZero();
        bt[d]=1-t;
        bt[(d+1)%3]=t;
        bl=PT2(1-t2,t2);
      }
    }
    //X
    calcPointDist(l._x,sqrDist,cp,abt);
    if(sqrDist<sqrDistance) {
      sqrDistance=sqrDist;
      bt=abt;
      bl=PT2(1,0);
    }
    //Y
    calcPointDist(l._y,sqrDist,cp,abt);
    if(sqrDist<sqrDistance) {
      sqrDistance=sqrDist;
      bt=abt;
      bl=PT2(0,1);
    }
    return false;
  }
}
template <typename T>
template <typename PT_BARY>
EIGEN_DEVICE_FUNC void TriangleTpl<T>::calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT_BARY& b) const
{
  PT diff = _a - pt;
  PT edge0 = _b - _a;
  PT edge1 = _c - _a;
  T a00 = edge0.dot(edge0);
  T a01 = edge0.dot(edge1);
  T a11 = edge1.dot(edge1);
  T b0 = diff.dot(edge0);
  T b1 = diff.dot(edge1);
  T c = diff.dot(diff);
  T det = std::abs(a00*a11 - a01*a01);
  T s = a01*b1 - a11*b0;
  T t = a01*b0 - a00*b1;

  if (s + t <= det) {
    if (s < (T)0) {
      if (t < (T)0) { // region 4
        if (b0 < (T)0) {
          t = (T)0;
          if (-b0 >= a00) {
            s = (T)1;
            sqrDistance = a00 + ((T)2)*b0 + c;
          } else {
            s = -b0/a00;
            sqrDistance = b0*s + c;
          }
        } else {
          s = (T)0;
          if (b1 >= (T)0) {
            t = (T)0;
            sqrDistance = c;
          } else if (-b1 >= a11) {
            t = (T)1;
            sqrDistance = a11 + ((T)2)*b1 + c;
          } else {
            t = -b1/a11;
            sqrDistance = b1*t + c;
          }
        }
      } else { // region 3
        s = (T)0;
        if (b1 >= (T)0) {
          t = (T)0;
          sqrDistance = c;
        } else if (-b1 >= a11) {
          t = (T)1;
          sqrDistance = a11 + ((T)2)*b1 + c;
        } else {
          t = -b1/a11;
          sqrDistance = b1*t + c;
        }
      }
    } else if (t < (T)0) { // region 5
      t = (T)0;
      if (b0 >= (T)0) {
        s = (T)0;
        sqrDistance = c;
      } else if (-b0 >= a00) {
        s = (T)1;
        sqrDistance = a00 + ((T)2)*b0 + c;
      } else {
        s = -b0/a00;
        sqrDistance = b0*s + c;
      }
    } else { // region 0
      // minimum at interior point
      T invDet = ((T)1)/det;
      s *= invDet;
      t *= invDet;
      sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                    t*(a01*s + a11*t + ((T)2)*b1) + c;
    }
  } else {
    T tmp0, tmp1, numer, denom;

    if (s < (T)0) { // region 2
      tmp0 = a01 + b0;
      tmp1 = a11 + b1;
      if (tmp1 > tmp0) {
        numer = tmp1 - tmp0;
        denom = a00 - ((T)2)*a01 + a11;
        if (numer >= denom) {
          s = (T)1;
          t = (T)0;
          sqrDistance = a00 + ((T)2)*b0 + c;
        } else {
          s = numer/denom;
          t = (T)1 - s;
          sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                        t*(a01*s + a11*t + ((T)2)*b1) + c;
        }
      } else {
        s = (T)0;
        if (tmp1 <= (T)0) {
          t = (T)1;
          sqrDistance = a11 + ((T)2)*b1 + c;
        } else if (b1 >= (T)0) {
          t = (T)0;
          sqrDistance = c;
        } else {
          t = -b1/a11;
          sqrDistance = b1*t + c;
        }
      }
    } else if (t < (T)0) { // region 6
      tmp0 = a01 + b1;
      tmp1 = a00 + b0;
      if (tmp1 > tmp0) {
        numer = tmp1 - tmp0;
        denom = a00 - ((T)2)*a01 + a11;
        if (numer >= denom) {
          t = (T)1;
          s = (T)0;
          sqrDistance = a11 + ((T)2)*b1 + c;
        } else {
          t = numer/denom;
          s = (T)1 - t;
          sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                        t*(a01*s + a11*t + ((T)2)*b1) + c;
        }
      } else {
        t = (T)0;
        if (tmp1 <= (T)0) {
          s = (T)1;
          sqrDistance = a00 + ((T)2)*b0 + c;
        } else if (b0 >= (T)0) {
          s = (T)0;
          sqrDistance = c;
        } else {
          s = -b0/a00;
          sqrDistance = b0*s + c;
        }
      }
    } else { // region 1
      numer = a11 + b1 - a01 - b0;
      if (numer <= (T)0) {
        s = (T)0;
        t = (T)1;
        sqrDistance = a11 + ((T)2)*b1 + c;
      } else {
        denom = a00 - ((T)2)*a01 + a11;
        if (numer >= denom) {
          s = (T)1;
          t = (T)0;
          sqrDistance = a00 + ((T)2)*b0 + c;
        } else {
          s = numer/denom;
          t = (T)1 - s;
          sqrDistance = s*(a00*s + a01*t + ((T)2)*b0) +
                        t*(a01*s + a11*t + ((T)2)*b1) + c;
        }
      }
    }
  }

  // Account for numerical round-off error.
  if (sqrDistance < (T)0) {
    sqrDistance = (T)0;
  }

  cp = _a + s*edge0 + t*edge1;
  b(1) = s;
  b(2) = t;
  b(0) = (T)1 - s - t;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TriangleTpl<T>::calcTriangleDist(const TriangleTpl<T>& t2,T& sqrDistance,PT& bt,PT& bt2) const
{
  if(intersect(t2)) {
    sqrDistance=0;
    return true;
  } else {
    T sqrDist;
    PT btCurr;
    PT2 blCurr;
    sqrDistance=ScalarUtil<scalarD>::scalar_max();
    for(sizeType d=0; d<3; d++) {
      LineSegTpl<T> l(t2[d],t2[(d+1)%3]);
      calcLineDist(l,sqrDist,btCurr,blCurr);
      if(sqrDist < sqrDistance) {
        bt=btCurr;
        bt2.setZero();
        bt2[d]=blCurr[0];
        bt2[(d+1)%3]=blCurr[1];
        sqrDistance=sqrDist;
      }
    }
    for(sizeType d=0; d<3; d++) {
      LineSegTpl<T> l((*this)[d],(*this)[(d+1)%3]);
      t2.calcLineDist(l,sqrDist,btCurr,blCurr);
      if(sqrDist < sqrDistance) {
        bt2=btCurr;
        bt.setZero();
        bt[d]=blCurr[0];
        bt[(d+1)%3]=blCurr[1];
        sqrDistance=sqrDist;
      }
    }
    return false;
  }
}
template <typename T>
EIGEN_DEVICE_FUNC typename TriangleTpl<T>::PT2 TriangleTpl<T>::project(const PT& d) const
{
  T a=_a.dot(d);
  T b=_b.dot(d);
  T c=_c.dot(d);
  return PT2(min(a,min(b,c)),max(a,max(b,c)));
}
template <typename T>
EIGEN_DEVICE_FUNC void TriangleTpl<T>::writeVTK(VTKWriter<T>& os) const
{
  os.setRelativeIndex();
  std::vector<PT,Eigen::aligned_allocator<PT> > vss;
  vss.push_back(_a);
  vss.push_back(_b);
  vss.push_back(_c);
  os.appendPoints(vss.begin(),vss.end());

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
  iss.push_back(Vec3i(0,1,0));
  iss.push_back(Vec3i(1,2,0));
  iss.push_back(Vec3i(2,0,0));
  os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
}
template <typename T>
EIGEN_DEVICE_FUNC const typename TriangleTpl<T>::PT& TriangleTpl<T>::operator[](sizeType d) const
{
  return d==0?_a:d==1?_b:_c;
}
//TetrahedronTpl
template <typename T>
EIGEN_DEVICE_FUNC TetrahedronTpl<T>::TetrahedronTpl() {}
template <typename T>
EIGEN_DEVICE_FUNC TetrahedronTpl<T>::TetrahedronTpl(const PT& a,const PT& b,const PT& c,const PT& d):_a(a),_b(b),_c(c),_d(d)
{
  _swap=false;
  if(volume() < 0.0f) {
    std::swap(_c,_d);
    _swap=true;
  }
}
template <typename T>
bool TetrahedronTpl<T>::write(std::ostream& os) const
{
  writeBinaryData(_a,os);
  writeBinaryData(_b,os);
  writeBinaryData(_c,os);
  writeBinaryData(_d,os);
  return os.good();
}
template <typename T>
bool TetrahedronTpl<T>::read(std::istream& is)
{
  readBinaryData(_a,is);
  readBinaryData(_b,is);
  readBinaryData(_c,is);
  readBinaryData(_d,is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC typename TetrahedronTpl<T>::PT TetrahedronTpl<T>::circumcenter() const
{
  MAT3 A;
  A.row(0)=_b-_a;
  A.row(1)=_c-_a;
  A.row(2)=_d-_a;
  PT B=PT(A.row(0).squaredNorm(),
          A.row(1).squaredNorm(),
          A.row(2).squaredNorm())*(T)0.5f;
  MAT3 AI=InvTraits<MAT3,Eigen::NumTraits<T>::IsInteger>::inverse(A);
  return AI*B+_a;
}
template <typename T>
EIGEN_DEVICE_FUNC typename TetrahedronTpl<T>::PT TetrahedronTpl<T>::masscenter() const
{
  return (_a+_b+_c+_d)/(T)4.0f;
}
template <typename T>
EIGEN_DEVICE_FUNC typename TetrahedronTpl<T>::PT4 TetrahedronTpl<T>::bary(const PT& pt) const
{
  MAT3 A;
  A.col(0)=_a-_d;
  A.col(1)=_b-_d;
  A.col(2)=_c-_d;
  PT abc=InvTraits<MAT3,Eigen::NumTraits<T>::IsInteger>::inverse(A)*(pt-_d);
  return PT4(abc.x(),abc.y(),abc.z(),1.0f-abc.sum());
}
template <typename T>
EIGEN_DEVICE_FUNC bool TetrahedronTpl<T>::isInside(const PT& pt) const
{
  PT4 bc=bary(pt);
  return bc.x() >= 0 && bc.y() >= 0 && bc.z() >= 0 && bc.w() >= 0;
}
template <typename T>
EIGEN_DEVICE_FUNC T TetrahedronTpl<T>::volume() const
{
  return (_b-_a).cross(_c-_a).dot(_d-_a)/6.0f;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TetrahedronTpl<T>::dualCellVolume(PT4& dv)
{
  dv=PT4::Zero();
  PT cc=circumcenter();
  bool safe=isInside(cc);
  safe=safe&&dualFaceVolume(dv[0],dv[1],dv[2],_a,_b,_c,cc);
  safe=safe&&dualFaceVolume(dv[0],dv[2],dv[3],_a,_c,_d,cc);
  safe=safe&&dualFaceVolume(dv[0],dv[3],dv[1],_a,_d,_b,cc);
  safe=safe&&dualFaceVolume(dv[1],dv[2],dv[3],_b,_c,_d,cc);
  return safe;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TetrahedronTpl<T>::dualFaceVolume(T& va,T& vb,T& vc,const PT& a,const PT& b,const PT& c,const PT& cc)
{
  TriangleTpl<T> tri(a,b,c);
  PT ccf=tri.circumcenter();
  va+=TetrahedronTpl(cc,ccf,a,(a+b)*(T)0.5f).volume()+TetrahedronTpl(cc,ccf,a,(a+c)*(T)0.5f).volume();
  vb+=TetrahedronTpl(cc,ccf,b,(b+c)*(T)0.5f).volume()+TetrahedronTpl(cc,ccf,b,(b+a)*(T)0.5f).volume();
  vc+=TetrahedronTpl(cc,ccf,c,(c+a)*(T)0.5f).volume()+TetrahedronTpl(cc,ccf,c,(c+b)*(T)0.5f).volume();
  return tri.isPlanePointInside(ccf);
}
template <typename T>
EIGEN_DEVICE_FUNC T TetrahedronTpl<T>::dihedralAngle(const PT& a,const PT& b,const PT& c,const PT& d)
{
  PT n1=(a-c).cross(b-c).normalized();
  PT n2=(a-d).cross(b-d).normalized();
  return M_PI*0.5f-asin(n1.dot(n2));
}
template <typename T>
EIGEN_DEVICE_FUNC typename TetrahedronTpl<T>::PT6 TetrahedronTpl<T>::dihedralAngleTet()
{
  PT6 ret;
  ret[0]=dihedralAngle(_a,_b,_c,_d);
  ret[1]=dihedralAngle(_a,_c,_b,_d);
  ret[2]=dihedralAngle(_a,_d,_b,_c);
  ret[3]=dihedralAngle(_b,_c,_a,_d);
  ret[4]=dihedralAngle(_b,_d,_a,_c);
  ret[5]=dihedralAngle(_c,_d,_a,_b);
  return ret;
}
template <typename T>
EIGEN_DEVICE_FUNC bool TetrahedronTpl<T>::calcLineDist(const LineSegTpl<T>& l,PT4& bt,PT2& bl) const
{
  bt=bary(l._x);
  PT4 dbt=bary(l._y)-bt;

  T b=0.0f,t=1.0f;
  for(char i=0; i<4; i++) {
    if(bt[i] < 0.0f && dbt[i] < 0.0f)
      return false;
    if(bt[i] < 0.0f && dbt[i] > 0.0f)
      b=max((T)b,T((T)(-bt[i])/(T)  dbt[i] ));
    if(bt[i] > 0.0f && dbt[i] < 0.0f)
      t=min((T)t,T((T)  bt[i] /(T)(-dbt[i])));
  }
  if(b>t)return false;
  bt+=dbt*b;
  bl=PT2(1.0f-b,b);
  return true;
}
template <typename T>
EIGEN_DEVICE_FUNC void TetrahedronTpl<T>::calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT4& bc) const
{
  T sqrDistanceTmp;
  PT cpTmp;
  PT b;

  // Construct the planes for the faces of the tetrahedron.The normals
  // are outer pointing, but specified not to be unit length.We only need
  // to know sidedness of the query point, so we will save cycles by not
  // computing unit-length normals.
  PlaneTpl<T> plane;
  Vec3i indices;

  // Determine which faces are visible to the query point.Only these
  // need to be processed by point-to-triangle distance queries.
  sqrDistance=ScalarUtil<T>::scalar_max();
  //PT minTetraClosest=PT::Zero();
  for(sizeType i=0; i<4; ++i) {
    getPlane(i,plane,indices);
    if(plane.side(pt) >= 0) {
      TriangleTpl<T> tri(getNode(indices[0]),getNode(indices[1]),getNode(indices[2]));
      tri.calcPointDist(pt,sqrDistanceTmp,cpTmp,b);
      if (sqrDistanceTmp < sqrDistance) {
        sqrDistance=sqrDistanceTmp;
        cp=cpTmp;
        bc=PT4::Zero();
        bc[indices[0]]=b[0];
        bc[indices[1]]=b[1];
        bc[indices[2]]=b[2];
      }
    }
  }

  if (sqrDistance == ScalarUtil<T>::scalar_max()) {
    sqrDistance=(T)0.0f;
    cp=pt;
    bc=bary(pt);
  }
}
template <typename T>
EIGEN_DEVICE_FUNC void TetrahedronTpl<T>::getPlane(const sizeType& i,PlaneTpl<T>& p,Vec3i& ind) const
{
  //a,c,b
  //a,d,c
  //a,b,d
  //b,c,d
  switch(i) {
  case 0:
    p=PlaneTpl<T>(_a,_c,_b);
    ind=Vec3i(0,2,1);
    break;
  case 1:
    p=PlaneTpl<T>(_a,_d,_c);
    ind=Vec3i(0,3,2);
    break;
  case 2:
    p=PlaneTpl<T>(_a,_b,_d);
    ind=Vec3i(0,1,3);
    break;
  case 3:
    p=PlaneTpl<T>(_b,_c,_d);
    ind=Vec3i(1,2,3);
    break;
  }
}
template <typename T>
EIGEN_DEVICE_FUNC const typename TetrahedronTpl<T>::PT& TetrahedronTpl<T>::getNode(const sizeType& i) const
{
  return (&_a)[i];
}
template <typename T>
EIGEN_DEVICE_FUNC void TetrahedronTpl<T>::writeVTK(VTKWriter<T>& os) const
{
  os.setRelativeIndex();
  std::vector<PT,Eigen::aligned_allocator<PT> > vss;
  vss.push_back(_a);
  vss.push_back(_b);
  vss.push_back(_c);
  vss.push_back(_d);
  os.appendPoints(vss.begin(),vss.end());

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
  iss.push_back(Vec3i(0,1,0));
  iss.push_back(Vec3i(0,2,0));
  iss.push_back(Vec3i(0,3,0));
  iss.push_back(Vec3i(1,2,0));
  iss.push_back(Vec3i(2,3,0));
  iss.push_back(Vec3i(3,0,0));
  os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
}
//OBBTpl<T,2>
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,2>::OBBTpl() {}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,2>::OBBTpl(const BBox<T,2>& bb)
{
  _ext=bb.getExtent()*(T)0.5f;
  _rot.setIdentity();
  _trans=bb._maxC-_ext;
}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,2>::OBBTpl(const BBox<T,3>& bb)
{
  new (this)OBBTpl(BBox<T,2>(PT(bb._minC[0],bb._minC[1]),PT(bb._maxC[0],bb._maxC[1])));
}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,2>::OBBTpl(const MAT& rot,const PT& trans,const BBox<T,2>& bb)
{
  new (this)OBBTpl(bb);
  _rot=rot;
  _trans=_rot*_trans+trans;
}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,2>::OBBTpl(const MAT3& rot,const PT3& trans,const BBox<T,3>& bb)
{
  new (this)OBBTpl(rot.block(0,0,2,2),trans.block(0,0,2,1),BBox<T,2>(PT(bb._minC[0],bb._minC[1]),PT(bb._maxC[0],bb._maxC[1])));
}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,2>::OBBTpl(const MAT& rot,const PT& trans,const PT& ext):_rot(rot),_trans(trans),_ext(ext) {}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,2>::OBBTpl(const MAT3& rot,const PT3& trans,const PT3& ext)
{
  _rot=rot.template block<2,2>(0,0);
  _trans=trans.template block<2,1>(0,0);
  _ext=ext.template block<2,1>(0,0);
}
template <typename T>
bool OBBTpl<T,2>::write(std::ostream& os) const
{
  writeBinaryData(_rot,os);
  writeBinaryData(_trans,os);
  writeBinaryData(_ext,os);
  return os.good();
}
template <typename T>
bool OBBTpl<T,2>::read(std::istream& is)
{
  readBinaryData(_rot,is);
  readBinaryData(_trans,is);
  readBinaryData(_ext,is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,2>::closest(PT pt,PT& n,PT* normal) const
{
  pt=(_rot.transpose()*(pt-_trans)).eval();
  bool inside=closestInner(pt,n,normal);
  n=(_rot*n).eval();
  if(normal)
    *normal=(_rot**normal).eval();
  return inside;
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,2>::closestInner(const PT& pt,PT& n,PT* normal) const
{
  BBox<T,2> box(-_ext,_ext);
  if(box.template containDim<2>(pt.template segment<2>(0))) {
    T minDist=ScalarUtil<T>::scalar_max(),dist;
    for(sizeType d=0; d<2; d++) {
      dist=_ext[d]+pt[d];
      if(dist < minDist) {
        minDist=dist;
        n=-PT::Unit(d)*dist;
        if(normal)*normal=-PT::Unit(d);
      }
      dist=_ext[d]-pt[d];
      if(dist < minDist) {
        minDist=dist;
        n=PT::Unit(d)*dist;
        if(normal)*normal=PT::Unit(d);
      }
    }
    return true;
  } else {
    PT cp=box.closestTo(pt,2);
    n=cp-pt;
    if(normal)
      *normal=-n/max(n.norm(),(T)1E-6f);
    return false;
  }
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,2>::intersect(const OBBTpl<T,2>& other) const
{
  MAT rotB2A=_rot.transpose()*other._rot;
  PT transB2A=_rot.transpose()*(other._trans-_trans);
  PT dirX=rotB2A.col(0)*other._ext.x();
  PT dirY=rotB2A.col(1)*other._ext.y();
  PT dir(std::abs(dirX.x())+std::abs(dirY.x()),std::abs(dirX.y())+std::abs(dirY.y()));
  if(SAT<T>::testSegment(PT(-_ext.x(),_ext.x()),PT(transB2A.x()-dir.x(),transB2A.x()+dir.x())))
    return false;
  if(SAT<T>::testSegment(PT(-_ext.y(),_ext.y()),PT(transB2A.y()-dir.y(),transB2A.y()+dir.y())))
    return false;

  transB2A=other._rot.transpose()*(_trans-other._trans);
  dirX=rotB2A.row(0)*_ext.x();
  dirY=rotB2A.row(1)*_ext.y();
  dir=PT(std::abs(dirX.x())+std::abs(dirY.x()),std::abs(dirX.y())+std::abs(dirY.y()));
  if(SAT<T>::testSegment(PT(-other._ext.x(),other._ext.x()),PT(transB2A.x()-dir.x(),transB2A.x()+dir.x())))
    return false;
  if(SAT<T>::testSegment(PT(-other._ext.y(),other._ext.y()),PT(transB2A.y()-dir.y(),transB2A.y()+dir.y())))
    return false;
  return true;
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,2>::intersect(const BBox<T,2>& other) const
{
  return intersect(OBBTpl<T,2>(other));
}
template <typename T>
EIGEN_DEVICE_FUNC void OBBTpl<T,2>::writeVTK(VTKWriter<T>& os) const
{
  os.setRelativeIndex();
  std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
  PT pt1=_rot*PT(-_ext.x(),-_ext.y())+_trans;
  PT pt2=_rot*PT( _ext.x(),-_ext.y())+_trans;
  PT pt3=_rot*PT( _ext.x(), _ext.y())+_trans;
  PT pt4=_rot*PT(-_ext.x(), _ext.y())+_trans;
  vss.push_back(PT3(pt1.x(),pt1.y(),0.0f));
  vss.push_back(PT3(pt2.x(),pt2.y(),0.0f));
  vss.push_back(PT3(pt3.x(),pt3.y(),0.0f));
  vss.push_back(PT3(pt4.x(),pt4.y(),0.0f));
  os.appendPoints(vss.begin(),vss.end());

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
  iss.push_back(Vec3i(0,1,0));
  iss.push_back(Vec3i(1,2,0));
  iss.push_back(Vec3i(2,3,0));
  iss.push_back(Vec3i(3,0,0));
  os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
}
//OBBTpl<T,3>
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,3>::OBBTpl() {}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,3>::OBBTpl(const BBox<T,3>& bb)
{
  _ext=bb.getExtent()*(T)0.5f;
  _rot.setIdentity();
  _trans=bb._maxC-_ext;
}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,3>::OBBTpl(const MAT& rot,const PT& trans,const BBox<T,3>& bb)
{
  new (this)OBBTpl(bb);
  _rot=rot;
  _trans=_rot*_trans+trans;
}
template <typename T>
EIGEN_DEVICE_FUNC OBBTpl<T,3>::OBBTpl(const MAT& rot,const PT& trans,const PT& ext):_rot(rot),_trans(trans),_ext(ext) {}
template <typename T>
bool OBBTpl<T,3>::write(std::ostream& os) const
{
  writeBinaryData(_rot,os);
  writeBinaryData(_trans,os);
  writeBinaryData(_ext,os);
  return os.good();
}
template <typename T>
bool OBBTpl<T,3>::read(std::istream& is)
{
  readBinaryData(_rot,is);
  readBinaryData(_trans,is);
  readBinaryData(_ext,is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,3>::closest(PT pt,PT& n,PT* normal) const
{
  pt=(_rot.transpose()*(pt-_trans)).eval();
  bool inside=closestInner(pt,n,normal);
  n=(_rot*n).eval();
  if(normal)
    *normal=(_rot**normal).eval();
  return inside;
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,3>::closestInner(const PT& pt,PT& n,PT* normal) const
{
  BBox<T,3> box(-_ext,_ext);
  if(box.template containDim<3>(pt)) {
    T minDist=ScalarUtil<T>::scalar_max(),dist;
    for(sizeType d=0; d<3; d++) {
      dist=_ext[d]+pt[d];
      if(dist < minDist) {
        minDist=dist;
        n=-PT::Unit(d)*dist;
        if(normal)*normal=-PT::Unit(d);
      }
      dist=_ext[d]-pt[d];
      if(dist < minDist) {
        minDist=dist;
        n=PT::Unit(d)*dist;
        if(normal)*normal=PT::Unit(d);
      }
    }
    return true;
  } else {
    PT cp=box.closestTo(pt,3);
    n=cp-pt;
    if(normal)
      *normal=-n/max(n.norm(),(T)1E-6f);
    return false;
  }
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,3>::intersect(const OBBTpl<T,3>& other) const
{
  MAT rotB2A=_rot.transpose()*other._rot;
  PT transB2A=_rot.transpose()*(other._trans-_trans);
  PT dirs[3];

  //set two
  PT ctrAInB=other._rot.transpose()*(_trans-other._trans);
  dirs[0]=rotB2A.row(0)*_ext.x();
  dirs[1]=rotB2A.row(1)*_ext.y();
  dirs[2]=rotB2A.row(2)*_ext.z();
  PT dir=PT(std::abs(dirs[0].x())+std::abs(dirs[1].x())+std::abs(dirs[2].x()),
            std::abs(dirs[0].y())+std::abs(dirs[1].y())+std::abs(dirs[2].y()),
            std::abs(dirs[0].z())+std::abs(dirs[1].z())+std::abs(dirs[2].z()));
  if(SAT<T>::testSegment(PT2(-other._ext.x(),other._ext.x()),PT2(ctrAInB.x()-dir.x(),ctrAInB.x()+dir.x())))
    return false;
  if(SAT<T>::testSegment(PT2(-other._ext.y(),other._ext.y()),PT2(ctrAInB.y()-dir.y(),ctrAInB.y()+dir.y())))
    return false;
  if(SAT<T>::testSegment(PT2(-other._ext.z(),other._ext.z()),PT2(ctrAInB.z()-dir.z(),ctrAInB.z()+dir.z())))
    return false;

  //set two
  dirs[0]=rotB2A.col(0)*other._ext.x();
  dirs[1]=rotB2A.col(1)*other._ext.y();
  dirs[2]=rotB2A.col(2)*other._ext.z();
  dir=PT(std::abs(dirs[0].x())+std::abs(dirs[1].x())+std::abs(dirs[2].x()),
         std::abs(dirs[0].y())+std::abs(dirs[1].y())+std::abs(dirs[2].y()),
         std::abs(dirs[0].z())+std::abs(dirs[1].z())+std::abs(dirs[2].z()));
  if(SAT<T>::testSegment(PT2(-_ext.x(),_ext.x()),PT2(transB2A.x()-dir.x(),transB2A.x()+dir.x())))
    return false;
  if(SAT<T>::testSegment(PT2(-_ext.y(),_ext.y()),PT2(transB2A.y()-dir.y(),transB2A.y()+dir.y())))
    return false;
  if(SAT<T>::testSegment(PT2(-_ext.z(),_ext.z()),PT2(transB2A.z()-dir.z(),transB2A.z()+dir.z())))
    return false;

#define CROSS0(A) 0.0f,-A[2],A[1]
#define CROSS1(A) A[2],0.0f,-A[0]
#define CROSS2(A) -A[1],A[0],0.0f

#define ABSDOT0(A,B) std::abs(A[1])*B[1]+std::abs(A[2])*B[2]
#define ABSDOT1(A,B) std::abs(A[0])*B[0]+std::abs(A[2])*B[2]
#define ABSDOT2(A,B) std::abs(A[0])*B[0]+std::abs(A[1])*B[1]

#define ABSDOTB0(A,B) std::abs(A[1].dot(B))+std::abs(A[2].dot(B))
#define ABSDOTB1(A,B) std::abs(A[0].dot(B))+std::abs(A[2].dot(B))
#define ABSDOTB2(A,B) std::abs(A[0].dot(B))+std::abs(A[1].dot(B))

#define CROSS_TEST(AXIS)\
for(sizeType j=0;j<3;j++)\
{\
PT axis(CROSS##AXIS(rotB2A.col(j)));\
if(axis.squaredNorm() > ScalarUtil<T>::scalar_eps())\
{\
T extAPrj=ABSDOT##AXIS(axis,_ext);\
T ctrBPrj=transB2A.dot(axis);\
T extBPrj=ABSDOTB##AXIS(dirs,axis);\
if(SAT<T>::testSegment(PT2(-extAPrj,extAPrj),PT2(ctrBPrj-extBPrj,ctrBPrj+extBPrj)))\
return false;\
}\
}

#undef CROSS0
#undef CROSS1
#undef CROSS2

#undef ABSDOT0
#undef ABSDOT1
#undef ABSDOT2

#undef ABSDOTB0
#undef ABSDOTB1
#undef ABSDOTB2
  return true;
}
template <typename T>
EIGEN_DEVICE_FUNC bool OBBTpl<T,3>::intersect(const BBox<T,3>& other) const
{
  return intersect(OBBTpl<T,3>(other));
}
template <typename T>
EIGEN_DEVICE_FUNC void OBBTpl<T,3>::writeVTK(VTKWriter<T>& os) const
{
  os.setRelativeIndex();
  std::vector<PT,Eigen::aligned_allocator<PT> > vss;
  vss.push_back(_rot*PT(-_ext.x(),-_ext.y(),-_ext.z())+_trans);
  vss.push_back(_rot*PT( _ext.x(),-_ext.y(),-_ext.z())+_trans);
  vss.push_back(_rot*PT( _ext.x(), _ext.y(),-_ext.z())+_trans);
  vss.push_back(_rot*PT(-_ext.x(), _ext.y(),-_ext.z())+_trans);
  vss.push_back(_rot*PT(-_ext.x(),-_ext.y(), _ext.z())+_trans);
  vss.push_back(_rot*PT( _ext.x(),-_ext.y(), _ext.z())+_trans);
  vss.push_back(_rot*PT( _ext.x(), _ext.y(), _ext.z())+_trans);
  vss.push_back(_rot*PT(-_ext.x(), _ext.y(), _ext.z())+_trans);
  os.appendPoints(vss.begin(),vss.end());

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > iss;
  iss.push_back(Vec3i(0,1,0));
  iss.push_back(Vec3i(1,2,0));
  iss.push_back(Vec3i(2,3,0));
  iss.push_back(Vec3i(3,0,0));
  iss.push_back(Vec3i(4,5,0));
  iss.push_back(Vec3i(5,6,0));
  iss.push_back(Vec3i(6,7,0));
  iss.push_back(Vec3i(7,4,0));
  iss.push_back(Vec3i(0,4,0));
  iss.push_back(Vec3i(1,5,0));
  iss.push_back(Vec3i(2,6,0));
  iss.push_back(Vec3i(3,7,0));
  os.appendCells(iss.begin(),iss.end(),VTKWriter<T>::LINE,true);
}
//KDOP18
template <typename T>
EIGEN_DEVICE_FUNC KDOP18<T>::KDOP18()
{
  empty();
}
template <typename T>
EIGEN_DEVICE_FUNC KDOP18<T>::KDOP18(const PT& v)
{
  _dist[0] = _dist[9]= v[0];
  _dist[1] = _dist[10] = v[1];
  _dist[2] = _dist[11] = v[2];

  T d3, d4, d5, d6, d7, d8;
  getDistances(v, d3, d4, d5, d6, d7, d8);
  _dist[3] = _dist[12] = d3;
  _dist[4] = _dist[13] = d4;
  _dist[5] = _dist[14] = d5;
  _dist[6] = _dist[15] = d6;
  _dist[7] = _dist[16] = d7;
  _dist[8] = _dist[17] = d8;
}
template <typename T>
EIGEN_DEVICE_FUNC KDOP18<T>::KDOP18(const PT& a,const PT& b)
{
  _dist[0]= min(a[0], b[0]);
  _dist[9]= max(a[0], b[0]);
  _dist[1]= min(a[1], b[1]);
  _dist[10] = max(a[1], b[1]);
  _dist[2]= min(a[2], b[2]);
  _dist[11] = max(a[2], b[2]);

  T ad3, ad4, ad5, ad6, ad7, ad8;
  getDistances(a, ad3, ad4, ad5, ad6, ad7, ad8);
  T bd3, bd4, bd5, bd6, bd7, bd8;
  getDistances(b, bd3, bd4, bd5, bd6, bd7, bd8);
  _dist[3]= min(ad3, bd3);
  _dist[12] = max(ad3, bd3);
  _dist[4]= min(ad4, bd4);
  _dist[13] = max(ad4, bd4);
  _dist[5]= min(ad5, bd5);
  _dist[14] = max(ad5, bd5);
  _dist[6]= min(ad6, bd6);
  _dist[15] = max(ad6, bd6);
  _dist[7]= min(ad7, bd7);
  _dist[16] = max(ad7, bd7);
  _dist[8]= min(ad8, bd8);
  _dist[17] = max(ad8, bd8);
}
template <typename T>
bool KDOP18<T>::write(std::ostream& os) const
{
  for(sizeType i=0; i<18; i++)
    writeBinaryData(_dist[i],os);
  return os.good();
}
template <typename T>
bool KDOP18<T>::read(std::istream& is)
{
  for(sizeType i=0; i<18; i++)
    readBinaryData(_dist[i],is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::reset()
{
  empty();
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::empty()
{
  for (int i=0; i<9; i++) {
    _dist[i] = ScalarUtil<T>::scalar_max();
    _dist[i+9] = -ScalarUtil<T>::scalar_max();
  }
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::enlarged(T len)
{
  for(int i=0; i<3; i++) {
    _dist[i]-=len;
    _dist[i+9]+=len;
  }
  for(int i=3; i<9; i++) {
    _dist[i]-=len*2.0f;
    _dist[i+9]+=len*2.0f;
  }
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::enlarged(T len,sizeType dim)
{
  ASSERT(dim == 3)
  enlarged(len);
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::setPoints(const PT& a,const PT& b,const PT& c)
{
  new(this) KDOP18(a,b);
  setUnion(c);
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::setUnion(const PT& p)
{
  _dist[0]= min(p[0], _dist[0]);
  _dist[9]= max(p[0], _dist[9]);
  _dist[1]= min(p[1], _dist[1]);
  _dist[10] = max(p[1], _dist[10]);
  _dist[2]= min(p[2], _dist[2]);
  _dist[11] = max(p[2], _dist[11]);

  T d3, d4, d5, d6, d7, d8;
  getDistances(p, d3, d4, d5, d6, d7, d8);
  _dist[3]= min(d3, _dist[3]);
  _dist[12] = max(d3, _dist[12]);
  _dist[4]= min(d4, _dist[4]);
  _dist[13] = max(d4, _dist[13]);
  _dist[5]= min(d5, _dist[5]);
  _dist[14] = max(d5, _dist[14]);
  _dist[6]= min(d6, _dist[6]);
  _dist[15] = max(d6, _dist[15]);
  _dist[7]= min(d7, _dist[7]);
  _dist[16] = max(d7, _dist[16]);
  _dist[8]= min(d8, _dist[8]);
  _dist[17] = max(d8, _dist[17]);
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::setUnion(const KDOP18& b)
{
  _dist[0]= min(b._dist[0], _dist[0]);
  _dist[9]= max(b._dist[9], _dist[9]);
  _dist[1]= min(b._dist[1], _dist[1]);
  _dist[10] = max(b._dist[10], _dist[10]);
  _dist[2]= min(b._dist[2], _dist[2]);
  _dist[11] = max(b._dist[11], _dist[11]);
  _dist[3]= min(b._dist[3], _dist[3]);
  _dist[12] = max(b._dist[12], _dist[12]);
  _dist[4]= min(b._dist[4], _dist[4]);
  _dist[13] = max(b._dist[13], _dist[13]);
  _dist[5]= min(b._dist[5], _dist[5]);
  _dist[14] = max(b._dist[14], _dist[14]);
  _dist[6]= min(b._dist[6], _dist[6]);
  _dist[15] = max(b._dist[15], _dist[15]);
  _dist[7]= min(b._dist[7], _dist[7]);
  _dist[16] = max(b._dist[16], _dist[16]);
  _dist[8]= min(b._dist[8], _dist[8]);
  _dist[17] = max(b._dist[17], _dist[17]);
}
template <typename T>
EIGEN_DEVICE_FUNC KDOP18<T> KDOP18<T>::getUnion(const KDOP18& b) const
{
  KDOP18<T> ret=*this;
  ret.setUnion(b);
  return ret;
}
template <typename T>
EIGEN_DEVICE_FUNC typename KDOP18<T>::PT KDOP18<T>::minCorner() const
{
  return PT(_dist[0],_dist[1 ],_dist[2 ]);
}
template <typename T>
EIGEN_DEVICE_FUNC typename KDOP18<T>::PT KDOP18<T>::maxCorner() const
{
  return PT(_dist[9],_dist[10],_dist[11]);
}
template <typename T>
EIGEN_DEVICE_FUNC bool KDOP18<T>::intersect(const KDOP18& b,sizeType DIM) const
{
  ASSERT(DIM == 3)
  for (int i=0; i<9; i++) {
    if (_dist[i] > b._dist[i+9]) return false;
    if (_dist[i+9] < b._dist[i]) return false;
  }
  return true;
}
template <typename T>
EIGEN_DEVICE_FUNC bool KDOP18<T>::intersect(const BBox<T>& b,sizeType DIM) const
{
  ASSERT(DIM == 3)
  BBox<T> bb(minCorner(),maxCorner());
  return bb.intersect(b);
}
template <typename T>
EIGEN_DEVICE_FUNC bool KDOP18<T>::contain(const PT& p) const
{
  for (int i=0; i<3; i++) {
    if (p[i] < _dist[i] || p[i] > _dist[i+9])
      return false;
  }

  T d[6];
  getDistances(p, d);
  for (int i=3; i<9; i++) {
    if (d[i-3] < _dist[i] || d[i-3] > _dist[i+9])
      return false;
  }

  return true;
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::getDistances(const PT& p,T &d3, T &d4, T &d5, T &d6, T &d7, T &d8)
{
  d3 = p[0] + p[1];
  d4 = p[0] + p[2];
  d5 = p[1] + p[2];
  d6 = p[0] - p[1];
  d7 = p[0] - p[2];
  d8 = p[1] - p[2];
}
template <typename T>
EIGEN_DEVICE_FUNC void KDOP18<T>::getDistances(const PT& p, T d[])
{
  d[0] = p[0] + p[1];
  d[1] = p[0] + p[2];
  d[2] = p[1] + p[2];
  d[3] = p[0] - p[1];
  d[4] = p[0] - p[2];
  d[5] = p[1] - p[2];
}
template <typename T>
EIGEN_DEVICE_FUNC T KDOP18<T>::getDistances(const PT &p, int i)
{
  if (i == 0) return p[0]+p[1];
  if (i == 1) return p[0]+p[2];
  if (i == 2) return p[1] + p[2];
  if (i == 3) return p[0] - p[1];
  if (i == 4) return p[0] - p[2];
  if (i == 5) return p[1] - p[2];
  return 0;
}
//Sphere
template <typename T>
EIGEN_DEVICE_FUNC Sphere<T>::Sphere() {}
template <typename T>
EIGEN_DEVICE_FUNC Sphere<T>::Sphere(const PT& ctr,const T& rad):_ctr(ctr),_rad(rad) {}
template <typename T>
bool Sphere<T>::write(std::ostream& os) const
{
  writeBinaryData(_ctr,os);
  writeBinaryData(_rad,os);
  return os.good();
}
template <typename T>
bool Sphere<T>::read(std::istream& is)
{
  readBinaryData(_ctr,is);
  readBinaryData(_rad,is);
  return is.good();
}
template <typename T>
EIGEN_DEVICE_FUNC bool Sphere<T>::closest(const PT& pt,PT& n,PT* normal) const
{
  PT dir=pt-_ctr;
  T len=dir.norm();
  n=dir*(_rad-len)/max(len,(T)1E-6f);
  if(normal)*normal=dir/max(len,(T)1E-6f);
  return len < _rad;
}
template <typename T>
EIGEN_DEVICE_FUNC bool Sphere<T>::intersect(const PT& a,const PT& b) const
{
  T A=(b-a).squaredNorm();
  if(A < ScalarUtil<T>::scalar_eps())
    return (_ctr-a).norm() < _rad;

  T B=2.0f*(b-a).dot(a-_ctr);
  T C=(a-_ctr).squaredNorm()-_rad*_rad;
  T delta=B*B-4.0f*A*C;
  if(delta < 0.0f)
    return false;

  T L=(-B-std::sqrt(delta))/(2*A);
  T H=(-B+std::sqrt(delta))/(2*A);
  return L < 1 && H > 0;
}

#define DEF_BBOX(T,DIM) \
template struct BBox<T,DIM>; \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<sizeType>(const BBox<sizeType,DIM>& other);  \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<unsigned char>(const BBox<unsigned char,DIM>& other);  \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<scalarF>(const BBox<scalarF,DIM>& other);  \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<scalarD>(const BBox<scalarD,DIM>& other);

DEF_BBOX(scalarD,2)
DEF_BBOX(scalarD,3)
template void TriangleTpl<scalarD>::calcPointDist(const Vec3d& pt,scalarD& sqrDistance,Vec3d& cp,Vec3d& b) const;
template void TriangleTpl<scalarD>::calcPointDist(const Vec3d& pt,scalarD& sqrDistance,Vec3d& cp,Vec4d& b) const;
template class LineSegTpl<scalarD>;
template class PlaneTpl<scalarD>;
template class TriangleTpl<scalarD>;
template class TetrahedronTpl<scalarD>;
template class OBBTpl<scalarD,2>;
template class OBBTpl<scalarD,3>;
template class KDOP18<scalarD>;
template class Sphere<scalarD>;

DEF_BBOX(scalarF,2)
DEF_BBOX(scalarF,3)
template void TriangleTpl<scalarF>::calcPointDist(const Vec3f& pt,scalarF& sqrDistance,Vec3f& cp,Vec3f& b) const;
template void TriangleTpl<scalarF>::calcPointDist(const Vec3f& pt,scalarF& sqrDistance,Vec3f& cp,Vec4f& b) const;
template class LineSegTpl<scalarF>;
template class PlaneTpl<scalarF>;
template class TriangleTpl<scalarF>;
template class TetrahedronTpl<scalarF>;
template class OBBTpl<scalarF,2>;
template class OBBTpl<scalarF,3>;
template class KDOP18<scalarF>;
template class Sphere<scalarF>;

DEF_BBOX(sizeType,2)
DEF_BBOX(sizeType,3)
template void TriangleTpl<sizeType>::calcPointDist(const Vec3i& pt,sizeType& sqrDistance,Vec3i& cp,Vec3i& b) const;
template void TriangleTpl<sizeType>::calcPointDist(const Vec3i& pt,sizeType& sqrDistance,Vec3i& cp,Vec4i& b) const;
template class LineSegTpl<sizeType>;
template class PlaneTpl<sizeType>;
template class TriangleTpl<sizeType>;
template class TetrahedronTpl<sizeType>;
template class OBBTpl<sizeType,2>;
template class OBBTpl<sizeType,3>;
template class KDOP18<sizeType>;
template class Sphere<sizeType>;

DEF_BBOX(char,2)
DEF_BBOX(char,3)
template void TriangleTpl<char>::calcPointDist(const Vec3c& pt,char& sqrDistance,Vec3c& cp,Vec3c& b) const;
template void TriangleTpl<char>::calcPointDist(const Vec3c& pt,char& sqrDistance,Vec3c& cp,Vec4c& b) const;
template class LineSegTpl<char>;
template class PlaneTpl<char>;
template class TriangleTpl<char>;
template class TetrahedronTpl<char>;
template class OBBTpl<char,2>;
template class OBBTpl<char,3>;
template class KDOP18<char>;
template class Sphere<char>;

DEF_BBOX(unsigned char,2)
DEF_BBOX(unsigned char,3)
template void TriangleTpl<unsigned char>::calcPointDist(const Vec3uc& pt,unsigned char& sqrDistance,Vec3uc& cp,Vec3uc& b) const;
template void TriangleTpl<unsigned char>::calcPointDist(const Vec3uc& pt,unsigned char& sqrDistance,Vec3uc& cp,Vec4uc& b) const;
template class LineSegTpl<unsigned char>;
template class PlaneTpl<unsigned char>;
template class TriangleTpl<unsigned char>;
template class TetrahedronTpl<unsigned char>;
template class OBBTpl<unsigned char,2>;
template class OBBTpl<unsigned char,3>;
template class KDOP18<unsigned char>;
template class Sphere<unsigned char>;

PRJ_END
