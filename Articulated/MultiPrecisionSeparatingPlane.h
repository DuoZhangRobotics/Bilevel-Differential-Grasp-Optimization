#ifndef MULTI_PRECISION_SEPARATING_PLANE_H
#define MULTI_PRECISION_SEPARATING_PLANE_H

#include "MultiPrecisionLQP.h"

PRJ_BEGIN

template <typename T>
class MultiPrecisionSeparatingPlane : public MultiPrecisionLQP<T>
{
public:
  DECL_MAP_TYPES_T
  typedef std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> PSS;
  MultiPrecisionSeparatingPlane(Options& ops,Vec4T& plane,T d0);
  void clearPoints();
  void resetPoints(const PSS& pss);
  void resetPoints(const Mat3XT& pssL,const Mat3XT& pssR);
  virtual T computeFGH(T mu,const Vec& w,Vec* g=NULL,DMat* h=NULL) const override;
  virtual Vec sampleValidW() const override;
  Vec4T solve(bool& succ);
  static T clog(T dist,T* D,T* DD,T d0,T coef);
protected:
  virtual void initialGuess(Vec& w) const override;
  virtual void limitAlpha(T& alpha,const Vec& w,const Vec& dx) const override;
  PSS _pssL,_pssR;
  Vec4T& _plane;
  T _d0;
};

PRJ_END

#endif
