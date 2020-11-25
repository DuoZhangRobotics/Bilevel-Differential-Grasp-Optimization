#ifndef IMPLICIT_FUNC_INTERFACE_H
#define IMPLICIT_FUNC_INTERFACE_H

#include "MathBasic.h"
#include "CollisionDetection.h"

PRJ_BEGIN

template <typename T,typename TI,typename TG>
struct Grid;

template <typename T>
class ImplicitFunc
{
public:
  template <typename T2>
  struct ScalarType
  {
    typedef T2 Type;
  };
  template <typename T2,int R,int C>
  struct ScalarType<Eigen::Matrix<T2,R,C> >
  {
    typedef T2 Type;
  };
  typedef typename ScalarType<T>::Type ST;
  typedef typename ScalarUtil<ST>::ScalarVec3 PT;
  virtual ~ImplicitFunc() {}
  virtual void reset() {}
  virtual void beginSampleSet(const Grid<T,ST,std::vector<T,Eigen::aligned_allocator<T> > >& toBeSampled) {}
  virtual void endSampleSet(const Grid<T,ST,std::vector<T,Eigen::aligned_allocator<T> > >& toBeSampled) {}
  virtual T operator()(const PT& pos) const =0;
  virtual PT grad(const PT& pos) const
  {
    FUNCTION_NOT_IMPLEMENTED
	return PT();
  }
  virtual PT gradN(const PT& pos,T delta) const
  {
    PT val;
    for(sizeType d=0; d<3; d++)
      val[d]=(operator()(pos+PT::Unit(d)*delta)-operator()(pos-PT::Unit(d)*delta))/delta/2;
    return val;
  }
  virtual BBox<ST> getBB() const {
    return BBox<ST>();
  }
};

template <typename T>
class VelFunc
{
public:
  typedef typename ImplicitFunc<T>::template ScalarType<T>::Type ST;
  virtual ~VelFunc() {}
  virtual typename ScalarUtil<ST>::ScalarVec3 operator()(const typename ScalarUtil<ST>::ScalarVec3& pos) const =0;
};

PRJ_END

#endif
