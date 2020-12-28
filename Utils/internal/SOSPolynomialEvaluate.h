template <typename T>
class ParallelEvaluate
{
public:
  typedef typename Eigen::Matrix<T,-1,1> VEC;
  typedef typename Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,sizeType>> STRIPS;
  typedef Eigen::Triplet<T,sizeType> STRIP;
  //utility
  static void gradient(VEC& g,sizeType r,const T& v)
  {
    OMP_CRITICAL_
    {
      T& val=g[r];
      val+=v;
    }
  }
  static void hessian(MAT& m,sizeType r,sizeType c,const T& v)
  {
    OMP_CRITICAL_
    {
      T& val=m(r,c);
      val+=v;
    }
  }
  static void hessianSym(MAT& m,sizeType r,sizeType c,const T& v)
  {
    OMP_CRITICAL_
    {
      T& val=m(r,c);
      val+=v;
    }
    OMP_CRITICAL_
    {
      T& val=m(c,r);
      val+=v;
    }
  }
  static void hessian(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
  }
  static void hessianSym(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
    m.push_back(STRIP(c,r,v));
  }
  template <typename POLY,typename HESS>
  static T eval(POLY& p,const COLD& x,VEC* grad,HESS* hess) {
    T ret=0;
    for(sizeType i=0; i<(sizeType)p._terms.size(); i++)
      ret+=p._terms[i].eval(x,grad,hess);
    return ret;
  }
};
template <>
class ParallelEvaluate<scalarD>
{
public:
  typedef scalarD T;
  typedef typename Eigen::Matrix<T,-1,1> VEC;
  typedef typename Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,sizeType>> STRIPS;
  typedef Eigen::Triplet<T,sizeType> STRIP;
  //utility
  static void gradient(VEC& g,sizeType r,const T& v)
  {
    T& val=g[r];
    OMP_ATOMIC_
    val+=v;
  }
  static void hessian(MAT& m,sizeType r,sizeType c,const T& v)
  {
    T& val=m(r,c);
    OMP_ATOMIC_
    val+=v;
  }
  static void hessianSym(MAT& m,sizeType r,sizeType c,const T& v)
  {
    {
      T& val=m(r,c);
      OMP_ATOMIC_
      val+=v;
    }
    {
      T& val=m(c,r);
      OMP_ATOMIC_
      val+=v;
    }
  }
  static void hessian(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
  }
  static void hessianSym(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
    m.push_back(STRIP(c,r,v));
  }
  template <typename POLY,typename HESS>
  static T eval(POLY& p,const COLD& x,VEC* grad,HESS* hess) {
    T ret=0;
    OMP_PARALLEL_FOR_I(OMP_ADD(ret))
    for(sizeType i=0; i<(sizeType)p._terms.size(); i++)
      ret+=p._terms[i].eval(x,grad,hess);
    return ret;
  }
};
template <typename TT>
class ParallelEvaluate<COMMON::SMatNumber<TT>>
{
public:
  typedef typename COMMON::SMatNumber<TT> T;
  typedef typename Eigen::Matrix<T,-1,1> VEC;
  typedef typename Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,sizeType>> STRIPS;
  typedef Eigen::Triplet<T,sizeType> STRIP;
  //utility
  static void gradient(VEC& g,sizeType r,const T& v)
  {
    T& val=g[r];
    val+=v;
  }
  static void hessian(MAT& m,sizeType r,sizeType c,const T& v)
  {
    T& val=m(r,c);
    val+=v;
  }
  static void hessianSym(MAT& m,sizeType r,sizeType c,const T& v)
  {
    {
      T& val=m(r,c);
      val+=v;
    }
    {
      T& val=m(c,r);
      val+=v;
    }
  }
  static void hessian(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
  }
  static void hessianSym(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
    m.push_back(STRIP(c,r,v));
  }
  template <typename POLY,typename HESS>
  static T eval(POLY& p,const COLD& x,VEC* grad,HESS* hess) {
    T ret=0;
    for(sizeType i=0; i<(sizeType)p._terms.size(); i++)
      ret+=p._terms[i].eval(x,grad,hess);
    return ret;
  }
};
template <typename T2,char LABEL>
class ParallelEvaluate<SOSPolynomial<T2,LABEL>>
{
public:
  typedef SOSPolynomial<T2,LABEL> T;
  typedef typename Eigen::Matrix<T,-1,1> VEC;
  typedef typename Eigen::Matrix<T,-1,-1> MAT;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,1> COLD;
  typedef Eigen::Matrix<typename ScalarOfT<T>::Type,-1,-1> MATD;
  typedef ParallelVector<Eigen::Triplet<T,sizeType>> STRIPS;
  typedef Eigen::Triplet<T,sizeType> STRIP;
  //utility
  static void gradient(VEC& g,sizeType r,const T& v)
  {
    g[r]+=v;
  }
  static void hessian(MAT& m,sizeType r,sizeType c,const T& v)
  {
    m(r,c)+=v;
  }
  static void hessianSym(MAT& m,sizeType r,sizeType c,const T& v)
  {
    m(r,c)+=v;
    m(c,r)+=v;
  }
  static void hessian(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
  }
  static void hessianSym(STRIPS& m,sizeType r,sizeType c,const T& v)
  {
    m.push_back(STRIP(r,c,v));
    m.push_back(STRIP(c,r,v));
  }
  template <typename POLY,typename HESS>
  static T eval(POLY& p,const COLD& x,VEC* grad,HESS* hess) {
    T ret=ScalarOfT<T>::convert(0);
    for(sizeType i=0; i<(sizeType)p._terms.size(); i++)
      ret+=p._terms[i].eval(x,grad,hess);
    return ret;
  }
};
