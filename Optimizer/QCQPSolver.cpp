#include "QCQPSolver.h"

USE_PRJ_NAMESPACE

template <typename T>
QCQPSolver<T>::~QCQPSolver() {}
template <typename T>
bool QCQPSolver<T>::solveQP(Vec& x,MatT H,const Vec& g,const MatT& cjac,const Vec& w,T beta)
{
  H.diagonal().array()+=beta;
  Vec lb=-w,ubA=Vec::Ones(g.size())-cjac*w;
  return solveQP(x,H,g,&cjac,&lb,NULL,NULL,&ubA,std::vector<Coli,Eigen::aligned_allocator<Coli>>());
}
//instance
PRJ_BEGIN
template struct QCQPSolver<double>;
#ifdef ALL_TYPES
template struct QCQPSolver<__float128>;
template struct QCQPSolver<mpfr::mpreal>;
#endif
PRJ_END
