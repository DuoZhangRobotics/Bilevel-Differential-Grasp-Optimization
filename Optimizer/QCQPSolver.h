#ifndef QCQP_SOLVER_H
#define QCQP_SOLVER_H

#include <Utils/Scalar.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct QCQPSolver
{
  DECL_MAP_TYPES_T
  enum QCQP_RETURN_CODE
  {
    SOLVED,
    NOT_POSITIVE_DEFINITE,
    INFEASIBLE,
    UNBOUNDED,
    UNKNOWN,
    API_OK,
  };
  virtual ~QCQPSolver();
  bool solveQP(Vec& x,MatT H,const Vec& g,const MatT& cjac,const Vec& w,T beta);
  //consistent QP
  virtual QCQP_RETURN_CODE solveQP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false)=0;
  virtual QCQP_RETURN_CODE solveQP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false)=0;
  //inconsistent QP
  virtual QCQP_RETURN_CODE solveL1QP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false)=0;
  virtual QCQP_RETURN_CODE solveL1QP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false)=0;
};

PRJ_END

#endif
