#ifndef QCQP_SOLVER_MOSEK_H
#define QCQP_SOLVER_MOSEK_H

#include "QCQPSolver.h"

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct QCQPSolverMosek : public QCQPSolver<T>
{
  using typename QCQPSolver<T>::QCQP_RETURN_CODE;
  DECL_MAP_TYPES_T
  QCQPSolverMosek();
  virtual ~QCQPSolverMosek();
  //consistent QP
  template <typename MAT>
  QCQP_RETURN_CODE solveQPTpl(Vec& x,const MAT& H,const Vec& g,const MAT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false);
  QCQP_RETURN_CODE solveQP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
  QCQP_RETURN_CODE solveQP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
  //inconsistent QP
  template <typename MAT>
  QCQP_RETURN_CODE solveL1QPTpl(Vec& x,const MAT& H,const Vec& g,const MAT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false);
  QCQP_RETURN_CODE solveL1QP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
  QCQP_RETURN_CODE solveL1QP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
private:
  //helper
  QCQP_RETURN_CODE putHess(void* task,const MatT& H,bool isA,bool callback) const;
  QCQP_RETURN_CODE putHess(void* task,const SMat& H,bool isA,bool callback) const;
  void* _env;
};

PRJ_END

#endif
