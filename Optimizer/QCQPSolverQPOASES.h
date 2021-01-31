#ifndef QCQP_SOLVER_QPOASES_H
#define QCQP_SOLVER_QPOASES_H

#include "QCQPSolver.h"
#include <qpOASES.hpp>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct QCQPSolverQPOASES : public QCQPSolver<T>
{
  using typename QCQPSolver<T>::QCQP_RETURN_CODE;
  DECL_MAP_TYPES_T
  QCQPSolverQPOASES();
  //consistent QP
  QCQP_RETURN_CODE solveQP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
  QCQP_RETURN_CODE solveQP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
  //inconsistent QP
  QCQP_RETURN_CODE solveL1QP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
  QCQP_RETURN_CODE solveL1QP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback=false) override;
private:
  qpOASES::SQProblem _prob;
};

PRJ_END

#endif
