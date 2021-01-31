#ifndef QCQP_SOLVER_GUROBI_H
#define QCQP_SOLVER_GUROBI_H

#include "QCQPSolver.h"

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct QCQPSolverGurobi : public QCQPSolver<T>
{
  using typename QCQPSolver<T>::QCQP_RETURN_CODE;
  DECL_MAP_TYPES_T
  QCQPSolverGurobi();
  virtual ~QCQPSolverGurobi();
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
  QCQP_RETURN_CODE putHess(void* task,const MatT& H,bool callback) const;
  QCQP_RETURN_CODE putHess(void* task,const SMat& H,bool callback) const;
  QCQP_RETURN_CODE putCons(void* m,const MatT& cjac,const Vec* lbA,const Vec* ubA,std::vector<Eigen::Triplet<double,int>>& trips,int nrVar,bool callback) const;
  QCQP_RETURN_CODE putCons(void* m,const SMat& cjac,const Vec* lbA,const Vec* ubA,std::vector<Eigen::Triplet<double,int>>& trips,int nrVar,bool callback) const;
  QCQP_RETURN_CODE putQCone(void* m,const Vec& x,const Vec& g,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,int nrVar,bool callback) const;
  void* _env;
};

PRJ_END

#endif
