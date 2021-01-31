#include "DSSQPObjective.h"
#include "QCQPSolverQPOASES.h"

USE_PRJ_NAMESPACE

template <typename T>
QCQPSolverQPOASES<T>::QCQPSolverQPOASES() {}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverQPOASES<T>::solveQP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  ASSERT_MSG(QCones.empty(),"QPOASES does not support second order cone with QP interface")
  qpOASES::int_t nWSR=10000;
  qpOASES::returnValue ret;
  Eigen::Matrix<scalarD,-1,-1,Eigen::RowMajor> Hd,cjacd;
  Eigen::Matrix<scalarD,-1,1> gd,lbd,ubd,lbAd,ubAd;
  Hd=H.unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  gd=g.unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(cjac)
    cjacd=cjac->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(lb)
    lbd=lb->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(ub)
    ubd=ub->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(lbA)
    lbAd=lbA->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(ubA)
    ubAd=ubA->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(!_prob.isInitialised() || _prob.getNV()!=x.size() || _prob.getNC()!=(cjac?cjac->rows():0)) {
    _prob=qpOASES::SQProblem(x.size(),cjac?cjac->rows():0);
    ret=_prob.init(Hd.data(),
                   gd.data(),
                   cjac?cjacd.data():NULL,
                   lb?lbd.data():NULL,
                   ub?ubd.data():NULL,
                   lbA?lbAd.data():NULL,
                   ubA?ubAd.data():NULL,
                   nWSR);
  } else {
    ret=_prob.hotstart(Hd.data(),
                       gd.data(),
                       cjac?cjacd.data():NULL,
                       lb?lbd.data():NULL,
                       ub?ubd.data():NULL,
                       lbA?lbAd.data():NULL,
                       ubA?ubAd.data():NULL,
                       nWSR);
  }
  if(ret==qpOASES::SUCCESSFUL_RETURN && _prob.isSolved()) {
    Cold xd=Cold::Zero(g.size());
    _prob.getPrimalSolution(xd.data());
    x=xd.template cast<T>();
    return QCQPSolver<T>::SOLVED;
  } else if(_prob.isInfeasible())
    return QCQPSolver<T>::INFEASIBLE;
  else if(_prob.isUnbounded())
    return QCQPSolver<T>::UNBOUNDED;
  else if(_prob.getHessianType()==qpOASES::HessianType::HST_INDEF)
    return QCQPSolver<T>::NOT_POSITIVE_DEFINITE;
  else return QCQPSolver<T>::UNKNOWN;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverQPOASES<T>::solveQP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  ASSERT_MSG(QCones.empty(),"QPOASES does not support second order cone with QP interface")
  qpOASES::int_t nWSR=10000;
  qpOASES::returnValue ret;
  qpOASES::SymSparseMat HQP,cjacQP;
  Eigen::SparseMatrix<scalarD,0,qpOASES::sparse_int_t> Hd,cjacd;
  Eigen::Matrix<scalarD,-1,1> gd,lbd,ubd,lbAd,ubAd;
  {
    Hd=H.unaryExpr([&](T in) {
      return (scalarD)std::to_double(in);
    });
    HQP=qpOASES::SymSparseMat(Hd.rows(),Hd.cols(),Hd.innerIndexPtr(),Hd.outerIndexPtr(),Hd.valuePtr());
    HQP.createDiagInfo();
  }
  gd=g.unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(cjac) {
    cjacd=cjac->unaryExpr([&](T in) {
      return (scalarD)std::to_double(in);
    });
    cjacQP=qpOASES::SymSparseMat(cjacd.rows(),cjacd.cols(),cjacd.innerIndexPtr(),cjacd.outerIndexPtr(),cjacd.valuePtr());
  }
  if(lb)
    lbd=lb->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(ub)
    ubd=ub->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(lbA)
    lbAd=lbA->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(ubA)
    ubAd=ubA->unaryExpr([&](T in) {
    return (scalarD)std::to_double(in);
  });
  if(!_prob.isInitialised() || _prob.getNV()!=x.size() || _prob.getNC()!=(cjac?cjac->rows():0)) {
    _prob=qpOASES::SQProblem(x.size(),cjac?cjac->rows():0);
    ret=_prob.init(&HQP,
                   gd.data(),
                   cjac?&cjacQP:NULL,
                   lb?lbd.data():NULL,
                   ub?ubd.data():NULL,
                   lbA?lbAd.data():NULL,
                   ubA?ubAd.data():NULL,
                   nWSR);
  } else {
    ret=_prob.hotstart(&HQP,
                       gd.data(),
                       cjac?&cjacQP:NULL,
                       lb?lbd.data():NULL,
                       ub?ubd.data():NULL,
                       lbA?lbAd.data():NULL,
                       ubA?ubAd.data():NULL,
                       nWSR);
  }
  if(ret==qpOASES::SUCCESSFUL_RETURN && _prob.isSolved()) {
    Cold xd=Cold::Zero(g.size());
    _prob.getPrimalSolution(xd.data());
    x=xd.template cast<T>();
    return QCQPSolver<T>::SOLVED;
  } else if(_prob.isInfeasible())
    return QCQPSolver<T>::INFEASIBLE;
  else if(_prob.isUnbounded())
    return QCQPSolver<T>::UNBOUNDED;
  else if(_prob.getHessianType()==qpOASES::HessianType::HST_INDEF)
    return QCQPSolver<T>::NOT_POSITIVE_DEFINITE;
  else return QCQPSolver<T>::UNKNOWN;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverQPOASES<T>::solveL1QP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverQPOASES<T>::solveL1QP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
//instance
PRJ_BEGIN
template struct QCQPSolverQPOASES<double>;
#ifdef ALL_TYPES
template struct QCQPSolverQPOASES<__float128>;
template struct QCQPSolverQPOASES<mpfr::mpreal>;
#endif
PRJ_END
