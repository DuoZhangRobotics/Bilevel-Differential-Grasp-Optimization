#include "DSSQPObjective.h"
#include "QCQPSolverMosek.h"

USE_PRJ_NAMESPACE

#ifdef USE_MOSEK_9
#include <mosek.h>
#define MOSEK_SAFE_CALL(FUNC)\
if((trmcode=FUNC)!=MSK_RES_OK) {\
  if(callback) {\
    WARNINGV("Mosek failed (code=%d) while calling: %s",trmcode,STRINGIFY_OMP(FUNC))\
  }\
  if(task)\
    MSK_deletetask(&task);\
  if(trmcode==MSK_RES_ERR_OBJ_Q_NOT_PSD)\
    return QCQPSolver<T>::NOT_POSITIVE_DEFINITE;\
  else return QCQPSolver<T>::UNKNOWN;\
}
static void MSKAPI printstr(void*,const char str[]) {
  printf("%s", str);
}
template <typename T>
QCQPSolverMosek<T>::QCQPSolverMosek():_env(NULL) {
  MSKrescodee trmcode;
  if((trmcode=MSK_makeenv(&_env,NULL))!=MSK_RES_OK) {
    ASSERT_MSGV(false,"Mosek failed (code=%d) while calling: %s",trmcode,STRINGIFY_OMP(MSK_makeenv))
  }
}
template <typename T>
QCQPSolverMosek<T>::~QCQPSolverMosek() {
  if(_env)
    MSK_deleteenv(&_env);
}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveQPTpl(Vec& x,const MAT& H,const Vec& g,const MAT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  ASSERT_MSG(QCones.empty(),"Mosek does not support second order cone with QP interface")
  QCQP_RETURN_CODE ret=QCQPSolver<T>::UNKNOWN;
  MSKtask_t task=NULL;
  MSKrescodee trmcode;
  MSKsolstae solsta;
  MOSEK_SAFE_CALL(MSK_maketask(_env,cjac?cjac->rows():0,g.size(),&task))
  MOSEK_SAFE_CALL(MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,callback?printstr:NULL))
  MOSEK_SAFE_CALL(MSK_appendvars(task,g.size()))
  //constraint
  if(cjac) {
    MOSEK_SAFE_CALL(MSK_appendcons(task,cjac->rows()))
    for(sizeType i=0; i<cjac->rows(); i++) {
      MOSEK_SAFE_CALL(MSK_putconbound(task,i,
                                      lbA&&ubA?MSK_BK_RA:lbA?MSK_BK_LO:ubA?MSK_BK_UP:MSK_BK_FR,
                                      std::to_double(lbA?lbA->coeff(i):0),
                                      std::to_double(ubA?ubA->coeff(i):0)))
    }
    if((ret=putHess(task,*cjac,true,callback))!=QCQPSolver<T>::API_OK)
      return ret;
  }
  //objective
  for(sizeType i=0; i<g.size(); i++) {
    MOSEK_SAFE_CALL(MSK_putvarbound(task,i,
                                    lb&&ub?MSK_BK_RA:lb?MSK_BK_LO:ub?MSK_BK_UP:MSK_BK_FR,
                                    std::to_double(lb?lb->coeff(i):0),
                                    std::to_double(ub?ub->coeff(i):0)))
    MOSEK_SAFE_CALL(MSK_putcj(task,i,std::to_double(g[i])))
  }
  if((ret=putHess(task,H,false,callback))!=QCQPSolver<T>::API_OK)
    return ret;
  MOSEK_SAFE_CALL(MSK_putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE))
  MOSEK_SAFE_CALL(MSK_optimizetrm(task,&trmcode))
  MOSEK_SAFE_CALL(MSK_solutionsummary(task,MSK_STREAM_LOG))
  MOSEK_SAFE_CALL(MSK_getsolsta(task,MSK_SOL_ITR,&solsta))

  ret=QCQPSolver<T>::UNKNOWN;
  Cold xx=Cold::Zero(g.size());
  switch (solsta)
  {
  case MSK_SOL_STA_OPTIMAL:
    MOSEK_SAFE_CALL(MSK_getxx(task,MSK_SOL_ITR,xx.data()))
    x=xx.template cast<T>();
    ret=QCQPSolver<T>::SOLVED;
    break;
  case MSK_SOL_STA_DUAL_INFEAS_CER:
  case MSK_SOL_STA_PRIM_INFEAS_CER:
    if(callback) {
      INFO("Primal or dual infeasibility certificate found")
    }
    ret=QCQPSolver<T>::INFEASIBLE;
    break;
  case MSK_SOL_STA_UNKNOWN:
    if(callback) {
      WARNINGV("The status of the solution could not be determined. Termination code: %d\n",trmcode)
    }
    break;
  default:
    if(callback) {
      WARNING("Other solution status")
    }
    break;
  }
  if(task)
    MSK_deletetask(&task);
  return ret;
}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveL1QPTpl(Vec& x,const MAT& H,const Vec& g,const MAT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  ASSERT_MSG(QCones.empty(),"Mosek does not support second order cone with QP interface")
  QCQP_RETURN_CODE ret=QCQPSolver<T>::UNKNOWN;
  MSKtask_t task=NULL;
  MSKrescodee trmcode;
  MSKsolstae solsta;
  MOSEK_SAFE_CALL(MSK_maketask(_env,cjac?cjac->rows():0,g.size(),&task))
  MOSEK_SAFE_CALL(MSK_linkfunctotaskstream(task,MSK_STREAM_LOG,NULL,callback?printstr:NULL))
  sizeType nrVar=g.size(),nrCon=0;
  //constraint
  if(cjac) {
    for(sizeType i=0; i<cjac->rows(); i++) {
      if(lbA && lbA->coeff(i)>-DSSQPObjective<T>::infty()/2) {
        //cjac*d>=lbA
        //cjac*d-lbA>=-SLACK
        nrVar++;
        nrCon++;
      }
      if(ubA && ubA->coeff(i)< DSSQPObjective<T>::infty()/2) {
        //ubA>=cjac*d
        //ubA-cjac*d>=-SLACK
        nrVar++;
        nrCon++;
      }
    }
    MOSEK_SAFE_CALL(MSK_appendvars(task,nrVar))
    MOSEK_SAFE_CALL(MSK_appendcons(task,nrCon+1))   //the last one is trust region

    SMat m;
    STrips trips;
    m.resize(nrCon,cjac->rows());
    nrVar=g.size(),nrCon=0;
    for(sizeType i=0; i<cjac->rows(); i++) {
      if(lbA && lbA->coeff(i)>-DSSQPObjective<T>::infty()/2) {
        //cjac*d>=lbA
        //cjac*d-lbA>=-SLACK
        trips.push_back(STrip(nrCon,i,1));
        MOSEK_SAFE_CALL(MSK_putcj(task,nrVar,std::to_double(rho)))
        MOSEK_SAFE_CALL(MSK_putaij(task,nrCon,nrVar,1))
        MOSEK_SAFE_CALL(MSK_putvarbound(task,nrVar,MSK_BK_LO,0,0))
        MOSEK_SAFE_CALL(MSK_putconbound(task,nrCon,MSK_BK_LO,std::to_double(lbA->coeff(i)),0))
        nrCon++;
        nrVar++;
      }
      if(ubA && ubA->coeff(i)< DSSQPObjective<T>::infty()/2) {
        //ubA>=cjac*d
        //ubA-cjac*d>=-SLACK
        trips.push_back(STrip(nrCon,i,-1));
        MOSEK_SAFE_CALL(MSK_putcj(task,nrVar,std::to_double(rho)))
        MOSEK_SAFE_CALL(MSK_putaij(task,nrCon,nrVar,1))
        MOSEK_SAFE_CALL(MSK_putvarbound(task,nrVar,MSK_BK_LO,0,0))
        MOSEK_SAFE_CALL(MSK_putconbound(task,nrCon,MSK_BK_LO,std::to_double(-ubA->coeff(i)),0))
        nrCon++;
        nrVar++;
      }
    }
    m.setFromTriplets(trips.begin(),trips.end());
    if((ret=putHess(task,MAT(m**cjac),true,callback))!=QCQPSolver<T>::API_OK)
      return ret;
  } else {
    MOSEK_SAFE_CALL(MSK_appendvars(task,nrVar))
    MOSEK_SAFE_CALL(MSK_appendcons(task,nrCon+1))   //the last one is trust region
  }
  //objective
  for(sizeType i=0; i<g.size(); i++) {
    MOSEK_SAFE_CALL(MSK_putvarbound(task,i,
                                    lb&&ub?MSK_BK_RA:lb?MSK_BK_LO:ub?MSK_BK_UP:MSK_BK_FR,
                                    std::to_double(lb?lb->coeff(i):0),
                                    std::to_double(ub?ub->coeff(i):0)))
    MOSEK_SAFE_CALL(MSK_putcj(task,i,std::to_double(g[i])))
  }
  if((ret=putHess(task,H,false,callback))!=QCQPSolver<T>::API_OK)
    return ret;
  //trust region
  std::vector<MSKrealt> qcval;
  std::vector<MSKint32t> qcsubi,qcsubj;
  for(sizeType i=0; i<g.size(); i++) {
    qcsubi.push_back(i);
    qcsubj.push_back(i);
    qcval.push_back(2);
  }
  MOSEK_SAFE_CALL(MSK_putqconk(task,nrCon,g.size(),&qcsubi[0],&qcsubj[0],&qcval[0]))
  MOSEK_SAFE_CALL(MSK_putconbound(task,nrCon,MSK_BK_UP,std::to_double(TR),std::to_double(TR)))

  MOSEK_SAFE_CALL(MSK_putobjsense(task,MSK_OBJECTIVE_SENSE_MINIMIZE))
  MOSEK_SAFE_CALL(MSK_optimizetrm(task,&trmcode))
  MOSEK_SAFE_CALL(MSK_solutionsummary(task,MSK_STREAM_LOG))
  MOSEK_SAFE_CALL(MSK_getsolsta(task,MSK_SOL_ITR,&solsta))

  ret=QCQPSolver<T>::UNKNOWN;
  Cold xx=Cold::Zero(nrVar);
  switch (solsta)
  {
  case MSK_SOL_STA_OPTIMAL:
    MOSEK_SAFE_CALL(MSK_getxx(task,MSK_SOL_ITR,xx.data()))
    x=xx.template cast<T>();
    ret=QCQPSolver<T>::SOLVED;
    break;
  case MSK_SOL_STA_DUAL_INFEAS_CER:
  case MSK_SOL_STA_PRIM_INFEAS_CER:
    if(callback) {
      INFO("Primal or dual infeasibility certificate found")
    }
    ret=QCQPSolver<T>::INFEASIBLE;
    break;
  case MSK_SOL_STA_UNKNOWN:
    if(callback) {
      WARNINGV("The status of the solution could not be determined. Termination code: %d\n",trmcode)
    }
    break;
  default:
    if(callback) {
      WARNING("Other solution status")
    }
    break;
  }
  if(task)
    MSK_deletetask(&task);
  return ret;
}
//helper
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::putHess(void* task,const MatT& H,bool isA,bool callback) const
{
  MSKrescodee trmcode;
  std::vector<MSKint32t> qosubi;
  std::vector<MSKint32t> qosubj;
  std::vector<MSKrealt> qoval;
  for(sizeType r=0; r<H.rows(); r++)
    for(sizeType c=0; c<H.cols(); c++) {
      if(!isA && r<c)
        continue;
      qosubi.push_back(r);
      qosubj.push_back(c);
      qoval.push_back(std::to_double(H(r,c)));
    }
  if(isA) {
    MOSEK_SAFE_CALL(MSK_putaijlist(task,qosubi.size(),&qosubi[0],&qosubj[0],&qoval[0]))
  } else {
    MOSEK_SAFE_CALL(MSK_putqobj(task,qosubi.size(),&qosubi[0],&qosubj[0],&qoval[0]))
  }
  return QCQPSolver<T>::API_OK;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::putHess(void* task,const SMat& H,bool isA,bool callback) const
{
  MSKrescodee trmcode;
  std::vector<MSKint32t> qosubi;
  std::vector<MSKint32t> qosubj;
  std::vector<MSKrealt> qoval;
  for(sizeType k=0; k<H.outerSize(); ++k)
    for(typename SMat::InnerIterator it(H,k); it; ++it) {
      if(!isA && it.row()<it.col())
        continue;
      qosubi.push_back(it.row());
      qosubj.push_back(it.col());
      qoval.push_back(std::to_double(it.value()));
    }
  if(isA) {
    MOSEK_SAFE_CALL(MSK_putaijlist(task,qoval.size(),&qosubi[0],&qosubj[0],&qoval[0]))
  } else {
    MOSEK_SAFE_CALL(MSK_putqobj(task,qoval.size(),&qosubi[0],&qosubj[0],&qoval[0]))
  }
  return QCQPSolver<T>::API_OK;
}
#else
template <typename T>
QCQPSolverMosek<T>::QCQPSolverMosek() {}
template <typename T>
QCQPSolverMosek<T>::~QCQPSolverMosek() {}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveQPTpl(Vec&,const MAT&,const Vec&,const MAT*,const Vec*,const Vec*,const Vec*,const Vec*,const std::vector<Coli,Eigen::aligned_allocator<Coli>>&,bool)
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveL1QPTpl(Vec&,const MAT&,const Vec&,const MAT*,const Vec*,const Vec*,const Vec*,const Vec*,T,T,const std::vector<Coli,Eigen::aligned_allocator<Coli>>&,bool)
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
//helper
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::putHess(void* task,const MatT& H,bool isA,bool callback) const
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
// template <typename T>
// typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::putHess(void* task,const MatT& H,bool isA,bool callback) const
// {
//   FUNCTION_NOT_IMPLEMENTED
//   return QCQPSolver<T>::UNKNOWN;
// }
#endif
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveQP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveQPTpl<MatT>(x,H,g,cjac,lb,ub,lbA,ubA,QCones,callback);
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveQP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveQPTpl<SMat>(x,H,g,cjac,lb,ub,lbA,ubA,QCones,callback);
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveL1QP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveL1QPTpl<MatT>(x,H,g,cjac,lb,ub,lbA,ubA,TR,rho,QCones,callback);
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverMosek<T>::solveL1QP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveL1QPTpl<SMat>(x,H,g,cjac,lb,ub,lbA,ubA,TR,rho,QCones,callback);
}
//instance
PRJ_BEGIN
template struct QCQPSolverMosek<double>;
#ifdef ALL_TYPES
template struct QCQPSolverMosek<__float128>;
template struct QCQPSolverMosek<mpfr::mpreal>;
#endif
PRJ_END
