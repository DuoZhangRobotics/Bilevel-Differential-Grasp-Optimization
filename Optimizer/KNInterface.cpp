#include "KNInterface.h"
#ifdef KNITRO_SUPPORT
#include <CommonFile/Timing.h>
#include <Utils/DebugGradient.h>
#include <Utils/Scalar.h>

USE_PRJ_NAMESPACE

//KTRInterface
template <typename T>
KNInterface<T>::KNInterface(DSSQPObjective<T>& g,bool timing)
  :knitro::KNProblem(g.inputs(),g.values()),_obj(g),_timing(timing)
{
  setObjectiveProperties();
  setVariableProperties();
  setConstraintProperties();
  //set derivatives
  Vec x;
  Vec fvec;
  SMat fjac;
  x.setZero(g.inputs());
  std::vector<int> rows,cols,ids;
  g.DSSQPObjective<T>::operator()(x,fvec,&fjac);
  for(sizeType k=0,off=0; k<fjac.outerSize(); ++k)
    for(typename DSSQPObjective<T>::SMat::InnerIterator it(fjac,k); it; ++it,off++) {
      //setJacIndexCons(off,it.row());
      //setJacIndexVars(off,it.col());
      rows.push_back(it.row());
      cols.push_back(it.col());
    }
  for(sizeType i=0; i<fjac.rows(); i++)
    ids.push_back(i);
  setMainCallbackCstIndexes(ids);
  setJacNnzPattern(knitro::KNSparseMatrixStructure(rows,cols));
  setObjEvalCallback(&KNInterface<T>::evaluateFC);
  setGradEvalCallback(&KNInterface<T>::evaluateGA);
  KNProblem::getMainNonLinearStructure().getEvalCallback().setParams(this);
}
template <typename T>
int KNInterface<T>::evaluateFC(KN_context_ptr,CB_context_ptr,KN_eval_request_ptr const evalRequest,KN_eval_result_ptr const evalResult,void* const userParams)
{
  KNInterface<T>* I=(KNInterface<T>*)userParams;
  if(I->_timing)
    TBEG("KTRInterface<T>::evaluateFC");
  Vec C=Vec::Zero(I->_obj.values());
  Eigen::Map<const Eigen::Matrix<double,-1,1>> xMap(evalRequest->x,I->_obj.inputs());
  evalResult->obj[0]=std::to_double(I->_obj(xMap.template cast<T>(),(Vec*)NULL));
  I->_obj(xMap.template cast<T>(),C,(STrips*)NULL);
  Eigen::Map<Eigen::Matrix<double,-1,1>>(evalResult->c,I->_obj.values())=C.unaryExpr([](T in) {
    return (double)std::to_double(in);
  });
  if(I->_timing)
    TEND();
  return 0;
}
template <typename T>
int KNInterface<T>::evaluateGA(KN_context_ptr,CB_context_ptr,KN_eval_request_ptr const evalRequest,KN_eval_result_ptr const evalResult,void* const userParams)
{
  KNInterface<T>* I=(KNInterface<T>*)userParams;
  if(I->_timing)
    TBEG("KTRInterface<T>::evaluateGA");
  Vec C=Vec::Zero(I->_obj.values());
  Vec G=Vec::Zero(I->_obj.inputs());
  typename DSSQPObjective<T>::SMat JAC;
  Eigen::Map<const Eigen::Matrix<double,-1,1>> xMap(evalRequest->x,I->_obj.inputs());
  I->_obj(xMap.cast<T>(),&G);
  Eigen::Map<Eigen::Matrix<double,-1,1>>(evalResult->objGrad,I->_obj.inputs())=G.unaryExpr([](T in) {
    return (double)std::to_double(in);
  });
  I->_obj.DSSQPObjective<T>::operator()(xMap.template cast<T>(),C,&JAC);
  for (sizeType k=0, off=0; k<JAC.outerSize(); ++k)
    for (typename DSSQPObjective<T>::SMat::InnerIterator it(JAC,k); it; ++it, off++)
      evalResult->jac[off]=std::to_double(it.value());
  if(I->_timing)
    TEND();
  return 0;
}
template <typename T>
bool KNInterface<T>::solve(bool callback,scalarD checkDeriv,scalarD ftol,scalarD xtol,bool useDirect,sizeType maxIter,scalarD muInit,Vec& x,Vec* fvec)
{
  knitro::KNSolver solver(this,KN_GRADOPT_EXACT,KN_HESSOPT_LBFGS);
  solver.setParam(KN_PARAM_ALGORITHM,useDirect?KN_ALG_BAR_DIRECT:KN_ALG_BAR_CG);
  //solver.setParam(KN_PARAM_ALGORITHM,useDirect?KN_ALG_ACT_SQP:KN_ALG_ACT_CG);
  solver.setParam(KN_PARAM_BAR_FEASIBLE,KN_BAR_FEASIBLE_GET_STAY);
  solver.setParam(KN_PARAM_OUTMODE,callback?KN_OUTMODE_SCREEN:KN_OUTMODE_FILE);
  solver.setParam(KN_PARAM_CG_MAXIT,50);
  //solver.setParam(KN_PARAM_CG_PRECOND,KN_CG_PRECOND_CHOL);
  solver.setParam(KN_PARAM_OUTLEV,KN_OUTLEV_ITER);
  if(checkDeriv>0) {
    solver.setParam(KN_PARAM_DERIVCHECK,KN_DERIVCHECK_FIRST);
    solver.setParam(KN_PARAM_DERIVCHECK_TOL,(double)checkDeriv);
  }
  if(ftol>0) {
    solver.setParam(KN_PARAM_FTOL,(double)ftol);
    solver.setParam(KN_PARAM_XTOL,(double)xtol);
  }
  if(maxIter>0)
    solver.setParam(KN_PARAM_MAXIT,(int)maxIter);
  if(muInit>0) {
    solver.setParam(KN_PARAM_BAR_PENRULE,KN_BAR_PENRULE_SINGLE);
    solver.setParam(KN_PARAM_INITPENALTY,(double)muInit);
  }
  solver.initProblem();
  int solveStatus=solver.solve();
  bool succ=(-103<=solveStatus&&solveStatus<=-100) || (-406<=solveStatus&&solveStatus<=-400) || (solveStatus>=0);
  std::vector<double> primal,dual;
  solver.getSolution(primal,dual);
  x=Eigen::Map<const Eigen::Matrix<double,-1,1>>(&primal[0],primal.size()).template cast<T>();
  if(fvec)
    _obj.operator()(x,*fvec,(STrips*)NULL);
  return succ;
}
//helper
template <typename T>
void KNInterface<T>::setObjectiveProperties() {
  setObjType(KN_OBJTYPE_GENERAL);
  setObjGoal(KN_OBJGOAL_MINIMIZE);
}
template <typename T>
void KNInterface<T>::setVariableProperties() {
  const Vec& x=_obj.init();
  Vec LB=_obj.lb(), UB=_obj.ub();
  getXInitial().clear();
  getVarLoBnds().clear();
  getVarUpBnds().clear();
  for (sizeType i=0; i<x.size(); i++) {
    getXInitial().add(i,std::to_double(x[i]));
    getVarLoBnds().add(i,std::to_double(LB[i]));
    getVarUpBnds().add(i,std::to_double(UB[i]));
  }
}
template <typename T>
void KNInterface<T>::setConstraintProperties() {
  Vec GL=_obj.gl(), GU=_obj.gu();
  getConTypes().clear();
  getConLoBnds().clear();
  getConUpBnds().clear();
  for (sizeType i=0; i<GL.size(); i++) {
    getConTypes().add(i,KN_CONTYPE_GENERAL);
    getConLoBnds().add(i,std::to_double(GL[i]));
    getConUpBnds().add(i,std::to_double(GU[i]));
  }
}
//instance
PRJ_BEGIN
template struct KNInterface<double>;
#ifdef ALL_TYPES
template struct KNInterface<__float128>;
template struct KNInterface<mpfr::mpreal>;
#endif
PRJ_END

#endif
