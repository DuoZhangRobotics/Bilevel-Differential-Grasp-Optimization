#ifndef KNITRO_INTERFACE_H
#define KNITRO_INTERFACE_H
#ifdef KNITRO_SUPPORT

#include "DSSQPObjective.h"
#include "KNSolver.h"

PRJ_BEGIN

template <typename T>
class KNInterface : public knitro::KNProblem
{
public:
  typedef typename DSSQPObjective<T>::Vec Vec;
  typedef typename DSSQPObjective<T>::DMat DMat;
  typedef typename DSSQPObjective<T>::SMat SMat;
  typedef typename DSSQPObjective<T>::STrip STrip;
  typedef typename DSSQPObjective<T>::STrips STrips;
  KNInterface(DSSQPObjective<T>& obj,bool timing=false);
  static int evaluateFC(KN_context_ptr,CB_context_ptr,KN_eval_request_ptr const evalRequest,KN_eval_result_ptr const evalResult,void* const userParams);
  static int evaluateGA(KN_context_ptr,CB_context_ptr,KN_eval_request_ptr const evalRequest,KN_eval_result_ptr const evalResult,void* const userParams);
  bool solve(bool callback,scalarD checkDeriv,scalarD ftol,scalarD xtol,bool useDirect,sizeType maxIter,scalarD muInit,Vec& x,Vec* fvec=NULL);
protected:
  void setObjectiveProperties();
  void setVariableProperties();
  void setConstraintProperties();
  DSSQPObjective<T>& _obj;
  bool _timing;
};

PRJ_END

#endif
#endif
