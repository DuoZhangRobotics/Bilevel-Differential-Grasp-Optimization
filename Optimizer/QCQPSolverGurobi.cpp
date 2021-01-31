#include "DSSQPObjective.h"
#include "QCQPSolverGurobi.h"

USE_PRJ_NAMESPACE

#ifdef GUROBI_SUPPORT
#include <gurobi_c.h>
#define GRB_SAFE_CALL(FUNC)\
if((err=FUNC)!=0) {\
  if(callback){\
    WARNINGV("Gurobi failed (code=%d) while calling: %s",err,STRINGIFY_OMP(FUNC))\
  }\
  if(model)\
    GRBfreemodel(model);\
  if(err==GRB_ERROR_Q_NOT_PSD)\
    return QCQPSolver<T>::NOT_POSITIVE_DEFINITE;\
  else return QCQPSolver<T>::UNKNOWN;\
}
template <typename T>
QCQPSolverGurobi<T>::QCQPSolverGurobi():_env(NULL) {
  GRBenv* env=(GRBenv*)_env;
  int err=GRBloadenv(&env,"");
  if(err) {
    ASSERT_MSGV(false,"Gurobi failed (code=%d) while calling: %s",err,STRINGIFY_OMP(GRBloadenv))
  }
  _env=env;
}
template <typename T>
QCQPSolverGurobi<T>::~QCQPSolverGurobi() {
  if(_env)
    GRBfreeenv((GRBenv*)_env);
}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveQPTpl(Vec& x,const MAT& H,const Vec& g,const MAT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  int err;
  GRBmodel* model=NULL;
  typename QCQPSolver<T>::QCQP_RETURN_CODE ret=QCQPSolver<T>::UNKNOWN;
  GRB_SAFE_CALL(GRBnewmodel((GRBenv*)_env,&model,"qcp",0,NULL,NULL,NULL,NULL,NULL))
  //build objective
  Cold gd=g.unaryExpr([&](T in) {
    return std::to_double(in);
  }),lbd,ubd;
  if(lb)
    lbd=lb->unaryExpr([&](T in) {
    return std::to_double(in);
  });
  if(ub)
    ubd=ub->unaryExpr([&](T in) {
    return std::to_double(in);
  });
  GRB_SAFE_CALL(GRBaddvars(model,g.size(),0,NULL,NULL,NULL,gd.data(),lb?lbd.data():NULL,ub?ubd.data():NULL,NULL,NULL))
  if((ret=putHess(model,H,callback))!=QCQPSolver<T>::API_OK)
    return QCQPSolver<T>::UNKNOWN;
  //build constraint
  std::vector<Eigen::Triplet<double,int>> tripsRelaxed;
  if(cjac && (ret=putCons(model,*cjac,lbA,ubA,tripsRelaxed,g.size(),callback))!=QCQPSolver<T>::API_OK)
    return QCQPSolver<T>::UNKNOWN;
  //quadratic cone
  if((ret=putQCone(model,x,g,QCones,g.size(),callback))!=QCQPSolver<T>::API_OK)
    return QCQPSolver<T>::UNKNOWN;

  //optimize
  int optimstatus;
  GRB_SAFE_CALL(GRBsetintparam((GRBenv*)_env,"OutputFlag",callback?1:0))
  GRB_SAFE_CALL(GRBsetintattr(model,GRB_INT_ATTR_MODELSENSE,GRB_MINIMIZE))
  GRB_SAFE_CALL(GRBoptimize(model))
  GRB_SAFE_CALL(GRBgetintattr(model,GRB_INT_ATTR_STATUS,&optimstatus))

  ret=QCQPSolver<T>::UNKNOWN;
  if(optimstatus == GRB_OPTIMAL) {
    ret=QCQPSolver<T>::SOLVED;
    Cold sol;
    sol.resize(g.size());
    GRB_SAFE_CALL(GRBgetdblattrarray(model,GRB_DBL_ATTR_X,0,x.size(),sol.data()))
    //INFOV("Optimal objective: %.4e x=%.2f, y=%.2f, z=%.2f",objval,sol[0],sol[1],sol[2])
    x=sol.template cast<T>();
  } else if(optimstatus == GRB_INF_OR_UNBD) {
    ret=QCQPSolver<T>::INFEASIBLE;
    WARNING("Model is infeasible or unbounded");
  } else {
    WARNING("Optimization was stopped early");
  }
  if(model)
    GRBfreemodel(model);
  return ret;
}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveL1QPTpl(Vec& x,const MAT& H,const Vec& g,const MAT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  int err;
  GRBmodel* model=NULL;
  typename QCQPSolver<T>::QCQP_RETURN_CODE ret=QCQPSolver<T>::UNKNOWN;
  GRB_SAFE_CALL(GRBnewmodel((GRBenv*)_env,&model,"qcp",0,NULL,NULL,NULL,NULL,NULL))
  //build objective
  Cold gd=g.unaryExpr([&](T in) {
    return std::to_double(in);
  }),lbd,ubd;
  if(lb)
    lbd=lb->unaryExpr([&](T in) {
    return std::to_double(in);
  });
  if(ub)
    ubd=ub->unaryExpr([&](T in) {
    return std::to_double(in);
  });
  GRB_SAFE_CALL(GRBaddvars(model,g.size(),0,NULL,NULL,NULL,gd.data(),lb?lbd.data():NULL,ub?ubd.data():NULL,NULL,NULL))
  if((ret=putHess(model,H,callback))!=QCQPSolver<T>::API_OK)
    return QCQPSolver<T>::UNKNOWN;
  //build constraint
  int nrVar=g.size(),nrCon=0;
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

    SMat m;
    STrips trips;
    Vec lbARelaxed=Vec::Zero(nrCon);
    Cold rhos=Cold::Constant(nrCon,std::to_double(rho));
    std::vector<Eigen::Triplet<double,int>> tripsRelaxed;
    m.resize(nrCon,cjac->rows());
    nrVar=g.size(),nrCon=0;
    for(sizeType i=0; i<cjac->rows(); i++) {
      if(lbA && lbA->coeff(i)>-DSSQPObjective<T>::infty()/2) {
        //cjac*d>=lbA
        //cjac*d-lbA>=-SLACK
        trips.push_back(STrip(nrCon,i,1));
        tripsRelaxed.push_back(Eigen::Triplet<double,int>(nrCon,nrVar,1));
        lbARelaxed[nrCon]= lbA->coeff(i);
        nrCon++;
        nrVar++;
      }
      if(ubA && ubA->coeff(i)< DSSQPObjective<T>::infty()/2) {
        //ubA>=cjac*d
        //ubA-cjac*d>=-SLACK
        trips.push_back(STrip(nrCon,i,-1));
        tripsRelaxed.push_back(Eigen::Triplet<double,int>(nrCon,nrVar,1));
        lbARelaxed[nrCon]=-ubA->coeff(i);
        nrCon++;
        nrVar++;
      }
    }
    m.setFromTriplets(trips.begin(),trips.end());
    GRB_SAFE_CALL(GRBaddvars(model,nrCon,0,NULL,NULL,NULL,rhos.data(),NULL,NULL,NULL,NULL))
    if((ret=putCons(model,MAT(m**cjac),&lbARelaxed,NULL,tripsRelaxed,nrVar,callback))!=QCQPSolver<T>::API_OK)
      return QCQPSolver<T>::UNKNOWN;
  }
  //trust region
  std::vector<int> qrow(g.size()),qcol(g.size());
  Cold qval=Cold::Ones(g.size());
  for(sizeType i=0; i<g.size(); i++)
    qrow[i]=qcol[i]=i;
  GRB_SAFE_CALL(GRBaddqconstr(model,0,NULL,NULL,qrow.size(),&qrow[0],&qcol[0],qval.data(),GRB_LESS_EQUAL,std::to_double(TR),NULL));
  //quadratic cone
  if((ret=putQCone(model,x,g,QCones,nrVar,callback))!=QCQPSolver<T>::API_OK)
    return QCQPSolver<T>::UNKNOWN;

  //optimize
  int optimstatus;
  GRB_SAFE_CALL(GRBsetintparam((GRBenv*)_env,"OutputFlag",callback?1:0))
  GRB_SAFE_CALL(GRBsetintattr(model,GRB_INT_ATTR_MODELSENSE,GRB_MINIMIZE))
  GRB_SAFE_CALL(GRBoptimize(model))
  GRB_SAFE_CALL(GRBgetintattr(model,GRB_INT_ATTR_STATUS,&optimstatus))

  ret=QCQPSolver<T>::UNKNOWN;
  if(optimstatus == GRB_OPTIMAL) {
    ret=QCQPSolver<T>::SOLVED;
    Cold sol;
    sol.resize(g.size());
    GRB_SAFE_CALL(GRBgetdblattrarray(model,GRB_DBL_ATTR_X,0,x.size(),sol.data()))
    //INFOV("Optimal objective: %.4e x=%.2f, y=%.2f, z=%.2f",objval,sol[0],sol[1],sol[2])
    x=sol.template cast<T>();
  } else if(optimstatus == GRB_INF_OR_UNBD) {
    ret=QCQPSolver<T>::INFEASIBLE;
    WARNING("Model is infeasible or unbounded");
  } else {
    WARNINGV("Optimization was stopped early code=%d",optimstatus);
  }
  if(model)
    GRBfreemodel(model);
  return ret;
}
//helper
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::putHess(void* m,const MatT& H,bool callback) const
{
  int err;
  GRBmodel* model=(GRBmodel*)m;
  std::vector<int> qosubi;
  std::vector<int> qosubj;
  std::vector<double> qoval;
  for(sizeType r=0; r<H.rows(); r++)
    for(sizeType c=0; c<H.cols(); c++) {
      qosubi.push_back(r);
      qosubj.push_back(c);
      qoval.push_back(std::to_double(H(r,c))*0.5f);
    }
  GRB_SAFE_CALL(GRBaddqpterms((GRBmodel*)model,qosubi.size(),&qosubi[0],&qosubj[0],&qoval[0]))
  return QCQPSolver<T>::API_OK;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::putHess(void* m,const SMat& H,bool callback) const
{
  int err;
  GRBmodel* model=(GRBmodel*)m;
  std::vector<int> qosubi;
  std::vector<int> qosubj;
  std::vector<double> qoval;
  for(sizeType k=0; k<H.outerSize(); ++k)
    for(typename SMat::InnerIterator it(H,k); it; ++it) {
      qosubi.push_back(it.row());
      qosubj.push_back(it.col());
      qoval.push_back(std::to_double(it.value()*0.5f));
    }
  GRB_SAFE_CALL(GRBaddqpterms((GRBmodel*)model,qosubi.size(),&qosubi[0],&qosubj[0],&qoval[0]))
  return QCQPSolver<T>::API_OK;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::putCons(void* m,const MatT& cjac,const Vec* lbA,const Vec* ubA,std::vector<Eigen::Triplet<double,int>>& trips,int nrVar,bool callback) const
{
  int err;
  GRBmodel* model=(GRBmodel*)m;
  for(sizeType r=0; r<cjac.rows(); r++)
    for(sizeType c=0; c<cjac.cols(); c++)
      trips.push_back(Eigen::Triplet<double,int>(r,c,std::to_double(cjac(r,c))));
  Eigen::SparseMatrix<double,Eigen::RowMajor,int> cjacGRB;
  cjacGRB.resize(cjac.rows(),nrVar);
  cjacGRB.setFromTriplets(trips.begin(),trips.end());
  if(lbA) {
    std::vector<char> sense(cjac.rows(),GRB_GREATER_EQUAL);
    Cold lbAd=lbA->unaryExpr([&](T in) {
      return std::to_double(in);
    });
    GRB_SAFE_CALL(GRBaddconstrs(model,cjacGRB.rows(),cjacGRB.nonZeros(),cjacGRB.outerIndexPtr(),cjacGRB.innerIndexPtr(),cjacGRB.valuePtr(),&sense[0],lbAd.data(),NULL))
  }
  if(ubA) {
    std::vector<char> sense(cjac.rows(),GRB_LESS_EQUAL);
    Cold ubAd=ubA->unaryExpr([&](T in) {
      return std::to_double(in);
    });
    GRB_SAFE_CALL(GRBaddconstrs(model,cjacGRB.rows(),cjacGRB.nonZeros(),cjacGRB.outerIndexPtr(),cjacGRB.innerIndexPtr(),cjacGRB.valuePtr(),&sense[0],ubAd.data(),NULL))
  }
  return QCQPSolver<T>::API_OK;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::putCons(void* m,const SMat& cjac,const Vec* lbA,const Vec* ubA,std::vector<Eigen::Triplet<double,int>>& trips,int nrVar,bool callback) const
{
  int err;
  GRBmodel* model=(GRBmodel*)m;
  for(sizeType k=0; k<cjac.outerSize(); ++k)
    for(typename SMat::InnerIterator it(cjac,k); it; ++it)
      trips.push_back(Eigen::Triplet<double,int>(it.row(),it.col(),std::to_double(it.value())));
  Eigen::SparseMatrix<double,Eigen::RowMajor,int> cjacGRB;
  cjacGRB.resize(cjac.rows(),nrVar);
  cjacGRB.setFromTriplets(trips.begin(),trips.end());
  if(lbA) {
    std::vector<char> sense(cjac.rows(),GRB_GREATER_EQUAL);
    Cold lbAd=lbA->unaryExpr([&](T in) {
      return std::to_double(in);
    });
    GRB_SAFE_CALL(GRBaddconstrs(model,cjacGRB.rows(),cjacGRB.nonZeros(),cjacGRB.outerIndexPtr(),cjacGRB.innerIndexPtr(),cjacGRB.valuePtr(),&sense[0],lbAd.data(),NULL))
  }
  if(ubA) {
    std::vector<char> sense(cjac.rows(),GRB_LESS_EQUAL);
    Cold ubAd=ubA->unaryExpr([&](T in) {
      return std::to_double(in);
    });
    GRB_SAFE_CALL(GRBaddconstrs(model,cjacGRB.rows(),cjacGRB.nonZeros(),cjacGRB.outerIndexPtr(),cjacGRB.innerIndexPtr(),cjacGRB.valuePtr(),&sense[0],ubAd.data(),NULL))
  }
  return QCQPSolver<T>::API_OK;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::putQCone(void* m,const Vec& x,const Vec& g,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,int nrVar,bool callback) const
{
  int err;
  GRBmodel* model=(GRBmodel*)m;
  std::vector<int> qrow,qcol;
  Cold qval=Cold::Ones(g.size()),cval=Vec2d(1,-1);
  std::vector<int> cind(2);

  //quadratic cone
  for(sizeType q=0; q<(sizeType)QCones.size(); q++) {
    Cold lb=Cold::Constant(QCones[q].size(),-GRB_INFINITY);
    Cold ub=Cold::Constant(QCones[q].size(), GRB_INFINITY);
    lb[0]=0;
    std::vector<char> types(QCones[q].size(),GRB_CONTINUOUS);
    GRB_SAFE_CALL(GRBaddvars(model,QCones[q].size(),0,NULL,NULL,NULL,NULL,lb.data(),ub.data(),&types[0],NULL))

    qrow.resize(QCones[q].size());
    qcol.resize(QCones[q].size());
    qval.setOnes(QCones[q].size());
    for(sizeType i=0; i<QCones[q].size(); i++) {
      qrow[i]=qcol[i]=i+nrVar;
      cind[0]=nrVar+i;
      cind[1]=QCones[q][i];
      GRB_SAFE_CALL(GRBaddconstr(model,2,&cind[0],cval.data(),GRB_EQUAL,std::to_double(x[QCones[q][i]]),NULL))
    }
    qval[0]=-1;
    GRB_SAFE_CALL(GRBaddqconstr(model,0,NULL,NULL,qrow.size(),&qrow[0],&qcol[0],qval.data(),GRB_LESS_EQUAL,0,NULL))
    nrVar+=QCones[q].size();
  }
  return QCQPSolver<T>::API_OK;
}
#else
template <typename T>
QCQPSolverGurobi<T>::QCQPSolverGurobi() {}
template <typename T>
QCQPSolverGurobi<T>::~QCQPSolverGurobi() {}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveQPTpl(Vec&,const MAT&,const Vec&,const MAT*,const Vec*,const Vec*,const Vec*,const Vec*,const std::vector<Coli,Eigen::aligned_allocator<Coli>>&,bool)
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
//consistent QP
template <typename T>
template <typename MAT>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveL1QPTpl(Vec&,const MAT&,const Vec&,const MAT*,const Vec*,const Vec*,const Vec*,const Vec*,T,T,const std::vector<Coli,Eigen::aligned_allocator<Coli>>&,bool)
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
//helper
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::putHess(void*,const MatT&,bool) const
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::putHess(void*,const SMat&,bool) const
{
  FUNCTION_NOT_IMPLEMENTED
  return QCQPSolver<T>::UNKNOWN;
}
#endif
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveQP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveQPTpl<MatT>(x,H,g,cjac,lb,ub,lbA,ubA,QCones,callback);
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveQP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveQPTpl<SMat>(x,H,g,cjac,lb,ub,lbA,ubA,QCones,callback);
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveL1QP(Vec& x,const MatT& H,const Vec& g,const MatT* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveL1QPTpl<MatT>(x,H,g,cjac,lb,ub,lbA,ubA,TR,rho,QCones,callback);
}
template <typename T>
typename QCQPSolver<T>::QCQP_RETURN_CODE QCQPSolverGurobi<T>::solveL1QP(Vec& x,const SMat& H,const Vec& g,const SMat* cjac,const Vec* lb,const Vec* ub,const Vec* lbA,const Vec* ubA,T TR,T rho,const std::vector<Coli,Eigen::aligned_allocator<Coli>>& QCones,bool callback)
{
  return solveL1QPTpl<SMat>(x,H,g,cjac,lb,ub,lbA,ubA,TR,rho,QCones,callback);
}
//instance
PRJ_BEGIN
template struct QCQPSolverGurobi<double>;
#ifdef ALL_TYPES
template struct QCQPSolverGurobi<__float128>;
template struct QCQPSolverGurobi<mpfr::mpreal>;
#endif
PRJ_END
