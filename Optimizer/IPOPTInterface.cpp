#ifdef IPOPT_SUPPORT
#include "IPOPTInterface.h"
#include <coin/IpIpoptApplication.hpp>
#include <CommonFile/Timing.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

template <typename T>
IPOPTInterface<T>::IPOPTInterface(DSSQPObjective<T>& obj,Cold& x,bool timing):_timing(timing),_obj(obj),_x(x) {}
template <typename T>
void IPOPTInterface<T>::set_xlu(const Cold& xl,const Cold& xu)
{
  if(xl.size()==0)
    _xl=NULL;
  else {
    _xl.reset(new Cold(xl));
    ASSERT(xl.size()==_obj.inputs())
  }
  if(xu.size()==0)
    _xu=NULL;
  else {
    _xu.reset(new Cold(xu));
    ASSERT(xu.size()==_obj.inputs())
  }
}
template <typename T>
void IPOPTInterface<T>::set_glu(const Cold& gl,const Cold& gu)
{
  if(gl.size()==0)
    _gl=NULL;
  else {
    _gl.reset(new Cold(gl));
    ASSERT(gl.size()==_obj.values())
  }
  if(gu.size()==0)
    _gu=NULL;
  else {
    _gu.reset(new Cold(gu));
    ASSERT(gu.size()==_obj.values());
  }
}
template <typename T>
bool IPOPTInterface<T>::get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                                     Index& nnz_h_lag, IndexStyleEnum& index_style)
{
  n=_obj.inputs();
  m=_obj.values();
  updateX(_x);
  if(_isDense) {
    nnz_jac_g=_gJacD.size();
    nnz_h_lag=_fHessD.size();
  } else {
    nnz_jac_g=_gJacS.nonZeros();
    nnz_h_lag=_fHessS.nonZeros();
  }
  index_style=TNLP::C_STYLE;
  return true;
}
template <typename T>
bool IPOPTInterface<T>::get_bounds_info(Index n, Number* x_l, Number* x_u,
                                        Index m, Number* g_l, Number* g_u)
{
  ASSERT(n==_obj.inputs() && m==_obj.values())
  if(x_l)
    Eigen::Map<Eigen::Matrix<Number,-1,1> >(x_l,n)=(_xl?*_xl:Cold::Constant(n,-getInfty())).template cast<Number>();
  if(x_u)
    Eigen::Map<Eigen::Matrix<Number,-1,1> >(x_u,n)=(_xu?*_xu:Cold::Constant(n, getInfty())).template cast<Number>();
  if(g_l)
    Eigen::Map<Eigen::Matrix<Number,-1,1> >(g_l,m)=(_gl?*_gl:Cold::Constant(m,0)).template cast<Number>();
  if(g_u)
    Eigen::Map<Eigen::Matrix<Number,-1,1> >(g_u,m)=(_gu?*_gu:Cold::Constant(m,0)).template cast<Number>();
  return true;
}
template <typename T>
bool IPOPTInterface<T>::get_starting_point(Index n, bool init_x, Number* x,
    bool init_z, Number* z_L, Number* z_U,
    Index m, bool init_lambda,
    Number* lambda)
{
  ASSERT(init_x && n==_x.size())
  ASSERT(!init_z)
  ASSERT(!init_lambda)
  if(x)
    Eigen::Map<Eigen::Matrix<Number,-1,1> >(x,n)=_x.cast<Number>();
  return true;
}
template <typename T>
bool IPOPTInterface<T>::eval_f(Index n, const Number* x, bool new_x, Number& obj_value)
{
  if(new_x && x)
    updateX(Eigen::Map<const Eigen::Matrix<Number,-1,1> >(x,n).cast<scalarD>());
  obj_value=_E;
  return true;
}
template <typename T>
bool IPOPTInterface<T>::eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f)
{
  if(new_x && x)
    updateX(Eigen::Map<const Eigen::Matrix<Number,-1,1> >(x,n).cast<scalarD>());
  if(grad_f)
    Eigen::Map<Eigen::Matrix<Number,-1,1> >(grad_f,n)=_fGrad.cast<Number>();
  return true;
}
template <typename T>
bool IPOPTInterface<T>::eval_g(Index n, const Number* x, bool new_x, Index m, Number* g)
{
  if(new_x && x)
    updateX(Eigen::Map<const Eigen::Matrix<Number,-1,1> >(x,n).cast<scalarD>());
  if(g)
    Eigen::Map<Eigen::Matrix<Number,-1,1> >(g,m)=_gVec.cast<Number>();
  return true;
}
template <typename T>
bool IPOPTInterface<T>::eval_jac_g(Index n, const Number* x, bool new_x,
                                   Index m, Index nele_jac, Index* iRow, Index *jCol,
                                   Number* values)
{
  if(new_x && x)
    updateX(Eigen::Map<const Eigen::Matrix<Number,-1,1> >(x,n).cast<scalarD>());
  sizeType nnz=0;
  if(_isDense) {
    for(sizeType r=0; r<_gJacD.rows(); r++)
      for(sizeType c=0; c<_gJacD.cols(); c++) {
        if(iRow) {
          *iRow=r;
          iRow++;
        }
        if(jCol) {
          *jCol=c;
          jCol++;
        }
        if(values) {
          *values=_gJacD(r,c);
          values++;
        }
        nnz++;
      }
  } else {
    for(sizeType k=0; k<_gJacS.outerSize(); ++k)
      for(DSSQPObjective<scalarD>::SMat::InnerIterator it(_gJacS,k); it; ++it) {
        if(iRow) {
          *iRow=it.row();
          iRow++;
        }
        if(jCol) {
          *jCol=it.col();
          jCol++;
        }
        if(values) {
          *values=it.value();
          values++;
        }
        nnz++;
      }
  }
  ASSERT_MSGV(nnz == nele_jac,"NNZ in hessian(%ld) does not match that declared in get_nlp_info(%ld)!",nnz,nele_jac)
  return true;
}
template <typename T>
bool IPOPTInterface<T>::eval_h(Index n, const Number* x, bool new_x,
                               Number obj_factor, Index m, const Number* lambda,
                               bool new_lambda, Index nele_hess, Index* iRow,
                               Index* jCol, Number* values)
{
  if(new_x && x)
    updateX(Eigen::Map<const Eigen::Matrix<Number,-1,1> >(x,n).cast<scalarD>());
  sizeType nnz=0;
  if(_isDense) {
    for(sizeType r=0; r<_fHessD.rows(); r++)
      for(sizeType c=0; c<_fHessD.cols(); c++) {
        if(iRow) {
          *iRow=r;
          iRow++;
        }
        if(jCol) {
          *jCol=c;
          jCol++;
        }
        if(values) {
          *values=_fHessD(r,c);
          values++;
        }
        nnz++;
      }
  } else {
    for(sizeType k=0; k<_fHessS.outerSize(); ++k)
      for(DSSQPObjective<scalarD>::SMat::InnerIterator it(_fHessS,k); it; ++it) {
        if(iRow) {
          *iRow=it.row();
          iRow++;
        }
        if(jCol) {
          *jCol=it.col();
          jCol++;
        }
        if(values) {
          *values=it.value();
          values++;
        }
        nnz++;
      }
  }
  ASSERT_MSGV(nnz == nele_hess,"NNZ in hessian(%ld) does not match that declared in get_nlp_info(%ld)!",nnz,nele_hess)
  return true;
}
template <typename T>
void IPOPTInterface<T>::finalize_solution(SolverReturn status,
    Index n, const Number* x, const Number* z_L, const Number* z_U,
    Index m, const Number* g, const Number* lambda,
    Number obj_value,
    const IpoptData* ip_data,
    IpoptCalculatedQuantities* ip_cq)
{
  _x=Eigen::Map<const Eigen::Matrix<Number,-1,1> >(x,n).cast<scalarD>();
}
template <typename T>
bool IPOPTInterface<T>::optimize(typename DSSQPObjective<T>::Vec& x, DSSQPObjective<T>& obj, scalarD tolG, sizeType maxIter, scalarD derivCheck, scalarD muInit, bool dense, int callback)
{
  Cold xd=x.unaryExpr([](const T& in) {
    return std::to_double(in);
  });
  Ipopt::SmartPtr<IPOPTInterface<T>> mynlp=new IPOPTInterface<T>(obj,xd);
  Ipopt::SmartPtr<Ipopt::IpoptApplication> app=IpoptApplicationFactory();
  mynlp->_isDense=dense;
  mynlp->set_xlu(obj.lb().unaryExpr([](const T& in) {
    return std::to_double(in);
  }),obj.ub().unaryExpr([](const T& in) {
    return std::to_double(in);
  }));
  mynlp->set_glu(obj.gl().unaryExpr([](const T& in) {
    return std::to_double(in);
  }),obj.gu().unaryExpr([](const T& in) {
    return std::to_double(in);
  }));
  app->RethrowNonIpoptException(true);
  if(derivCheck>0) {
    app->Options()->SetStringValue("derivative_test","first-order");
    app->Options()->SetNumericValue("derivative_test_tol",derivCheck);
  } else {
    app->Options()->SetStringValue("derivative_test","none");
  }
  app->Options()->SetNumericValue("tol",tolG);
  if(muInit>0)
    app->Options()->SetNumericValue("mu_init",muInit);
  app->Options()->SetIntegerValue("max_iter",maxIter);
  app->Options()->SetIntegerValue("print_level",callback<0?0:callback>12?12:callback);
  app->Options()->SetStringValue("hessian_approximation","limited-memory");
  //app->Options()->SetStringValue("mu_strategy","adaptive");
  app->Options()->SetStringValue("mu_strategy","monotone");
  Ipopt::ApplicationReturnStatus status=app->Initialize();
  status=app->OptimizeTNLP(mynlp);
  x=xd.template cast<T>();
  //ASSERT_MSG(status==Ipopt::Solve_Succeeded,"Error in Ipopt optimization!")
  return status==Ipopt::Solve_Succeeded;
}
template <typename T>
void IPOPTInterface<T>::updateX(const Cold& x)
{
  if(_timing)
    TBEG("IPOPTInterface<T>::updateX");
  typename DSSQPObjective<T>::DMat gJacD,fHessD;
  typename DSSQPObjective<T>::SMat gJacS,fHessS;
  typename DSSQPObjective<T>::Vec gVec,fGrad;
  if(_isDense) {
    _obj(x.template cast<T>(),gVec,&gJacD);
    _gVec=gVec.unaryExpr([](const T& in) {
      return std::to_double(in);
    });
    _gJacD=gJacD.unaryExpr([](const T& in) {
      return std::to_double(in);
    });
    _E=std::to_double(_obj(x.template cast<T>(),&fGrad));
    _fGrad=fGrad.unaryExpr([](const T& in) {
      return std::to_double(in);
    });
    _fHessD.resize(x.size(),x.size());
  } else {
    _obj(x.template cast<T>(),gVec,&gJacS);
    _gVec=gVec.unaryExpr([](const T& in) {
      return std::to_double(in);
    });
    _gJacS=gJacS.unaryExpr([](const T& in) {
      return std::to_double(in);
    });
    _E=std::to_double(_obj(x.template cast<T>(),&fGrad));
    _fGrad=fGrad.unaryExpr([](const T& in) {
      return std::to_double(in);
    });
    _fHessS.resize(x.size(),x.size());
  }
  if(_timing)
    TEND();
}
//instance
PRJ_BEGIN
template struct IPOPTInterface<double>;
#ifdef ALL_TYPES
template struct IPOPTInterface<__float128>;
template struct IPOPTInterface<mpfr::mpreal>;
#endif
PRJ_END

#endif
