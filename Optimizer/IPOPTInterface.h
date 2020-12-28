#ifndef IPOPT_INTERFACE_H
#define IPOPT_INTERFACE_H
#ifdef IPOPT_SUPPORT

#include <coin/IpTNLP.hpp>
#include "DSSQPObjective.h"

PRJ_BEGIN

template <typename T>
class IPOPTInterface : public Ipopt::TNLP
{
  typedef Ipopt::IpoptCalculatedQuantities IpoptCalculatedQuantities;
  typedef Ipopt::SolverReturn SolverReturn;
  typedef Ipopt::IpoptData IpoptData;
  typedef Ipopt::Number Number;
  typedef Ipopt::Index Index;
public:
  IPOPTInterface(DSSQPObjective<T>& obj,Cold& x,bool timing=false);
  virtual void set_xlu(const Cold& xl,const Cold& xu);
  virtual void set_glu(const Cold& gl,const Cold& gu);
  virtual bool get_nlp_info(Index& n, Index& m, Index& nnz_jac_g,
                            Index& nnz_h_lag, IndexStyleEnum& index_style) override;
  virtual bool get_bounds_info(Index n, Number* x_l, Number* x_u,
                               Index m, Number* g_l, Number* g_u) override;
  virtual bool get_starting_point(Index n, bool init_x, Number* x,
                                  bool init_z, Number* z_L, Number* z_U,
                                  Index m, bool init_lambda,
                                  Number* lambda) override;
  virtual bool eval_f(Index n, const Number* x, bool new_x, Number& obj_value) override;
  virtual bool eval_grad_f(Index n, const Number* x, bool new_x, Number* grad_f) override;
  virtual bool eval_g(Index n, const Number* x, bool new_x, Index m, Number* g) override;
  virtual bool eval_jac_g(Index n, const Number* x, bool new_x,
                          Index m, Index nele_jac, Index* iRow, Index *jCol,
                          Number* values) override;
  virtual bool eval_h(Index n, const Number* x, bool new_x,
                      Number obj_factor, Index m, const Number* lambda,
                      bool new_lambda, Index nele_hess, Index* iRow,
                      Index* jCol, Number* values) override;
  virtual void finalize_solution(SolverReturn status,
                                 Index n, const Number* x, const Number* z_L, const Number* z_U,
                                 Index m, const Number* g, const Number* lambda,
                                 Number obj_value,
                                 const IpoptData* ip_data,
                                 IpoptCalculatedQuantities* ip_cq) override;
  static bool optimize(typename DSSQPObjective<T>::Vec& x, DSSQPObjective<T>& obj, scalarD tolG, sizeType maxIter, scalarD derivCheck, scalarD muInit, bool dense=false, int callback=5);
  static scalarD getInfty() {
    return 2e19;
  }
protected:
  void updateX(const Cold& x);
  bool _isDense,_timing;
  DSSQPObjective<T>& _obj;
  DSSQPObjective<scalarD>::DMat _gJacD,_fHessD;
  DSSQPObjective<scalarD>::SMat _gJacS,_fHessS;
  std::shared_ptr<Cold> _xl,_xu;
  std::shared_ptr<Cold> _gl,_gu;
  Cold _gVec,_fGrad;
  scalarD _E;
  Cold& _x;
};

PRJ_END

#endif
#endif
