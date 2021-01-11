#include "MultiPrecisionLQP.h"
#include <Utils/Utils.h>

PRJ_BEGIN

//parameters
template <>
struct MultiPrecisionLQPTraits<scalarD>
{
  static void registerOptions(Options& ops)
  {
    REGISTER_FLOAT_TYPE("muInit",MultiPrecisionLQP<scalarD>,scalarD,t._muInit)
    REGISTER_FLOAT_TYPE("muDec",MultiPrecisionLQP<scalarD>,scalarD,t._muDec)
    REGISTER_FLOAT_TYPE("muFinal",MultiPrecisionLQP<scalarD>,scalarD,t._muFinal)
    REGISTER_FLOAT_TYPE("margin",MultiPrecisionLQP<scalarD>,scalarD,t._margin)

    REGISTER_FLOAT_TYPE("alphaInc",MultiPrecisionLQP<scalarD>,scalarD,t._alphaInc)
    REGISTER_FLOAT_TYPE("alphaDec",MultiPrecisionLQP<scalarD>,scalarD,t._alphaDec)
    REGISTER_FLOAT_TYPE("tolAlpha",MultiPrecisionLQP<scalarD>,scalarD,t._tolAlpha)

    REGISTER_FLOAT_TYPE("tolG",MultiPrecisionLQP<scalarD>,scalarD,t._tolG)
    REGISTER_FLOAT_TYPE("tolGFinal",MultiPrecisionLQP<scalarD>,scalarD,t._tolGFinal)
    REGISTER_FLOAT_TYPE("tolGFirstOrder",MultiPrecisionLQP<scalarD>,scalarD,t._tolGFirstOrder)
    REGISTER_FLOAT_TYPE("c1",MultiPrecisionLQP<scalarD>,scalarD,t._c1)
    REGISTER_FLOAT_TYPE("c2",MultiPrecisionLQP<scalarD>,scalarD,t._c2)

    REGISTER_BOOL_TYPE("callback",MultiPrecisionLQP<scalarD>,bool,t._callback)
    REGISTER_BOOL_TYPE("forceSPD",MultiPrecisionLQP<scalarD>,bool,t._forceSPD)
    REGISTER_BOOL_TYPE("highPrec",MultiPrecisionLQP<scalarD>,bool,t._highPrec)
    REGISTER_INT_TYPE("maxIter",MultiPrecisionLQP<scalarD>,sizeType,t._maxIter)
  }
  static void initOptions(MultiPrecisionLQP<scalarD>& sol)
  {
    sol._muInit=1e-1f;
    sol._muDec=0.1f;
    sol._muFinal=1e-3f;
    sol._margin=1e-10f;

    sol._alphaInc=1.1f;
    sol._alphaDec=0.5f;
    sol._tolAlpha=1e-20f;

    sol._tolG=1.0f;
    sol._tolGFinal=1e-8f;
    sol._tolGFirstOrder=1e10f;
    sol._c1=0.0f;
    sol._c2=1.0f;

    sol._callback=true;
    sol._forceSPD=true;
    sol._highPrec=true;
    sol._maxIter=10000;
  }
};
template <>
struct MultiPrecisionLQPTraits<__float128>
{
  static void registerOptions(Options& ops)
  {
    REGISTER_FLOAT128_TYPE("muInit",MultiPrecisionLQP<__float128>,__float128,t._muInit)
    REGISTER_FLOAT128_TYPE("muDec",MultiPrecisionLQP<__float128>,__float128,t._muDec)
    REGISTER_FLOAT128_TYPE("muFinal",MultiPrecisionLQP<__float128>,__float128,t._muFinal)
    REGISTER_FLOAT128_TYPE("margin",MultiPrecisionLQP<__float128>,__float128,t._margin)

    REGISTER_FLOAT128_TYPE("alphaInc",MultiPrecisionLQP<__float128>,__float128,t._alphaInc)
    REGISTER_FLOAT128_TYPE("alphaDec",MultiPrecisionLQP<__float128>,__float128,t._alphaDec)
    REGISTER_FLOAT128_TYPE("tolAlpha",MultiPrecisionLQP<__float128>,__float128,t._tolAlpha)

    REGISTER_FLOAT128_TYPE("tolG",MultiPrecisionLQP<__float128>,__float128,t._tolG)
    REGISTER_FLOAT128_TYPE("tolGFinal",MultiPrecisionLQP<__float128>,__float128,t._tolGFinal)
    REGISTER_FLOAT128_TYPE("tolGFirstOrder",MultiPrecisionLQP<__float128>,__float128,t._tolGFirstOrder)
    REGISTER_FLOAT128_TYPE("c1",MultiPrecisionLQP<__float128>,__float128,t._c1)
    REGISTER_FLOAT128_TYPE("c2",MultiPrecisionLQP<__float128>,__float128,t._c2)

    REGISTER_BOOL_TYPE("callback",MultiPrecisionLQP<__float128>,bool,t._callback)
    REGISTER_BOOL_TYPE("forceSPD",MultiPrecisionLQP<__float128>,bool,t._forceSPD)
    REGISTER_BOOL_TYPE("highPrec",MultiPrecisionLQP<__float128>,bool,t._highPrec)
    REGISTER_INT_TYPE("maxIter",MultiPrecisionLQP<__float128>,sizeType,t._maxIter)
  }
  static void initOptions(MultiPrecisionLQP<__float128>& sol)
  {
    sol._muInit=1e-1f;
    sol._muDec=0.1f;
    sol._muFinal=1e-5f;
    sol._margin=1e-20f;

    sol._alphaInc=1.1f;
    sol._alphaDec=0.5f;
    sol._tolAlpha=1e-20f;

    sol._tolG=1.0f;
    sol._tolGFinal=1e-20f;
    sol._tolGFirstOrder=1e20f;
    sol._c1=0.0f;
    sol._c2=1.0f;

    sol._callback=true;
    sol._forceSPD=true;
    sol._highPrec=false;
    sol._maxIter=10000;
  }
};
template <>
struct MultiPrecisionLQPTraits<mpfr::mpreal>
{
  static void registerOptions(Options& ops)
  {
    REGISTER_MPFR_TYPE("muInit",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._muInit)
    REGISTER_MPFR_TYPE("muDec",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._muDec)
    REGISTER_MPFR_TYPE("muFinal",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._muFinal)
    REGISTER_MPFR_TYPE("margin",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._margin)

    REGISTER_MPFR_TYPE("alphaInc",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._alphaInc)
    REGISTER_MPFR_TYPE("alphaDec",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._alphaDec)
    REGISTER_MPFR_TYPE("tolAlpha",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._tolAlpha)

    REGISTER_MPFR_TYPE("tolG",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._tolG)
    REGISTER_MPFR_TYPE("tolGFinal",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._tolGFinal)
    REGISTER_MPFR_TYPE("tolGFirstOrder",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._tolGFirstOrder)
    REGISTER_MPFR_TYPE("c1",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._c1)
    REGISTER_MPFR_TYPE("c2",MultiPrecisionLQP<mpfr::mpreal>,mpfr::mpreal,t._c2)

    REGISTER_BOOL_TYPE("callback",MultiPrecisionLQP<mpfr::mpreal>,bool,t._callback)
    REGISTER_BOOL_TYPE("forceSPD",MultiPrecisionLQP<mpfr::mpreal>,bool,t._forceSPD)
    REGISTER_BOOL_TYPE("highPrec",MultiPrecisionLQP<mpfr::mpreal>,bool,t._highPrec)
    REGISTER_INT_TYPE("maxIter",MultiPrecisionLQP<mpfr::mpreal>,sizeType,t._maxIter)
  }
  static void initOptions(MultiPrecisionLQP<mpfr::mpreal>& sol)
  {
    sol._muInit=1e-1f;
    sol._muDec=0.1f;
    sol._muFinal=1e-5f;
    sol._margin=1e-20f;

    sol._alphaInc=1.1f;
    sol._alphaDec=0.5f;
    sol._tolAlpha=DELTAPrecision<mpfr::mpreal>::delta();

    sol._tolG=1.0f;
    sol._tolGFinal=std::max<mpfr::mpreal>(mpfr::mpreal("1e-100"),DELTAPrecision<mpfr::mpreal>::delta());
    sol._tolGFirstOrder=mpfr::mpreal("1e100");
    sol._c1=0.0f;
    sol._c2=1.0f;

    sol._callback=true;
    sol._forceSPD=true;
    sol._highPrec=false;
    sol._maxIter=10000;
  }
};
//MultiPrecisionLQP
template <typename T>
MultiPrecisionLQP<T>::MultiPrecisionLQP(Options& ops)
{
  if(!ops.hasType<MultiPrecisionLQP<T>>())
    MultiPrecisionLQPTraits<T>::registerOptions(ops);
  reset(ops);
}
template <typename T>
MultiPrecisionLQP<T>::MultiPrecisionLQP(Options& ops,const Vec& c,const DMat& H,const Coli& foot)
{
  if(!ops.hasType<MultiPrecisionLQP<T>>())
    MultiPrecisionLQPTraits<T>::registerOptions(ops);
  reset(ops,c,H,foot);
}
template <typename T>
void MultiPrecisionLQP<T>::reset(Options& ops)
{
  MultiPrecisionLQPTraits<T>::initOptions(*this);
  ops.setOptions(this);
}
template <typename T>
void MultiPrecisionLQP<T>::reset(Options& ops,const Vec& c,const DMat& H,const Coli& foot)
{
  reset(ops);
  _foot=foot;
  _H=H;
  _c=c;
  ASSERT(_foot.sum()==_c.size())
}
template <typename T>
typename MultiPrecisionLQP<T>::Vec MultiPrecisionLQP<T>::solve(bool& succ,const Vec* initw) const
{
  succ=false;
  DMat h,hi;
  Vec dx,g,g2,w0,w;
  T alphaInit,f,vio,dTg,dTg2,alpha=1,f2=0,mu=_muInit;
  //feasible initial guess
  if(initw) {
    w=*initw;
    mu=_muFinal;    //this is a heuristic, if user provides init value, we set mu to the smallest possible value.
  } else initialGuess(w);
  if(_callback) {
    std::cout << "QPProblem:" << std::endl;
    printConstraints(w);
  }
  //main loop
  for(sizeType it=0;; it++) {
    f=computeFGH(mu,w,&g,&h,_forceSPD);
    vio=std::sqrt(g.squaredNorm());
    if(mu==_muFinal && vio<_tolGFinal) {
      if(_callback) {
        std::cout << "Solved with tolG(" << vio << ")<=" << _tolG << std::endl;
        printConstraints(w);
      }
      succ=true;
      break;
    }
    //newton
    if(!solveDx(vio,mu,w,dx,h,g)) {
      std::cout << "Iteration " << it << ": solver failed" << std::endl;
      return w;
    }
    //limit alpha
    limitAlpha(alpha,w,dx);
    //line-search
    alphaInit=alpha;
    {
      w0=w;
      dTg=dx.dot(g);
      while(alpha>_tolAlpha) {
        w=w0+dx*alpha;
        f2=computeFGH(mu,w,&g2,NULL);
        dTg2=dx.dot(g2);
        if(std::isfinite(f2) && f2<f+_c1*alpha*dTg && -dTg2<=-_c2*dTg) {
          //smooth-change of alpha
          alpha=std::min(alpha*_alphaInc,T(1));
          break;
        } else {
          w=w0;
          alpha*=_alphaDec;
        }
      }
    }
    //termination condition
    if(_callback)
      std::cout << "Iteration " << it << ": mu(" << mu << ") f(" << f << ") f2(" << f2 << ") tolG(" << vio << ")" << " alphaInit(" << alphaInit << ") alphaFinal(" << alpha << ")" << std::endl;
    if(alpha<=_tolAlpha) {
      if(_callback)
        std::cout << "Stopped with alpha(" << alpha << ")<=" << _tolAlpha << " tolG(" << vio << ")" << std::endl;
      if(mu>_muFinal)
        mu=std::max<T>(mu*_muDec,_muFinal);
      else {
        if(_callback)
          printConstraints(w);
        break;
      }
    }
    if(it>=_maxIter) {
      if(_callback) {
        std::cout << "Stopped with it(" << it << ")>=" << _maxIter << " tolG(" << _tolG << ")" << std::endl;
        printConstraints(w);
      }
      break;
    }
    if(vio<_tolG)
      mu=std::max<T>(mu*_muDec,_muFinal);
  }
  return w;
}
template <typename T>
typename MultiPrecisionLQP<T>::Vec MultiPrecisionLQP<T>::sampleValidW() const
{
  Vec w;
  while(true) {
    bool valid=true;
    w=Vec::Random(_c.size())*0.5f+Vec::Constant(_c.size(),0.5f);
    for(sizeType f=0,off=0; f<_foot.size(); off+=_foot[f],f++) {
      Eigen::Block<Vec> wSeg=w.block(off,0,_foot[f],1);
      for(sizeType r=0; r<wSeg.size(); r++)
        if(wSeg(r,0)<0)
          valid=false;
      while(wSeg.sum()>1)
        wSeg/=2;
    }
    if(valid)
      break;
  }
  return w;
}
template <typename T>
void MultiPrecisionLQP<T>::debugGradient()
{
  DEFINE_NUMERIC_DELTA_T(T)
  {
    INFO("-------------------------------------------------------------DebugSolveNewton")
    sizeType N=_c.size();
    Vec g=Vec::Random(N),x,w,ctmp;
    DMat h=DMat::Random(N,N),Htmp;
    SolveNewton<T>::solveLU(h,g,x,false);
    DEBUG_GRADIENT("SolveLU-Low",std::sqrt(g.squaredNorm()),std::sqrt((h*x+g).squaredNorm()))
    SolveNewton<T>::solveLU(h,g,x,true);
    DEBUG_GRADIENT("SolveLU-High",std::sqrt(g.squaredNorm()),std::sqrt((h*x+g).squaredNorm()))
    h=(h*h.transpose()).eval();
    SolveNewton<T>::solveNewton(h,g,x,false);
    DEBUG_GRADIENT("solveNewton-Low",std::sqrt(g.squaredNorm()),std::sqrt((h*x+g).squaredNorm()))
    SolveNewton<T>::solveNewton(h,g,x,true);
    DEBUG_GRADIENT("solveNewton-High",std::sqrt(g.squaredNorm()),std::sqrt((h*x+g).squaredNorm()))
    Htmp=_H;
    ctmp=_c;
    _H.resize(0,0);
    _c.setRandom(N);
    w=sampleValidW();
    computeFGH(_muFinal,w,NULL,&h);
    SolveNewtonLP<T>::solveNewton(_muFinal,w,g,x,_foot);
    DEBUG_GRADIENT("solveNewton-LP",std::sqrt(g.squaredNorm()),std::sqrt((h*x+g).squaredNorm()))
    _H=Htmp;
    _c=ctmp;
  }

  while(true) {
    DMat h,ih;
    Vec w=sampleValidW(),g,g2,dw=Vec::Random(_c.size());
    T f=computeFGH(_muInit,w,&g,&h,false);
    T f2=computeFGH(_muInit,w+dw*DELTA,&g2);
    if(!std::isfinite(f))
      continue;
    INFO("-------------------------------------------------------------DebugLQPGradient")
    DEBUG_GRADIENT("MultiPrecisionLQP-DF",g.dot(dw),g.dot(dw)-(f2-f)/DELTA)
    DEBUG_GRADIENT("MultiPrecisionLQP-DG",std::sqrt((h*dw).squaredNorm()),std::sqrt((h*dw-(g2-g)/DELTA).squaredNorm()))
    break;
  }
}
template <typename T>
bool MultiPrecisionLQP<T>::callback() const
{
  return _callback;
}
template <typename T>
bool MultiPrecisionLQP<T>::highPrec() const
{
  return _highPrec;
}
template <typename T>
T MultiPrecisionLQP<T>::tolGFinal() const
{
  return _tolGFinal;
}
template <typename T>
T MultiPrecisionLQP<T>::muFinal() const
{
  return _muFinal;
}
template <typename T>
T MultiPrecisionLQP<T>::tolAlpha() const
{
  return _tolAlpha;
}
template <typename T>
sizeType MultiPrecisionLQP<T>::maxIter() const
{
  return _maxIter;
}
template <typename T>
void MultiPrecisionLQP<T>::writeProb(const std::string& path) const
{
  std::ofstream os(path,std::ios::binary);
  writeBinaryData(_foot,os);
  Matd H=_H.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Cold c=_c.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(H,os);
  writeBinaryData(c,os);
}
template <typename T>
void MultiPrecisionLQP<T>::readAndTestProb(const std::string& path)
{
  bool succ;
  if(!exists(path))
    return;
  std::ifstream is(path,std::ios::binary);
  readBinaryData(_foot,is);
  Matd H;
  Cold c;
  readBinaryData(H,is);
  readBinaryData(c,is);
  _H=H.template cast<T>();
  _c=c.template cast<T>();
  solve(succ);
}
template <typename T>
T MultiPrecisionLQP<T>::computeFGH(T mu,const Vec& w,Vec* g,DMat* h,bool) const
{
  Vec Hw;
  if(_H.size()==0)
    Hw.setZero(_c.size());
  else Hw=_H*w;
  T obj=(_c+Hw/2).dot(w),wSumVio;
  if(g)
    *g=_c+Hw;
  if(h) {
    if(_H.size()==0)
      h->setZero(_c.size(),_c.size());
    else *h=_H;
  }
  //g
  for(sizeType f=0,off=0; f<_foot.size(); off+=_foot[f],f++) {
    Eigen::Block<const Vec> wSeg=w.block(off,0,_foot[f],1);
    wSumVio=1-wSeg.sum();
    obj-=std::log(wSumVio)*mu;
    for(sizeType i=off; i<off+_foot[f]; i++)
      obj-=std::log(w[i])*mu;
    if(g) {
      g->segment(off,_foot[f]).array()+=mu/wSumVio;
      g->segment(off,_foot[f]).array()-=mu/wSeg.array();
    }
    if(h) {
      h->block(off,off,_foot[f],_foot[f]).array()+=mu/wSumVio/wSumVio;
      h->block(off,off,_foot[f],_foot[f]).diagonal().array()+=mu/wSeg.array()/wSeg.array();
    }
  }
  return obj;
}
//helper
template <typename T>
void MultiPrecisionLQP<T>::initialGuess(Vec& w) const
{
  w.resize(_c.size());
  for(sizeType f=0,off=0; f<_foot.size(); off+=_foot[f],f++) {
    Eigen::Block<Vec> wSeg=w.block(off,0,_foot[f],1);
    wSeg.setConstant(1/(T)(wSeg.size()+1));
  }
}
template <typename T>
void MultiPrecisionLQP<T>::printConstraints(const Vec& w) const
{
  for(sizeType f=0,off=0; f<_foot.size(); off+=_foot[f],f++) {
    Eigen::Block<const Vec> wSeg=w.block(off,0,_foot[f],1);
    T wSumVio=1-wSeg.sum();
    std::cout << "Foot" << f << ": wsumVio=" << wSumVio << " ";
    for(sizeType i=0; i<wSeg.size(); i++)
      std::cout << "w" << i << "Vio=" << wSeg(i,0) << " ";
    std::cout << std::endl;
  }
}
template <typename T>
void MultiPrecisionLQP<T>::limitAlpha(T& alpha,const Vec& w,const Vec& dx) const
{
  T dxSum;
  for(sizeType i=0; i<_c.size(); i++)
    //make sure: w[i]+dx[i]*alpha>=margin
    if(dx[i]<0)
      alpha=std::min(alpha,(w[i]-_margin)/-dx[i]);
  //make sure: 1-sum(w)-sum(dx)*alpha>=margin
  for(sizeType f=0,off=0; f<_foot.size(); off+=_foot[f],f++) {
    Eigen::Block<const Vec> dxSeg=dx.block(off,0,_foot[f],1);
    Eigen::Block<const Vec> wSeg=w.block(off,0,_foot[f],1);
    dxSum=dxSeg.sum();
    if(dxSum>0)
      alpha=std::min(alpha,(1-wSeg.sum()-_margin)/dxSum);
  }
}
template <typename T>
bool MultiPrecisionLQP<T>::solveDx(T vio,T mu,const Vec& w,Vec& dx,const DMat& h,const Vec& g) const
{
  bool succ=true;
  if(vio>_tolGFirstOrder)
    dx=-g;
  else {
    if(_H.size()==0)
      succ=SolveNewtonLP<T>::solveNewton(mu,w,g,dx,_foot);
    else succ=SolveNewton<T>::solveNewton(h,g,dx,_highPrec);
  }
  return succ;
}
//instance
template class MultiPrecisionLQP<double>;
//#ifdef ALL_TYPES  //we have to comment this line out because we need ultra-high precision for LPForceSequence
template class MultiPrecisionLQP<__float128>;
template class MultiPrecisionLQP<mpfr::mpreal>;
//#endif
PRJ_END
