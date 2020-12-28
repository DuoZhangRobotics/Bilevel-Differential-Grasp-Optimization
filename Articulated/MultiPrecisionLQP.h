#ifndef MULTI_PRECISION_LQP_H
#define MULTI_PRECISION_LQP_H

#include <Utils/ArticulatedBodyPragma.h>
#include <Utils/DebugGradient.h>
#include <Utils/SparseUtils.h>
#include <Utils/Options.h>

PRJ_BEGIN

//Solve the problem of following form:
//argmin c*w+w^THw/2
//  s.t. w_i>=0
//       \sum_i w_i<=1
//
//We solve this by minimizing:
//argmin f(w)=c*w+w^THw/2-\mu\sum_i\log(w_i)-\mu\log(1-\sum_i w_i)
//
//We then take Newton's step:
//H(w)\Delta(w)+G(w)=0
//We can simply use a line-search so that w+\Delta(w) does not violate the constraint
template <typename T>
struct SolveNewton
{
  DECL_MAP_TYPES_T
  typedef Eigen::Matrix<double,-1,-1> Matd;
public:
  template <typename RHS>
  static bool solveLU(const DMat& h,const RHS& g,RHS& dx,bool highPrec) {
    if(highPrec==1) {
      DMat LU;
      Coli P;
      if(!LUP(h,LU,P))
        return false;
      LUPSolve<RHS>(LU,P,dx=-g);
      return true;
    } else {
      return solveRefinement<Eigen::PartialPivLU<Matd>,RHS>(h,g,dx);
    }
  }
  template <typename RHS>
  static bool solveNewton(const DMat& h,const RHS& g,RHS& dx,bool highPrec) {
    if(highPrec==1) {
      dx=-g;
      DMat hTmp=h;
      LTDL(hTmp);
      for(sizeType i=0; i<hTmp.rows(); i++)
        if(hTmp(i,i)<=0)
          return false;
      LTDLSolve<RHS>(hTmp,dx);
      return true;
    } else {
      return solveRefinement<Eigen::LDLT<Matd>,RHS>(h,g,dx);
    }
  }
  template <typename FACTOR,typename RHS>
  static bool solveRefinement(const DMat& h,const RHS& g,RHS& dx)
  {
    //scale
    T maxG=0,maxH=0;
    for(sizeType i=0; i<g.size(); i++) {
      //g
      for(sizeType j=0; j<g.cols(); j++) {
        T gAbs=std::abs(g(i,j));
        if(gAbs>maxG)
          maxG=gAbs;
      }
      //h
      for(sizeType j=0; j<h.cols(); j++) {
        T hAbs=std::abs(h(i,j));
        if(hAbs>maxH)
          maxH=hAbs;
      }
    }
    if(maxG>0) {
      //factorize low-precision
      Matd hL=(h/maxH).unaryExpr([&](const T& in) {
        return (double)std::to_double(in);
      });
      FACTOR factor(hL);
      //solve low-precision
      dx=-g/maxG;
      dx=factor.solve(dx.unaryExpr([&](const T& in) {
        return (double)std::to_double(in);
      })).template cast<T>()*(maxG/maxH);
    } else dx.setZero(g.rows(),g.cols());
    return true;
  }
  template <typename RHS>
  static void inverse(const DMat& h,RHS& hi) {
    DMat rhs;
    LTDL(rhs=h);
    hi.setIdentity(h.rows(),h.cols());
    LTDLSolve<RHS>(rhs,hi);
  }
private:
  static bool LUP(const DMat& A,DMat& LU,Coli& P)
  {
    sizeType N=A.rows();
    LU=A;
    //#Unit permutation matrix, P[N] initialized with N
    P.resize(N+1);
    for(sizeType i=0; i<P.size(); i++)
      P[i]=i;
    //#factorize
    for(sizeType i=0; i<N; i++) {
      T maxA=0;
      sizeType imax=i;
      for(sizeType k=i; k<N; k++) {
        T absA=std::abs(LU(k,i));
        if(absA>maxA) {
          maxA=absA;
          imax=k;
        }
      }
      if(maxA<DELTAPrecision<T>::delta()) {
        //#failure, matrix is degenerate
        return false;
      }
      if(imax!=i) {
        //#pivoting P
        sizeType j=P[i];
        P[i]=P[imax];
        P[imax]=j;
        //#pivoting rows of A
        Vec ptr=LU.row(i);
        LU.row(i)=LU.row(imax);
        LU.row(imax)=ptr;
        //#counting pivots starting from N (for determinant)
        P[N]+=1;
      }
      for(sizeType j=i+1; j<N; j++) {
        LU(j,i)/=LU(i,i);
        for(sizeType k=i+1; k<N; k++)
          LU(j,k)-=LU(j,i)*LU(i,k);
      }
    }
    return true;
  }
  template <typename RHS>
  static void LUPSolve(const DMat& LU,const Coli& P,RHS& b)
  {
    sizeType N=LU.rows();
    RHS x=RHS::Zero(b.rows(),b.cols());
    for(sizeType i=0; i<N; i++) {
      x.row(i)=b.row(P[i]);
      for(sizeType k=0; k<i; k++)
        x.row(i)-=LU(i,k)*x.row(k);
    }
    for(sizeType n_i=0; n_i<N; n_i++) {
      sizeType i=N-1-n_i;
      for(sizeType k=i+1; k<N; k++)
        x.row(i)-=LU(i,k)*x.row(k);
      x.row(i)/=LU(i,i);
    }
    b=x;
  }
  static void LTDL(DMat& h)
  {
    //factorize
    T a;
    for(sizeType k=h.rows()-1,i,j; k>=0; k--) {
      i=k-1;
      while(i>=0) {
        a=h(k,i)/h(k,k);
        j=i;
        while(j>=0) {
          h(i,j)-=a*h(k,j);
          j=j-1;
        }
        h(k,i)=a;
        i=i-1;
      }
    }
  }
  template <typename RHS>
  static void LTDLSolve(const DMat& h,RHS& x)
  {
    //solve LT
    for(sizeType i=h.rows()-1,j; i>=0; i--) {
      j=i-1;
      while(j>=0) {
        x.row(j)-=h(i,j)*x.row(i);
        j=j-1;
      }
    }
    //solve D
    for(sizeType i=0; i<h.rows(); i++)
      x.row(i)/=h(i,i);
    //solve L
    for(sizeType i=0,j; i<h.rows(); i++) {
      j=i-1;
      while(j>=0) {
        x.row(i)-=h(i,j)*x.row(j);
        j=j-1;
      }
    }
  }
};
template <typename T>
struct SolveNewtonLP
{
  DECL_MAP_TYPES_T
  typedef Eigen::Matrix<double,-1,-1> Matd;
public:
  template <typename RHS>
  static bool solveNewton(T mu,const Vec& w,const RHS& g,RHS& dx,const Coli& foot) {
    DMat hi;
    inverse<DMat>(mu,w,hi,foot);
    dx=-hi*g;
    return true;
  }
  template <typename RHS>
  static void inverse(T mu,const Vec& w,RHS& hi,const Coli& foot) {
    hi.setZero(w.size(),w.size());
    for(sizeType f=0,off=0; f<foot.size(); off+=foot[f],f++) {
      Eigen::Block<const Vec> wSeg=w.block(off,0,foot[f],1);
      Vec wSegSqr=(wSeg.array()*wSeg.array()).matrix();
      T wSumVio=1-wSeg.sum();
      Eigen::Block<RHS> hiBlk=hi.block(off,off,foot[f],foot[f]);
      hiBlk=wSegSqr.asDiagonal();
      hiBlk-=wSegSqr*wSegSqr.transpose()/(wSegSqr.sum()+wSumVio*wSumVio);
    }
    hi/=mu;
  }
};
template <>
struct SolveNewton<double>
{
  typedef double T;
  DECL_MAP_TYPES_T
public:
  template <typename RHS>
  static bool solveLU(const DMat& h,const RHS& g,RHS& dx,bool) {
    Eigen::PartialPivLU<DMat> sol=h.partialPivLu();
    //if(sol.info()!=Eigen::Success)
    //  return false;
    dx=-sol.solve(g);
    return true;
  }
  template <typename RHS>
  static bool solveNewton(const DMat& h,const RHS& g,RHS& dx,bool) {
    Eigen::LDLT<DMat> sol=h.ldlt();
    if(sol.info()!=Eigen::Success)
      return false;
    dx=-sol.solve(g);
    return true;
  }
  template <typename RHS>
  static bool inverse(const DMat& h,RHS& hi) {
    Eigen::LDLT<DMat> sol=h.ldlt();
    if(sol.info()!=Eigen::Success)
      return false;
    hi=h.ldlt().solve(DMat::Identity(h.rows(),h.cols()));
    return true;
  }
};
template <typename T>
class MDP;
template <typename T>
struct MultiPrecisionLQPTraits;
template <typename T>
class MultiPrecisionLQP
{
public:
  DECL_MAP_TYPES_T
  friend struct MDP<T>;
  friend struct MultiPrecisionLQPTraits<T>;
  MultiPrecisionLQP(Options& ops);
  MultiPrecisionLQP(Options& ops,const Vec& c,const DMat& H,const Coli& foot);
  void reset(Options& ops);
  void reset(Options& ops,const Vec& c,const DMat& H,const Coli& foot);
  Vec solve(bool& succ,const Vec* initw=NULL) const;
  virtual Vec sampleValidW() const;
  void debugGradient();
  bool callback() const;
  bool highPrec() const;
  T tolGFinal() const;
  T muFinal() const;
  T tolAlpha() const;
  sizeType maxIter() const;
  void writeProb(const std::string& path) const;
  void readAndTestProb(const std::string& path);
  virtual T computeFGH(T mu,const Vec& w,Vec* g=NULL,DMat* h=NULL) const;
protected:
  virtual void initialGuess(Vec& w) const;
  virtual void printConstraints(const Vec& w) const;
  virtual void limitAlpha(T& alpha,const Vec& w,const Vec& dx) const;
  //data
  Coli _foot;
  DMat _H;
  Vec _c;
  //param
  T _muInit,_muDec,_muFinal,_margin;
  T _alphaInc,_alphaDec,_tolAlpha;
  T _tolG,_tolGFinal,_tolGFirstOrder;
  T _c1,_c2;    //Wolfe-Condition
  bool _callback,_highPrec;
  sizeType _maxIter;
};

PRJ_END

#endif
