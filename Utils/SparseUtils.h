#ifndef SPARSE_UTILS_H
#define SPARSE_UTILS_H

#include <Eigen/Sparse>
#include <CommonFile/IO.h>
#include "ParallelVector.h"
#include "Scalar.h"

PRJ_BEGIN

template <typename T>
struct SparseTraits
{
  typedef Eigen::SparseMatrix<T,0,sizeType> SMat;
  typedef Eigen::Matrix<T,-1,-1> DMat;
  typedef Eigen::Matrix<T,-1,1> Vec;
  typedef Eigen::Triplet<T,sizeType> STrip;
  typedef ParallelVector<STrip> STrips;
};
//sparse matrix building
//dense
template <typename MAT,typename Derived>
static void addBlock(MAT& H,sizeType r,sizeType c,const Eigen::MatrixBase<Derived>& coef)
{
  H.block(r,c,coef.rows(),coef.cols())+=coef;
}
template <typename MAT>
static void addBlock(MAT& H,sizeType r,sizeType c,typename MAT::Scalar coef)
{
  H(r,c)+=coef;
}
template <typename MAT,typename T>
static void addBlockId(MAT& H,sizeType r,sizeType c,sizeType nr,T coefId)
{
  H.block(r,c,nr,nr).diagonal().array()+=coefId;
}
//sparse
template <typename Derived>
static void addBlock(ParallelVector<Eigen::Triplet<typename Derived::Scalar,sizeType>>& H,sizeType r,sizeType c,const Eigen::MatrixBase<Derived>& coef)
{
  sizeType nrR=coef.rows();
  sizeType nrC=coef.cols();
  for(sizeType i=0; i<nrR; i++)
    for(sizeType j=0; j<nrC; j++)
      H.push_back(Eigen::Triplet<typename Derived::Scalar,sizeType>(r+i,c+j,coef(i,j)));
}
template <typename T>
static void addBlock(ParallelVector<Eigen::Triplet<T,sizeType>>& H,sizeType r,sizeType c,T coef)
{
  H.push_back(Eigen::Triplet<T,sizeType>(r,c,coef));
}
template <typename Derived>
static void addBlockI(ParallelVector<Eigen::Triplet<typename Derived::Scalar,sizeType>>& H,sizeType r,sizeType c,const Eigen::MatrixBase<Derived>& coefI)
{
  sizeType nr=coefI.size();
  for(sizeType i=0; i<nr; i++)
    H.push_back(Eigen::Triplet<typename Derived::Scalar,sizeType>(r+i,c+i,coefI[i]));
}
template <typename T>
static void addBlockId(ParallelVector<Eigen::Triplet<T,sizeType>>& H,sizeType r,sizeType c,sizeType nr,T coefId)
{
  for(sizeType i=0; i<nr; i++)
    H.push_back(Eigen::Triplet<T,sizeType>(r+i,c+i,coefId));
}
template <typename T,int O,typename I>
static void addBlock(ParallelVector<Eigen::Triplet<T,sizeType>>& H,sizeType r,sizeType c,const Eigen::SparseMatrix<T,O,I>& coef)
{
  for(sizeType k=0; k<coef.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(coef,k); it; ++it)
      H.push_back(Eigen::Triplet<T,sizeType>(r+it.row(),c+it.col(),it.value()));
}
//sparseVec
template <typename Derived>
static void addBlock(Eigen::SparseVector<typename Derived::Scalar,0,sizeType>& H,sizeType r,const Eigen::MatrixBase<Derived>& coef)
{
  for(sizeType i=0; i<coef.rows(); i++)
    H.coeffRef(r+i)+=coef[i];
}
//maxAbs
template <typename T,int O,typename I>
T absMax(const Eigen::SparseMatrix<T,O,I>& h)
{
  T ret=0;
  for(sizeType k=0; k<h.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(h,k); it; ++it)
      ret=std::max(ret,std::abs(it.value()));
  return ret;
}
template <typename T,int O,typename I>
T absMaxRel(const Eigen::SparseMatrix<T,O,I>& h,const Eigen::SparseMatrix<T,O,I>& hRef,bool detail)
{
  sizeType row=-1,col=-1;
  T ret=0,num=0,denom=0;
  T EPS=ScalarUtil<T>::scalar_max();
  //check against h
  for(sizeType k=0; k<h.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(h,k); it; ++it) {
      T err=std::abs(it.value()-hRef.coeff(it.row(),it.col()));
      T ref=std::abs(hRef.coeff(it.row(),it.col()));
      T rel=err/std::max<T>(ref,EPS);
      if(err<EPS)
        continue;
      if(rel>ret) {
        num=err;
        denom=ref;
        row=it.row();
        col=it.col();
      }
      ret=std::max(ret,rel);
    }
  //check against hRef
  for(sizeType k=0; k<hRef.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(hRef,k); it; ++it) {
      T err=std::abs(h.coeff(it.row(),it.col())-it.value());
      T ref=std::abs(it.value());
      T rel=err/std::max<T>(ref,EPS);
      if(err<EPS)
        continue;
      if(rel>ret) {
        num=err;
        denom=ref;
        row=it.row();
        col=it.col();
      }
      ret=std::max(ret,rel);
    }
  if(detail) {
    INFOV("(%d,%d) Num: %f Denom: %f",row,col,num,denom)
  }
  return ret;
}
//build KKT matrix
template <typename MAT>
MAT buildKKT(const MAT& h,const MAT& a,typename MAT::Scalar shift)
{
  MAT kkt;
  kkt.setZero(h.rows()+a.rows(),h.rows()+a.rows());
  kkt.block(0,0,h.rows(),h.rows())=h;
  kkt.block(0,0,h.rows(),h.rows()).diagonal().array()+=shift;
  kkt.block(h.rows(),0,a.rows(),a.cols())=a;
  kkt.block(0,h.rows(),a.cols(),a.rows())=a.transpose();
  return kkt;
}
template <typename T,int O,typename I>
typename Eigen::SparseMatrix<T,O,I> buildKKT(const Eigen::SparseMatrix<T,O,I>& h,const Eigen::SparseMatrix<T,O,I>& a,T shift,bool at=false,const Eigen::Matrix<T,-1,1>* shiftOffDiag=NULL)
{
  typedef Eigen::SparseMatrix<T,O,I> SMat;
  SMat kkt;
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(sizeType k=0; k<h.outerSize(); ++k)
    for(typename SMat::InnerIterator it(h,k); it; ++it)
      trips.push_back(typename Eigen::Triplet<T,I>(it.row(),it.col(),it.value()));
  if(shift!=0)
    for(sizeType k=0; k<h.rows(); ++k)
      trips.push_back(typename Eigen::Triplet<T,I>(k,k,shift));
  for(sizeType k=0; k<a.outerSize(); ++k)
    for(typename SMat::InnerIterator it(a,k); it; ++it)
      if(at) {
        trips.push_back(typename Eigen::Triplet<T,I>(h.rows()+it.col(),it.row(),it.value()));
        trips.push_back(typename Eigen::Triplet<T,I>(it.row(),h.rows()+it.col(),it.value()));
      } else {
        trips.push_back(typename Eigen::Triplet<T,I>(h.rows()+it.row(),it.col(),it.value()));
        trips.push_back(typename Eigen::Triplet<T,I>(it.col(),h.rows()+it.row(),it.value()));
      }
  sizeType dualDim=at?a.cols():a.rows();
  if(shiftOffDiag>0) {
    ASSERT_MSG(shiftOffDiag->size()==dualDim,"BuiltKKT has incorrect shiftOffDiag size")
    for(sizeType i=0; i<dualDim; i++)
      trips.push_back(typename Eigen::Triplet<T,I>(h.rows()+i,h.cols()+i,shiftOffDiag->coeff(i)));
  }
  kkt.resize(h.rows()+dualDim,h.rows()+dualDim);
  kkt.setFromTriplets(trips.begin(),trips.end());
  return kkt;
}
//kronecker-product
template <typename T,int O,typename I>
typename Eigen::SparseMatrix<T,O,I> kronecker(const Eigen::SparseMatrix<T,O,I>& h,sizeType n)
{
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(sizeType k=0; k<h.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(h,k); it; ++it)
      for(sizeType d=0; d<n; d++)
        trips.push_back(typename Eigen::Triplet<T,I>(it.row()*n+d,it.col()*n+d,it.value()));

  Eigen::SparseMatrix<T,O,I> ret;
  ret.resize(h.rows()*n,h.cols()*n);
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
//concat-diag
template <typename T,int O,typename I>
typename Eigen::SparseMatrix<T,O,I> concatDiag(const Eigen::SparseMatrix<T,O,I>& a,const Eigen::SparseMatrix<T,O,I>& b)
{
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(sizeType k=0; k<a.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(a,k); it; ++it)
      trips.push_back(typename Eigen::Triplet<T,I>(it.row(),it.col(),it.value()));
  for(sizeType k=0; k<b.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(b,k); it; ++it)
      trips.push_back(typename Eigen::Triplet<T,I>(a.rows()+it.row(),a.cols()+it.col(),it.value()));

  typename Eigen::SparseMatrix<T,O,I> ret;
  ret.resize(a.rows()+b.rows(),a.cols()+b.cols());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
//sparseIO
template <typename T,int O,typename I,typename MT>
Eigen::SparseMatrix<T,O,I> toSparse(const MT& m,T eps=0)
{
  Eigen::SparseMatrix<T,O,I> ret;
  ParallelVector<Eigen::Triplet<T,I>> trips;
  for(sizeType r=0; r<m.rows(); r++)
    for(sizeType c=0; c<m.cols(); c++)
      if(std::abs(m(r,c))>eps)
        trips.push_back(Eigen::Triplet<T,I>(r,c,m(r,c)));
  ret.resize(m.rows(),m.cols());
  ret.setFromTriplets(trips.begin(),trips.end());
  return ret;
}
template <typename T,int O,typename I>
Eigen::SparseMatrix<T,O,I> concat(const Eigen::SparseMatrix<T,O,I>& A,const Eigen::SparseMatrix<T,O,I>& B)
{
  ASSERT(A.cols()==B.cols())
  Eigen::SparseMatrix<T,O,I> M(A.rows()+B.rows(),A.cols());
  M.reserve(A.nonZeros()+B.nonZeros());
  for(sizeType c=0; c<A.cols(); ++c) {
    M.startVec(c);
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator itA(A,c); itA; ++itA)
      M.insertBack(itA.row(),c)=itA.value();
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator itB(B,c); itB; ++itB)
      M.insertBack(itB.row()+A.rows(),c)=itB.value();
  }
  M.finalize();
  return M;
}
template <typename T,int O,typename I>
bool readBinaryData(Eigen::SparseMatrix<T,O,I>& m,std::istream& is,IOData* dat=NULL)
{
  sizeType rows,cols;
  std::vector<I> r;
  std::vector<I> c;
  std::vector<scalarD> v;
  readBinaryData(rows,is);
  readBinaryData(cols,is);
  readBinaryData(r,is);
  readBinaryData(c,is);
  readBinaryData(v,is);
  m.resize(rows,cols);
  m.reserve(v.size());
  for(sizeType ci=0; ci<m.cols(); ++ci) {
    m.startVec(ci);
    for(sizeType off=c[ci]; off<c[ci+1]; off++)
      m.insertBack(r[off],ci)=v[off];
  }
  m.finalize();
  return is.good();
}
template <typename T,int O,typename I>
bool writeBinaryData(const Eigen::SparseMatrix<T,O,I>& m,std::ostream& os,IOData* dat=NULL)
{
  std::vector<I> r(m.nonZeros());
  std::vector<I> c(m.cols()+1);
  std::vector<scalarD> v(m.nonZeros());
  for(sizeType k=0,offr=0; k<m.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<T,O,I>::InnerIterator it(m,k); it; ++it) {
      v[offr]=std::to_double(it.value());
      r[offr++]=it.row();
      c[k+1]++;
    }
  for(sizeType k=0; k<m.outerSize(); ++k)
    c[k+1]+=c[k];
  writeBinaryData(m.rows(),os);
  writeBinaryData(m.cols(),os);
  writeBinaryData(r,os);
  writeBinaryData(c,os);
  writeBinaryData(v,os);
  return os.good();
}
//svec indices
static sizeType IJToTri(sizeType N,sizeType r,sizeType c)
{
  if(r>c)std::swap(r,c);
  return N*(N-1)/2-(N-r)*(N-r-1)/2+c;
}
static void TriToIJ(sizeType N,sizeType k,sizeType& i,sizeType& j)
{
  i=N-1-sizeType(std::sqrt(-8*k+4*N*(N+1)-7)/2.0-0.5);
  j=k+i-N*(N+1)/2+(N-i)*(N+1-i)/2;
}
static sizeType NSVecToN(sizeType N)
{
  return std::sqrt(N*2);
}
static void debugTriToIJ(sizeType N)
{
  sizeType i,j,off=0;
  for(sizeType r=0; r<N; r++)
    for(sizeType c=r; c<N; c++,off++) {
      TriToIJ(N,off,i,j);
      ASSERT(i==r && j==c && IJToTri(N,i,j)==off)
    }
}
//svec
template <typename T>
static typename Eigen::Matrix<T,-1,1> toSVec(const Eigen::Matrix<T,-1,-1>& m)
{
  typedef Eigen::Matrix<T,-1,1> Vec;
  //typedef Eigen::Matrix<T,-1,-1> DMat;
  const T sqrt2=std::sqrt((T)2);
  Vec ret=Vec::Zero(m.rows()*(m.rows()+1)/2);
  for(sizeType r=0,off=0; r<m.rows(); r++)
    for(sizeType c=r; c<m.cols(); c++,off++)
      ret[off]=r==c?m(r,c):m(r,c)*sqrt2;
  return ret;
}
template <typename T>
static typename Eigen::Matrix<T,-1,-1> fromSVec(const Eigen::Matrix<T,-1,1>& s)
{
  //typedef Eigen::Matrix<T,-1,1> Vec;
  typedef Eigen::Matrix<T,-1,-1> DMat;
  const T sqrt2=std::sqrt((T)2);
  sizeType rows=NSVecToN(s.size());
  DMat ret=DMat::Zero(rows,rows);
  for(sizeType r=0,off=0; r<ret.rows(); r++)
    for(sizeType c=r; c<ret.cols(); c++,off++)
      if(r==c)
        ret(r,c)=s[off];
      else ret(r,c)=ret(c,r)=s[off]/sqrt2;
  return ret;
}
//visualize
template <typename S,int O,typename I>
static void writeSMatVTK(const typename Eigen::SparseMatrix<S,O,I>& coef,const std::string& path,const Coli* blkR,const Coli* blkC)
{
  VTKWriter<scalar> os("SMat",path,true);
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>> tss;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> lss;
  //element
  for(sizeType k=0; k<coef.outerSize(); ++k)
    for(typename Eigen::SparseMatrix<S,O,I>::InnerIterator it(coef,k); it; ++it) {
      sizeType off=(sizeType)vss.size();
      vss.push_back(Vec3(it.col()  ,coef.rows()-it.row()-1,1));
      vss.push_back(Vec3(it.col()+1,coef.rows()-it.row()-1,1));
      vss.push_back(Vec3(it.col()+1,coef.rows()-it.row()  ,1));
      vss.push_back(Vec3(it.col()  ,coef.rows()-it.row()  ,1));
      tss.push_back(Vec3i(off+0,off+1,off+2));
      tss.push_back(Vec3i(off+0,off+2,off+3));
    }
  //line
  {
    sizeType off=(sizeType)vss.size();
    vss.push_back(Vec3(0          ,0          ,1));
    vss.push_back(Vec3(coef.cols(),0          ,1));
    vss.push_back(Vec3(coef.cols(),coef.rows(),1));
    vss.push_back(Vec3(0          ,coef.rows(),1));
    lss.push_back(Vec2i(off+0,off+1));
    lss.push_back(Vec2i(off+1,off+2));
    lss.push_back(Vec2i(off+2,off+3));
    lss.push_back(Vec2i(off+3,off+0));
  }
  //blk lines
  if(blkR)
    for(sizeType i=0,Y=coef.rows(); i<blkR->size(); Y-=blkR->coeff(i),i++)
      if(Y>0 && Y<coef.rows()) {
        sizeType off=(sizeType)vss.size();
        vss.push_back(Vec3(0          ,Y,1));
        vss.push_back(Vec3(coef.cols(),Y,1));
        lss.push_back(Vec2i(off+0,off+1));
      }
  if(blkC)
    for(sizeType i=0,X=0; i<blkC->size(); X+=blkC->coeff(i),i++)
      if(X>0 && X<coef.cols()) {
        sizeType off=(sizeType)vss.size();
        vss.push_back(Vec3(X,0          ,1));
        vss.push_back(Vec3(X,coef.rows(),1));
        lss.push_back(Vec2i(off+0,off+1));
      }
  //feed
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(tss.begin(),tss.end(),VTKWriter<scalar>::TRIANGLE);
  os.appendCells(lss.begin(),lss.end(),VTKWriter<scalar>::LINE);
}
//3x3 multiply Left/Right flatten
template <typename T>
Eigen::SparseMatrix<T,0,sizeType> mat3x3MulRight(const Eigen::Matrix<T,3,3>& m)
{
  typedef ParallelVector<Eigen::Triplet<T,sizeType> > TRIPS;
  Eigen::SparseMatrix<T,0,sizeType> M;
  M.resize(9,9);
  TRIPS trips;
  for(sizeType r=0; r<3; r++)
    for(sizeType c=0; c<3; c++)
      for(sizeType k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,sizeType>(r+c*3,r+k*3,m(k,c)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}
template <typename T>
Eigen::SparseMatrix<T,0,sizeType> mat3x3MulTRight(const Eigen::Matrix<T,3,3>& m)
{
  typedef ParallelVector<Eigen::Triplet<T,sizeType> > TRIPS;
  Eigen::SparseMatrix<T,0,sizeType> M;
  M.resize(9,9);
  TRIPS trips;
  for(sizeType r=0; r<3; r++)
    for(sizeType c=0; c<3; c++)
      for(sizeType k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,sizeType>(r+c*3,k+r*3,m(k,c)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}
template <typename T>
Eigen::SparseMatrix<T,0,sizeType> mat3x3MulLeft(const Eigen::Matrix<T,3,3>& m)
{
  typedef ParallelVector<Eigen::Triplet<T,sizeType> > TRIPS;
  Eigen::SparseMatrix<T,0,sizeType> M;
  M.resize(9,9);
  TRIPS trips;
  for(sizeType r=0; r<3; r++)
    for(sizeType c=0; c<3; c++)
      for(sizeType k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,sizeType>(r+c*3,k+c*3,m(r,k)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}
template <typename T>
Eigen::SparseMatrix<T,0,sizeType> mat3x3MulTLeft(const Eigen::Matrix<T,3,3>& m)
{
  typedef ParallelVector<Eigen::Triplet<T,sizeType> > TRIPS;
  Eigen::SparseMatrix<T,0,sizeType> M;
  M.resize(9,9);
  TRIPS trips;
  for(sizeType r=0; r<3; r++)
    for(sizeType c=0; c<3; c++)
      for(sizeType k=0; k<3; k++)
        trips.push_back(Eigen::Triplet<T,sizeType>(r+c*3,k+c*3,m(k,r)));
  M.setFromTriplets(trips.begin(),trips.end());
  return M;
}

PRJ_END

#endif
