#include "FGTTreeNode.h"
#include <Utils/RotationUtil.h>
#include <Utils/DebugGradient.h>

USE_PRJ_NAMESPACE

template <typename T>
struct AtomicAdd
{
  static void eval(T& A,T B) {
    OMP_CRITICAL_ {
      A+=B;
    }
  }
};
template <>
struct AtomicAdd<scalarD>
{
  static void eval(scalarD& A,scalarD B) {
    OMP_ATOMIC_
    A+=B;
  }
};
#define SY_COND (Sy?Sy->coeff(idy):1)
template <typename T>
FGTTreeNode<T>::FGTTreeNode() {}
template <typename T>
FGTTreeNode<T>::FGTTreeNode(Vec* Sy,Mat3XT& xy,const Vec2i& range,sizeType leafThres,Mat3XT* xyn):_Fs(range[1]-range[0]) {
  _range=range;
  if(range[1]-range[0]<leafThres) {
    if(Sy)
      _Fs=Sy->segment(range[0],range[1]-range[0]).sum();
    for(sizeType i=range[0]; i<range[1]; i++)
      _bbl.setUnion(xy.col(i));
    _bb=_bbl;
    _spherel=_sphere=computeSphere(xy,(_bbl._minC+_bbl._maxC)/2);
  } else {
    //compute variance
    sizeType D;
    Vec3T var=variance(xy);
    var.maxCoeff(&D);
    //sort along D
    std::vector<sizeType> id;
    for(sizeType i=range[0]; i<range[1]; i++)
      id.push_back(i);
    std::sort(id.begin(),id.end(),[&](sizeType a,sizeType b) {
      return xy(D,a)<xy(D,b);
    });
    if(Sy)
      swapId(id,[&](sizeType i,sizeType j) {
      ASSERT(i>=range[0] && i<range[1] && j>=range[0] && j<range[1])
      T tmp=Sy->coeffRef(i);
      Sy->coeffRef(i)=Sy->coeffRef(j);
      Sy->coeffRef(j)=tmp;
    });
    swapId(id,[&](sizeType i,sizeType j) {
      ASSERT(i>=range[0] && i<range[1] && j>=range[0] && j<range[1])
      Vec3T tmp=xy.col(i);
      xy.col(i)=xy.col(j);
      xy.col(j)=tmp;
    });
    if(xyn)
      swapId(id,[&](sizeType i,sizeType j) {
      ASSERT(i>=range[0] && i<range[1] && j>=range[0] && j<range[1])
      Vec3T tmp=xyn->col(i);
      xyn->col(i)=xyn->col(j);
      xyn->col(j)=tmp;
    });
    //divide
    sizeType mid=(range[0]+range[1])/2;
    _l.reset(new FGTTreeNode(Sy,xy,Vec2i(range[0],mid),leafThres,xyn));
    _r.reset(new FGTTreeNode(Sy,xy,Vec2i(mid,range[1]),leafThres,xyn));
    _Fs=_l->_Fs+_r->_Fs;
    //merge
    _bbl=_l->_bbl;
    _bbl.setUnion(_r->_bbl);
    _bb=_bbl;
    _spherel=_sphere=mergeSphere(_l->_spherel,_r->_spherel);
  }
}
template <typename T>
void FGTTreeNode<T>::parityCheck(const Vec* Sy,const Mat3XT& xy,T thres) const
{
  for(sizeType i=_range[0]; i<_range[1]; i++)  {
    ASSERT_MSGV(_bb.distTo(xy.col(i))<thres,"Tree BBox err: dist=%f",_bb.distTo(xy.col(i)))
    ASSERT_MSGV(_sphere.distTo(Sphere<T>(xy.col(i),0))<thres,"Tree Sphere err: dist=%f",_sphere.distTo(Sphere<T>(xy.col(i),0)))
  }
  if(_l) {
    if(Sy) {
      ASSERT_MSGV(std::abs(_l->_Fs+_r->_Fs-_Fs)<thres,"Tree Fs err: totalFs=%f, leftFs=%f, rightFs=%f, totalFs=%f",_Fs,_l->_Fs,_r->_Fs,_l->_Fs+_r->_Fs)
    }
    ASSERT_MSGV(size()==_l->size()+_r->size(),"Tree size err: totalSz=%d, leftSz=%d, rightSz=%d, totalSz=%d",size(),_l->size(),_r->size(),_l->size()+_r->size())
    _l->parityCheck(Sy,xy,thres);
    _r->parityCheck(Sy,xy,thres);
  }
}
template <typename T>
void FGTTreeNode<T>::transform(const Mat3T& R,const Vec3T& t) {
  _sphere._ctr=R*_spherel._ctr+t;
  Vec3T c;
  _bb.reset();
  for(sizeType x=0; x<2; x++) {
    c[0]=x==0?_bbl._minC[0]:_bbl._maxC[0];
    for(sizeType y=0; y<2; y++) {
      c[1]=y==0?_bbl._minC[1]:_bbl._maxC[1];
      for(sizeType z=0; z<2; z++) {
        c[2]=z==0?_bbl._minC[2]:_bbl._maxC[2];
        _bb.setUnion(R*c+t);
      }
    }
  }
  if(_l) {
    _l->transform(R,t);
    _r->transform(R,t);
  }
}
template <typename T>
typename FGTTreeNode<T>::Vec3T FGTTreeNode<T>::variance(Mat3XT& y) const {
  Vec3T mean=Vec3T::Zero();
  for(sizeType i=_range[0]; i<_range[1]; i++)
    mean+=y.col(i);
  mean/=size();

  Vec3T ret=Vec3T::Zero();
  for(sizeType i=_range[0]; i<_range[1]; i++)
    ret.array()+=(y.col(i)-mean).array().square();
  return ret;
}
template <typename T>
Sphere<T> FGTTreeNode<T>::computeSphere(const Mat3XT& y,const Vec3T& ctr) const {
  Sphere<T> s;
  s._ctr=ctr;
  s._rad=0;
  for(sizeType i=_range[0]; i<_range[1]; i++)
    s._rad=std::max<T>(std::sqrt((y.col(i)-ctr).squaredNorm()),s._rad);
  return s;
}
template <typename T>
Sphere<T> FGTTreeNode<T>::mergeSphere(const Sphere<T>& l,const Sphere<T>& r) const {
  Vec3T l2r=l._ctr-r._ctr;
  T dist=std::sqrt(l2r.squaredNorm());
  if(dist+r._rad<=l._rad)
    return l;
  else if(dist+l._rad<r._rad)
    return r;
  else {
    Sphere<T> s;
    Vec3T l2rn=l2r/dist;
    Vec3T left=l._ctr+l2rn*l._rad;
    Vec3T right=r._ctr-l2rn*r._rad;
    s._ctr=(left+right)/2;
    s._rad=std::sqrt((left-right).squaredNorm())/2;
    return s;
  }
}
template <typename T>
T FGTTreeNode<T>::distTo(const FGTTreeNode<T>& other) const
{
  return std::max<T>(_sphere.distTo(other._sphere),_bb.distTo(other._bb));
}
template <typename T>
sizeType FGTTreeNode<T>::size() const
{
  return _range[1]-_range[0];
}
//FGT
template <typename T>
void FGTTreeNode<T>::closestYNode(const FGTTreeNode<T>** minLeaf,T& minDist,const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode)
{
  if(!yNode._l) {
    T minDistCurr=xNode.distTo(yNode);
    if(minDistCurr<minDist) {
      minDist=minDistCurr;
      *minLeaf=&yNode;
    }
  } else {
    T minDistCurr=xNode.distTo(yNode);
    if(minDistCurr>minDist)
      return;
    else {
      closestYNode(minLeaf,minDist,*(yNode._l),xNode);
      closestYNode(minLeaf,minDist,*(yNode._r),xNode);
    }
  }
}
template <typename T>
void FGTTreeNode<T>::initErrorBound(const Vec* Sy,const Mat3XT& y,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr)
{
  if(xNode._l) {
    initErrorBound(Sy,y,x,yNode,*(xNode._l),invHSqr);
    initErrorBound(Sy,y,x,yNode,*(xNode._r),invHSqr);
    xNode._tildeGMinInit=std::min(xNode._l->_tildeGMinInit,xNode._r->_tildeGMinInit);
    xNode._tildeGMin=0;
    xNode._FtSave=0;
  } else {
    const FGTTreeNode<T>* yLeaf=NULL;
    T minDist=ScalarUtil<T>::scalar_max();
    closestYNode(&yLeaf,minDist,yNode,xNode);
    xNode._tildeGMinInit=ScalarUtil<T>::scalar_max();
    for(sizeType idx=xNode._range[0]; idx<xNode._range[1]; idx++) {
      T GMin=0;
      for(sizeType idy=yLeaf->_range[0]; idy<yLeaf->_range[1]; idy++)
        GMin+=std::exp(-(y.col(idy)-x.col(idx)).squaredNorm()*invHSqr)*SY_COND;
      xNode._tildeGMinInit=std::min<T>(xNode._tildeGMinInit,GMin);
    }
    xNode._tildeGMin=0;
    xNode._FtSave=0;
  }
}
template <typename T>
void FGTTreeNode<T>::FGT(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr,T eps,Vec3i* profile)
{
  if(profile)
    profile->setZero();
  initErrorBound(Sy,y,x,yNode,xNode,invHSqr);
  contrib(G,DGDT,Sy,y,yl,x,yNode,xNode,yNode._Fs,invHSqr,eps,profile);
}
template <typename T>
void FGTTreeNode<T>::contrib(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T F,T invHSqr,T eps,Vec3i* profile)
{
  Vec2T minMax;
#define MEAN_CTR
#ifdef MEAN_CTR
  T errMean=errorMeanCtr(yNode,xNode,minMax,invHSqr);
#else
  T errMean=errorMean(yNode,xNode,minMax,invHSqr);
#endif
  T errBound=eps*std::max(xNode._tildeGMin,xNode._tildeGMinInit)*(yNode._Fs+xNode._FtSave)/F;
  if(errMean<errBound) {
#ifdef MEAN_CTR
    meanCtr(G,DGDT,yNode,xNode,F,invHSqr,errMean,eps);
#else
    mean(G,DGDT,yNode,xNode,minMax,F,invHSqr,errMean,eps);
#endif
    if(profile)
      profile->coeffRef(0)++;
  } else {
    if(!yNode._l || !xNode._l) {
      direct(G,DGDT,Sy,y,yl,x,yNode,xNode,invHSqr);
      if(profile)
        profile->coeffRef(1)++;
    } else {
      sizeType p;
      sizeType cDirect=costDirect(yNode,xNode);
      sizeType cTaylor=costTaylor(yNode,xNode,p,invHSqr,errBound);
      if(cDirect<=cTaylor || costChildren(yNode,xNode,invHSqr,errBound)<=cTaylor) {
        contrib(G,DGDT,Sy,y,yl,x,*(yNode._l),*(xNode._l),F,invHSqr,eps,profile);
        contrib(G,DGDT,Sy,y,yl,x,*(yNode._r),*(xNode._l),F,invHSqr,eps,profile);
        contrib(G,DGDT,Sy,y,yl,x,*(yNode._l),*(xNode._r),F,invHSqr,eps,profile);
        contrib(G,DGDT,Sy,y,yl,x,*(yNode._r),*(xNode._r),F,invHSqr,eps,profile);
      } else {
        taylor(G,DGDT,Sy,y,yl,x,yNode,xNode,invHSqr,p);
        if(profile)
          profile->coeffRef(2)++;
      }
    }
  }
}
template <typename T>
sizeType FGTTreeNode<T>::cost(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,sizeType& p,T invHSqr,T errBound)
{
  sizeType cDirect=costDirect(yNode,xNode);
  sizeType cTaylor=costTaylor(yNode,xNode,p,invHSqr,errBound);
  return std::min<sizeType>(cDirect,cTaylor);
}
template <typename T>
sizeType FGTTreeNode<T>::costChildren(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,T invHSqr,T errBound)
{
  sizeType p,ret=0;
  ret+=cost(*(yNode._l),*(xNode._l),p,invHSqr,errBound);
  ret+=cost(*(yNode._r),*(xNode._l),p,invHSqr,errBound);
  ret+=cost(*(yNode._l),*(xNode._r),p,invHSqr,errBound);
  ret+=cost(*(yNode._r),*(xNode._r),p,invHSqr,errBound);
  return ret;
}
//mean
template <typename T>
void FGTTreeNode<T>::mean(Vec& G,MatX4T* DGDT,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,const Vec2T& minMax,T F,T invHSqr,T errMean,T eps)
{
  ASSERT_MSG(!DGDT,"Mean does not support DGDT")
  T val=yNode._Fs*0.5f*(std::exp(-minMax[0]*minMax[0]*invHSqr)+std::exp(-minMax[1]*minMax[1]*invHSqr));
  G.segment(xNode._range[0],xNode.size()).array()+=val;
  xNode._tildeGMin=G.segment(xNode._range[0],xNode.size()).minCoeff();
  xNode._FtSave+=yNode._Fs-errMean*F/(eps*std::max(xNode._tildeGMin,xNode._tildeGMinInit));
}
template <typename T>
T FGTTreeNode<T>::errorMean(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,Vec2T& minMax,T invHSqr)
{
  minMax[0]=yNode.distTo(xNode);
  Vec3T distBox;
  for(sizeType d=0; d<3; d++) {
    distBox[d]=std::max<T>(yNode._bb._maxC[d],xNode._bb._maxC[d]);
    distBox[d]-=std::min<T>(yNode._bb._minC[d],xNode._bb._minC[d]);
  }
  T distSphere=std::sqrt((yNode._sphere._ctr-xNode._sphere._ctr).squaredNorm())+yNode._sphere._rad+xNode._sphere._rad;
  minMax[1]=std::min(distSphere,std::sqrt(distBox.squaredNorm()));
  return yNode._Fs*0.5f*(std::exp(-minMax[0]*minMax[0]*invHSqr)-std::exp(-minMax[1]*minMax[1]*invHSqr));
}
template <typename T>
sizeType FGTTreeNode<T>::costMean(const FGTTreeNode<T>&,const FGTTreeNode<T>& xNode)
{
  return 3+xNode.size();
}
//mean-ctr
template <typename T>
void FGTTreeNode<T>::meanCtr(Vec& G,MatX4T* DGDT,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T F,T invHSqr,T errMean,T eps)
{
  Vec3T dir=yNode._sphere._ctr-xNode._sphere._ctr;
  T coef=yNode._Fs*std::exp(-dir.squaredNorm()*invHSqr);
  G.segment(xNode._range[0],xNode.size()).array()+=coef;
  if(DGDT) {
    Mat3X4T DValDT=dir*Vec4T(yNode._spherel._ctr[0],yNode._spherel._ctr[1],yNode._spherel._ctr[2],1).transpose()*(coef*2*invHSqr);
    OMP_PARALLEL_FOR_
    for(sizeType i=xNode._range[0]; i<xNode._range[1]; i++)
      DGDT->template block<3,4>(i*3,0)-=DValDT;
  }
  xNode._tildeGMin=G.segment(xNode._range[0],xNode.size()).minCoeff();
  xNode._FtSave+=yNode._Fs-errMean*F/(eps*std::max(xNode._tildeGMin,xNode._tildeGMinInit));
}
template <typename T>
T FGTTreeNode<T>::errorMeanCtr(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,Vec2T& minMax,T invHSqr)
{
  return errorMean(yNode,xNode,minMax,invHSqr)*2;
}
template <typename T>
sizeType FGTTreeNode<T>::costMeanCtr(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode)
{
  return costMean(yNode,xNode);
}
//direct
template <typename T>
void FGTTreeNode<T>::direct(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr)
{
  T coef;
  Vec3T dir;
  OMP_PARALLEL_FOR_I(OMP_PRI(coef,dir))
  for(sizeType idx=xNode._range[0]; idx<xNode._range[1]; idx++) {
    for(sizeType idy=yNode._range[0]; idy<yNode._range[1]; idy++) {
      dir=y.col(idy)-x.col(idx);
      coef=SY_COND*std::exp(-dir.squaredNorm()*invHSqr);
      G[idx]+=coef;
      if(DGDT) {
        Eigen::Block<MatX4T,3,4> blk=DGDT->template block<3,4>(idx*3,0);
        blk-=dir*Vec4T(yl->coeff(0,idy),yl->coeff(1,idy),yl->coeff(2,idy),1).transpose()*(coef*2*invHSqr);
      }
    }
  }
  xNode._tildeGMin=G.segment(xNode._range[0],xNode.size()).minCoeff();
  xNode._FtSave+=yNode._Fs;
}
template <typename T>
sizeType FGTTreeNode<T>::costDirect(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode)
{
  return 3*yNode.size()*xNode.size();
}
//taylor
template <typename T>
void FGTTreeNode<T>::taylor(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr,T eps)
{
  sizeType p;
  costTaylor(yNode,xNode,p,invHSqr,eps);
  taylor(G,DGDT,Sy,y,yl,x,yNode,xNode,invHSqr,p);
}
template <typename T>
void FGTTreeNode<T>::taylor(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr,sizeType p)
{
  sizeType pSqr=p*p;
  T invH=std::sqrt(invHSqr);

  //build M
  Vec M[4];
  MatX4T DMDT;
  M[3].setZero(p*p*p);
  if(DGDT) {
    M[0]=M[1]=M[2]=M[3];
    DMDT.setZero(3*p*p*p,4);
  }
  OMP_PARALLEL_FOR_
  for(sizeType idy=yNode._range[0]; idy<yNode._range[1]; idy++) {
    //build coef
    Vec4T ylH=Vec4T((*yl)(0,idy),(*yl)(1,idy),(*yl)(2,idy),1);
    T coef=SY_COND*std::exp(-(y.col(idy)-yNode._sphere._ctr).squaredNorm()*invHSqr);
    coef*=std::exp(2*(y.col(idy)-yNode._sphere._ctr).dot(xNode._sphere._ctr-yNode._sphere._ctr)*invHSqr);
    //build cache
    Mat3XT cache;
    cache.setOnes(3,p);
    Vec3T dy=(y.col(idy)-yNode._sphere._ctr)*invH;
    for(sizeType i=1; i<p; i++)
      cache.col(i).array()=cache.col(i-1).array()*dy.array()*2/i;
    dy=(y.col(idy)-xNode._sphere._ctr)*(-2*invHSqr);
    //accumulate M
    for(sizeType a0=0,off0=0; a0<p; a0++,off0+=pSqr)
      for(sizeType a1=0,off1=off0; a0+a1<p; a1++,off1+=p)
        for(sizeType a2=0,off2=off1; a0+a1+a2<p; a2++,off2++) {
          T coefCache=cache(0,a0)*cache(1,a1)*cache(2,a2)*coef;
          if(DGDT) {
            for(sizeType r=0,offr=off2*3; r<3; r++,offr++) {
              AtomicAdd<T>::eval(M[r][off2],ylH[r]*coefCache);
              for(sizeType c=0; c<4; c++)
                AtomicAdd<T>::eval(DMDT(offr,c),dy[r]*ylH[c]*coefCache);
            }
          }
          AtomicAdd<T>::eval(M[3][off2],coefCache);
        }
  }

  //use M
  OMP_PARALLEL_FOR_
  for(sizeType idx=xNode._range[0]; idx<xNode._range[1]; idx++) {
    //build coef
    T coef=std::exp(-(x.col(idx)-yNode._sphere._ctr).squaredNorm()*invHSqr);
    //build cache
    Mat3XT cache;
    cache.setOnes(3,p);
    Vec3T dx=(x.col(idx)-xNode._sphere._ctr)*invH;
    for(sizeType i=1; i<p; i++)
      cache.col(i).array()=cache.col(i-1).array()*dx.array();
    dx=(x.col(idx)-xNode._sphere._ctr)*(2*invHSqr);
    //assumulate G
    Eigen::Block<MatX4T,3,4> blk=DGDT->template block<3,4>(idx*3,0);
    for(sizeType a0=0,off0=0; a0<p; a0++,off0+=pSqr)
      for(sizeType a1=0,off1=off0; a0+a1<p; a1++,off1+=p)
        for(sizeType a2=0,off2=off1; a0+a1+a2<p; a2++,off2++) {
          T coefCache=cache(0,a0)*cache(1,a1)*cache(2,a2)*coef;
          if(DGDT) {
            blk.col(0)+=dx*M[0][off2]*coefCache;
            blk.col(1)+=dx*M[1][off2]*coefCache;
            blk.col(2)+=dx*M[2][off2]*coefCache;
            blk.col(3)+=dx*M[3][off2]*coefCache;
            blk+=DMDT.template block<3,4>(off2*3,0)*coefCache;
          }
          G[idx]+=M[3][off2]*coefCache;
        }
  }
  xNode._tildeGMin=G.segment(xNode._range[0],xNode.size()).minCoeff();
}
template <typename T>
sizeType FGTTreeNode<T>::costTaylor(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,sizeType& p,T invHSqr,T errBound)
{
  T radSum=yNode._sphere._rad+xNode._sphere._rad;
  T term1=2*(yNode._sphere._rad*xNode._sphere._rad)*invHSqr;
  T term2=std::exp(radSum*radSum*invHSqr);
  T errTaylor=term1*term2*yNode._Fs;
  for(p=1;; p++) {
    if(errTaylor<errBound || p>100)
      break;
    errTaylor*=term1/(p+1);
  }
  //T errTaylorRef=yNode._Fs/factorial(p)*std::pow(term1,p)*term2;
  return combination(p-1+3,3)*(yNode.size()+xNode.size());
}
//heper
template <typename T>
void FGTTreeNode<T>::debugSwapId(sizeType base,sizeType N)
{
  std::vector<sizeType> id;
  for(sizeType i=0; i<N; i++)
    id.push_back(base+i);
  Vec data=Vec::Random(base+N),data2=data;
  //swayIdRef
  for(sizeType i=0; i<N; i++)
    data2[i+base]=data[id[i]];
  //swapId
  swapId(id,[&](sizeType i,sizeType j) {
    T tmp=data[i];
    data[i]=data[j];
    data[j]=tmp;
  });
  DEFINE_NUMERIC_DELTA_T(T)
  DEBUG_GRADIENT("DebugSwapId",std::sqrt(data.squaredNorm()),std::sqrt((data-data2).squaredNorm()))
}
template <typename T>
void FGTTreeNode<T>::swapId(std::vector<sizeType> id,std::function<void(sizeType,sizeType)> f)
{
  sizeType base=*(std::minmax_element(id.begin(),id.end()).first);
  for(sizeType i=0; i<(sizeType)id.size(); i++) {
    sizeType curr=i;
    while(id[curr]!=curr+base) {
      std::swap(id[curr],id[id[curr]-base]);
      f(id[curr],curr+base);
      curr=id[curr]-base;
    }
  }
}
template <typename T>
void FGTTreeNode<T>::randomTree(Vec* Sy,Mat3XT& y,Mat3XT& yl,Mat3XT& x,FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,sizeType N,sizeType leafThres,bool random)
{
  if(random) {
    yl=(Mat3XT::Random(3,N*N*N)+Mat3XT::Ones(3,N*N*N))/2;
    x=(Mat3XT::Random(3,N*N*N)+Mat3XT::Ones(3,N*N*N))/2;
  } else {
    yl=Mat3XT::Zero(3,N*N*N);
    for(sizeType X=0,off=0; X<N; X++)
      for(sizeType Y=0; Y<N; Y++)
        for(sizeType Z=0; Z<N; Z++,off++)
          yl.col(off)=Vec3T((X+0.5f)/N,(Y+0.5f)/N,(Z+0.5f)/N);

    x=Mat3XT::Zero(3,N*N*N);
    for(sizeType X=0,off=0; X<N; X++)
      for(sizeType Y=0; Y<N; Y++)
        for(sizeType Z=0; Z<N; Z++,off++)
          x.col(off)=Vec3T((X+0.5f)/N,(Y+0.5f)/N,(Z+0.5f)/N);
  }
  if(Sy)
    *Sy=(Vec::Random(yl.cols())+Vec::Ones(yl.cols()))/2;
  yNode=FGTTreeNode<T>(Sy,yl,Vec2i(0,yl.cols()),leafThres);
  xNode=FGTTreeNode<T>(NULL,x,Vec2i(0,yl.cols()),leafThres);
  yNode.parityCheck(Sy,yl);
  xNode.parityCheck(NULL,x);
  y=yl;
}
template <typename T>
void FGTTreeNode<T>::debugTaylor(sizeType N,sizeType leafThres,bool random,T invHSqr,T eps,bool useSy) {
  Vec Sy;
  Mat3XT y,x,yl;
  FGTTreeNode<T> yNode,xNode;
  randomTree(useSy?&Sy:NULL,y,yl,x,yNode,xNode,N,leafThres,random);
  std::string strSy=useSy?"-Sy":"";

  DEFINE_NUMERIC_DELTA_T(T)
  Vec GDirect=Vec::Zero(y.cols());
  Vec GTaylor=Vec::Zero(y.cols());
  MatX4T DGDTDirect=MatX4T::Zero(3*y.cols(),4);
  MatX4T DGDTTaylor=MatX4T::Zero(3*y.cols(),4);
  direct(GDirect,&DGDTDirect,useSy?&Sy:NULL,y,&yl,x,yNode,xNode,invHSqr);
  taylor(GTaylor,&DGDTTaylor,useSy?&Sy:NULL,y,&yl,x,yNode,xNode,invHSqr,eps);
  DEBUG_GRADIENT("Taylor-G"+strSy,std::sqrt(GDirect.squaredNorm()),std::sqrt((GDirect-GTaylor).squaredNorm()))
  DEBUG_GRADIENT("Taylor-DGDT"+strSy,std::sqrt(DGDTDirect.squaredNorm()),std::sqrt((DGDTDirect-DGDTTaylor).squaredNorm()))
}
template <typename T>
void FGTTreeNode<T>::debugFGT(sizeType N,sizeType leafThres,bool random,T invHSqr,T eps,bool useSy) {
  Vec Sy;
  Mat3XT y,x,yl;
  FGTTreeNode<T> yNode,xNode;
  randomTree(useSy?&Sy:NULL,y,yl,x,yNode,xNode,N,leafThres,random);
  std::string strSy=useSy?"-Sy":"";

  Vec3i profile;
  Vec GDirect=Vec::Zero(y.cols());
  Vec GTaylor=Vec::Zero(y.cols());
  MatX4T DGDTDirect=MatX4T::Zero(3*y.cols(),4);
  MatX4T DGDTTaylor=MatX4T::Zero(3*y.cols(),4);
  direct(GDirect,&DGDTDirect,useSy?&Sy:NULL,y,&yl,x,yNode,xNode,invHSqr);
  FGT(GTaylor,&DGDTTaylor,useSy?&Sy:NULL,y,&yl,x,yNode,xNode,invHSqr,eps,&profile);

  DEFINE_NUMERIC_DELTA_T(T)
  DEBUG_GRADIENT("FGT-G"+strSy,std::sqrt(GDirect.squaredNorm()),std::sqrt((GDirect-GTaylor).squaredNorm()))
  DEBUG_GRADIENT("FGT-DGDT"+strSy,std::sqrt(DGDTDirect.squaredNorm()),std::sqrt((DGDTDirect-DGDTTaylor).squaredNorm()))
  INFOV("#mean=%d #direct=%d #taylor=%d",profile[0],profile[1],profile[2])
}
template <typename T>
void FGTTreeNode<T>::debugTree(sizeType N,sizeType leafThres,bool random,bool useSy) {
  Vec Sy;
  Mat3XT y,x,yl;
  FGTTreeNode<T> yNode,xNode;
  randomTree(useSy?&Sy:NULL,y,yl,x,yNode,xNode,N,leafThres,random);

  Mat3T R=expWGradV<T,Vec3T>(Vec3T::Random()*M_PI);
  Vec3T t=Vec3T::Random();
  y=R*yl+t*Vec::Ones(yl.cols()).transpose();
  yNode.transform(R,t);
  yNode.parityCheck(&Sy,y);
}
template <typename T>
sizeType FGTTreeNode<T>::combination(sizeType d,sizeType k)
{
  ASSERT(d>=k)
  //return factorial(d)/factorial(k)/factorial(d-k);
  sizeType ret=1;
  for(sizeType i=d-k+1; i<=d; i++)
    ret*=i;
  return ret/factorial(k);
}
template <typename T>
sizeType FGTTreeNode<T>::factorial(sizeType p)
{
  if(p<=1)
    return 1;
  sizeType ret=1;
  for(sizeType i=2; i<=p; i++)
    ret*=i;
  return ret;
}
//instance
PRJ_BEGIN
template class FGTTreeNode<double>;
#ifdef ALL_TYPES
template class FGTTreeNode<__float128>;
template class FGTTreeNode<mpfr::mpreal>;
#endif
PRJ_END
