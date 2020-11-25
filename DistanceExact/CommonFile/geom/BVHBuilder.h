#ifndef BVH_BUILDER_H
#define BVH_BUILDER_H

#include "../MathBasic.h"
#include "../IO.h"
#include "BVHNode.h"
#include <deque>
#include <stack>
#include <type_traits>

PRJ_BEGIN

template <int D>struct SurfaceArea;
template <>struct SurfaceArea<3> {
  template <typename BoxType>
  static scalar area(const BoxType& bb) {
    typedef typename BoxType::PT PT;
    PT ext=compMax(PT::Zero(),(PT)(bb.maxCorner()-bb.minCorner()));
    return (ext[0]*ext[1]+ext[0]*ext[2]+ext[1]*ext[2])*2.0f;
  }
};
template <>struct SurfaceArea<2> {
  template <typename BoxType>
  static scalar area(const BoxType& bb) {
    typedef typename BoxType::PT PT;
    PT ext=compMax(PT::Zero(),(PT)(bb.maxCorner()-bb.minCorner()));
    return (ext[0]+ext[1])*2.0f;
  }
};
template <typename NODE,int DIM>
struct BVHBuilder {
public:
  static_assert(DIM == 2 || DIM == 3,"Incorrect BVH Dimension!");
  static const int dim=NODE::dim;
  typedef typename NODE::BoxType BoxType;
  typedef typename BoxType::PT PT;
  struct BVHHandle {
    BVHHandle():_left(false),_cost(0.0f) {}
    BVHHandle(scalar val,sizeType rid,bool left):_val(val),_cost(0.0f),_rid(rid),_left(left) {}
    bool operator<(const BVHHandle& other) const {
      return _val<other._val;
    }
    bool operator<=(const BVHHandle& other) const {
      return _val<=other._val;
    }
    bool operator>(const BVHHandle& other) const {
      return _val>other._val;
    }
    bool operator>=(const BVHHandle& other) const {
      return _val>=other._val;
    }
    bool operator==(const BVHHandle& other) const {
      return _val==other._val;
    }
    scalar _val,_cost;
    sizeType _rid;
    bool _left;
  };
  sizeType buildBVH(std::vector<NODE>& bvh,sizeType f,sizeType t) {
    if(bvh.empty())
      return -1;
    //find roots
    _roots.clear();
    for(sizeType i=f; i<t; i++) {
      ASSERT(bvh[i]._parent == -1)
      _roots.push_back(i);
    }
    return buildBVHInner(bvh);
  }
  sizeType buildBVH(std::vector<NODE>& bvh) {
    if(bvh.empty())
      return -1;
    //find roots
    _roots.clear();
    for(sizeType i=0; i<(sizeType)bvh.size(); i++)
      if(bvh[i]._parent == -1)
        _roots.push_back(i);
    return buildBVHInner(bvh);
  }
private:
  sizeType buildBVHInner(std::vector<NODE>& bvh) {
    //short bounding box alone each axis
    for(sizeType d=0; d<DIM; d++) {
      _hdls[d].clear();
      for(sizeType i=0; i<(sizeType)_roots.size(); i++) {
        _hdls[d].push_back(BVHHandle(bvh[_roots[i]]._bb.minCorner()[d],_roots[i],true));
        _hdls[d].push_back(BVHHandle(bvh[_roots[i]]._bb.maxCorner()[d],_roots[i],false));
      }
      std::sort(_hdls[d].begin(),_hdls[d].end());
    }
    //build BVH
    Vec3i f(0,0,0);
    Vec3i t((sizeType)_hdls[0].size(),
            (sizeType)_hdls[1].size(),
            (sizeType)_hdls[2].size());
    bvh.reserve(bvh.size()+_roots.size()-1);
    sizeType ret=buildBVH(bvh,0,(sizeType)_roots.size(),f,t);
    bvh.back()._parent=-1;
    return ret;
  }
//#define DEBUG_BVH
  sizeType buildBVH(std::vector<NODE>& bvh,sizeType fr,sizeType tr,Vec3i f,Vec3i t) {
    if(tr == fr+1)
      return _roots[fr];
    //find split location
    sizeType lt=fr,rf=tr;
    Vec3i ltd=f,rfd=t;
    {
      //calculate split cost
      sizeType minD=-1,minId=-1;
      scalar minCost=ScalarUtil<scalar>::scalar_max();
      for(sizeType d=0; d<DIM; d++) {
#ifdef DEBUG_BVH
        {
          std::vector<BVHHandle> hdls;
          for(sizeType i=fr; i<tr; i++) {
            hdls.push_back(BVHHandle(bvh[_roots[i]]._bb._minC[d],_roots[i],true));
            hdls.push_back(BVHHandle(bvh[_roots[i]]._bb._maxC[d],_roots[i],false));
          }
          std::sort(hdls.begin(),hdls.end());

          ASSERT((sizeType)hdls.size() == t[d]-f[d])
          ASSERT((sizeType)hdls.size() == (tr-fr)*2)
          for(sizeType i=0; i<(sizeType)hdls.size(); i++)
            ASSERT(hdls[i] == _hdls[d][i+f[d]]);
        }
#endif
        calcCost(bvh,_hdls[d],f[d],t[d]);
        for(sizeType i=f[d]; i<(sizeType)t[d]; i++) {
#ifdef DEBUG_BVH
          ASSERT(debugCost(bvh,_hdls[d],f[d],t[d],i) == _hdls[d][i]._cost);
#endif
          if(_hdls[d][i]._cost < minCost) {
            minCost=_hdls[d][i]._cost;
            minD=d;
            minId=i;
          }
        }
      }
      //find split set
      for(sizeType i=fr; i<tr; i++)
        bvh[_roots[i]]._parent=-1;
      const BVHHandle& minHdl=_hdls[minD][minId];
      for(sizeType i=f[minD]; i<t[minD]; i++) {
        const BVHHandle& hdl=_hdls[minD][i];
        if(bvh[hdl._rid]._parent != -1)
          continue;
        if(hdl._left && hdl <= minHdl) {
          _roots[lt++]=hdl._rid;
          bvh[hdl._rid]._parent=1;
        }
        if(!hdl._left && hdl >= minHdl) {
          _roots[--rf]=hdl._rid;
          bvh[hdl._rid]._parent=2;
        }
      }
      ASSERT_MSG(lt == rf,"Strange Error!");
      if(lt == tr || rf == fr) {
        lt=rf=(fr+tr)/2;
        for(sizeType i=fr; i<lt; i++)
          bvh[_roots[i]]._parent=1;
        for(sizeType i=rf; i<tr; i++)
          bvh[_roots[i]]._parent=2;
      }
      for(sizeType d=0; d<DIM; d++)
        split(bvh,_hdls[d],ltd[d],rfd[d]);
    }

    //split
    NODE n;
    n._l=buildBVH(bvh,fr,lt,f,ltd);
    n._r=buildBVH(bvh,rf,tr,rfd,t);
    bvh[n._l]._parent=(sizeType)bvh.size();
    bvh[n._r]._parent=(sizeType)bvh.size();

    //add new root
    n._bb=bvh[n._l]._bb;
    n._bb.setUnion(bvh[n._r]._bb);
    n._nrCell=bvh[n._l]._nrCell+bvh[n._r]._nrCell;
    bvh.push_back(n);
    return bvh.size()-1;
  }
  static void calcCost(const std::vector<NODE>& bvh,std::vector<BVHHandle>& X,sizeType f,sizeType t) {
    BoxType bb;
    sizeType nrC=0;
    for(sizeType i=f; i<t; i++) {
      X[i]._cost=SurfaceArea<DIM>::area(bb)*scalar((long)nrC);
      if(X[i]._left) {
        PT minC=bvh[X[i]._rid]._bb.minCorner();
        PT maxC=bvh[X[i]._rid]._bb.maxCorner();
        bb.setUnion(BoxType(minC,maxC));
        nrC+=bvh[X[i]._rid]._nrCell;
      }
    }
    bb.reset();
    nrC=0;
    for(sizeType i=t-1; i>=f; i--) {
      X[i]._cost+=SurfaceArea<DIM>::area(bb)*scalar((long)nrC);
      if(!X[i]._left) {
        PT minC=bvh[X[i]._rid]._bb.minCorner();
        PT maxC=bvh[X[i]._rid]._bb.maxCorner();
        bb.setUnion(BoxType(minC,maxC));
        nrC+=bvh[X[i]._rid]._nrCell;
      }
    }
  }
  static scalar debugCost(const std::vector<NODE>& bvh,std::vector<BVHHandle>& X,sizeType f,sizeType t,sizeType I) {
    scalar cost=0.0f;
    sizeType nrC=0;
    BoxType bb;
    for(sizeType i=f; i<I; i++)
      if(X[i]._left) {
        bb.setUnion(bvh[X[i]._rid]._bb);
        nrC+=bvh[X[i]._rid]._nrCell;
      }
    cost+=SurfaceArea<DIM>::area(bb)*(scalar)nrC;

    bb.reset();
    nrC=0;
    for(sizeType i=t-1; i>I; i--)
      if(!X[i]._left) {
        bb.setUnion(bvh[X[i]._rid]._bb);
        nrC+=bvh[X[i]._rid]._nrCell;
      }
    cost+=SurfaceArea<DIM>::area(bb)*(scalar)nrC;

    return cost;
  }
  static void split(const std::vector<NODE>& bvh,std::vector<BVHHandle>& X,sizeType& lt,sizeType& rf) {
    std::vector<BVHHandle> R;
    sizeType j=lt;
    for(sizeType i=lt; i<rf; i++) {
      sizeType p=bvh[X[i]._rid]._parent;
      ASSERT(p == 1 || p == 2)
      if(p == 1) X[j++]=X[i];
      else R.push_back(X[i]);
    }
    lt=rf=j;
    std::copy(R.begin(),R.end(),X.begin()+rf);
  }
  std::vector<BVHHandle> _hdls[3];
  std::vector<sizeType> _roots;
};

template <typename T,typename BBOX>
void buildBVH(std::vector<Node<T,BBOX> >& bvh,sizeType dim,T verbose)
{
  if(bvh.empty())
    return;
  ASSERT(dim == 2 || dim == 3)
  sizeType nrBVH=(sizeType)bvh.size();
  if(dim == 2)BVHBuilder<Node<T,BBOX>,2>().buildBVH(bvh);
  else BVHBuilder<Node<T,BBOX>,3>().buildBVH(bvh);
  for(sizeType i=nrBVH; i<(sizeType)bvh.size(); i++)
    bvh[i]._cell=verbose;
}
template <typename T,typename BBOX>
void buildBVH(std::vector<Node<T,BBOX> >& bvh,sizeType f,sizeType t,sizeType dim,T verbose)
{
  if(bvh.empty())
    return;
  ASSERT(dim == 2 || dim == 3)
  sizeType nrBVH=(sizeType)bvh.size();
  if(dim == 2)BVHBuilder<Node<T,BBOX>,2>().buildBVH(bvh,f,t);
  else BVHBuilder<Node<T,BBOX>,3>().buildBVH(bvh,f,t);
  for(sizeType i=nrBVH; i<(sizeType)bvh.size(); i++)
    bvh[i]._cell=verbose;
}

template <typename T,typename BBOX,typename TV=std::vector<Node<T,BBOX> > >
static void writeBVHByLevel(const TV& bvh,T verbose,const std::string& path=".")
{
  if(bvh.empty())
    return;
  std::deque<sizeType> lv;
  lv.push_back((sizeType)bvh.size()-1);
  for(sizeType i=0; !lv.empty(); i++) {
    std::ostringstream oss;
    oss << (path+"/lv") << i << ".vtk";
    VTKWriter<scalar> os("BB Level",oss.str(),true);

    sizeType nrN=(sizeType)lv.size();
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > bbs;
    for(sizeType i=0; i<nrN; i++) {
      const Node<T,BBOX>& n=bvh[lv.front()];
      BBOX bb=n._bb;
      bbs.push_back(bb.minCorner());
      bbs.push_back(bb.maxCorner());
      if(n._cell == verbose) {
        lv.push_back(n._l);
        lv.push_back(n._r);
      }
      lv.pop_front();
    }
    os.appendVoxels(bbs.begin(),bbs.end(),false);
  }
}
template <typename T,typename BBOX,typename TV=std::vector<Node<T,BBOX> > >
static void writeBVHKDTree(const TV& bvh,sizeType dim,const std::string& path=".")
{
  if(bvh.empty())
    return;
  std::stack<std::pair<sizeType,BBOX>> ss;
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>> iss;
  ss.push(std::make_pair(bvh.size()-1,bvh.back()._bb));

  while(!ss.empty()) {
    sizeType id=ss.top().first;
    if(bvh[id]._l==-1 && bvh[id]._r==-1) {
      ss.pop();
      continue;
    }

    const Node<T,BBOX>& l=bvh[bvh[id]._l];
    const Node<T,BBOX>& r=bvh[bvh[id]._r];
    BBOX bb=ss.top().second,bbl=bb,bbr=bb;
    typename BBOX::PT::Scalar ctr;
    for(sizeType d=0; d<dim; d++) {
      bool valid=false;
      if(l._bb._maxC[d]<=r._bb._minC[d]) {
        ctr=(l._bb._maxC[d]+r._bb._minC[d])/2;
        bbr._minC[d]=bbl._maxC[d]=ctr;
        valid=true;
      } else if(r._bb._maxC[d]<=l._bb._minC[d]) {
        ctr=(r._bb._maxC[d]+l._bb._minC[d])/2;
        bbl._minC[d]=bbr._maxC[d]=ctr;
        valid=true;
      }

      if(valid) {
        sizeType d1=(d+1)%3;
        sizeType d2=(d+2)%3;

        Vec3 x00;
        x00[d]=ctr;
        x00[d1]=bb._minC[d1];
        x00[d2]=bb._minC[d2];

        Vec3 x10;
        x10[d]=ctr;
        x10[d1]=bb._maxC[d1];
        x10[d2]=bb._minC[d2];

        Vec3 x11;
        x11[d]=ctr;
        x11[d1]=bb._maxC[d1];
        x11[d2]=bb._maxC[d2];

        Vec3 x01;
        x01[d]=ctr;
        x01[d1]=bb._minC[d1];
        x01[d2]=bb._maxC[d2];

        sizeType off=(sizeType)vss.size();
        vss.push_back(x00);
        vss.push_back(x10);
        vss.push_back(x11);
        vss.push_back(x01);
        iss.push_back(Vec3i(off+0,off+1,off+2));
        iss.push_back(Vec3i(off+0,off+2,off+3));
        break;
      }
    }
    ss.pop();
    ss.push(std::make_pair(bvh[id]._l,bbl));
    ss.push(std::make_pair(bvh[id]._r,bbr));
  }

  VTKWriter<scalar> os("BB Level",path,true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(iss.begin(),iss.end(),VTKWriter<scalar>::TRIANGLE);
}

template <typename T,typename BBOX=BBox<scalar>,typename TV=std::vector<Node<T,BBOX> > >
class BVHQuery
{
public:
  template <typename T2,typename CACHE>
  struct AddCache {
    AddCache(std::vector<CACHE>& cache):_cache(cache) {}
    void onCell(const Node<T,BBOX>& n,const Node<T2,BBOX>& n2) {
      _cache.push_back(CACHE(n._cell,n2._cell));
    }
    std::vector<CACHE>& _cache;
  };
  BVHQuery(const TV& bvh,sizeType dim,T verbose)
    :_active(NULL),_bvh(bvh),_dim(dim),_verbose(verbose) {}
  //update
  void compact() {
    TV tmp(_bvh.size()-nrEmpty());
    if(tmp.empty()) {
      const_cast<TV&>(_bvh).clear();
      return;
    }
    //assign
    std::stack<std::pair<sizeType,sizeType> > ss;
    tmp[(sizeType)tmp.size()-1]=_bvh.back();
    ss.push(std::make_pair((sizeType)_bvh.size()-1,(sizeType)tmp.size()-1));
    sizeType curr=ss.top().second;
    while(!ss.empty())
      if(_bvh[ss.top().first]._cell == _verbose) {
        sizeType l=_bvh[ss.top().first]._l;
        sizeType r=_bvh[ss.top().first]._r;
        tmp[curr-1]=_bvh[l];
        tmp[curr-2]=_bvh[r];
        tmp[ss.top().second]._l=curr-1;
        tmp[ss.top().second]._r=curr-2;
        tmp[curr-1]._parent=ss.top().second;
        tmp[curr-2]._parent=ss.top().second;
        ss.pop();
        ss.push(std::make_pair(l,curr-1));
        ss.push(std::make_pair(r,curr-2));
        curr-=2;
      } else ss.pop();
    const_cast<TV&>(_bvh)=tmp;
  }
  bool updateBVH(scalar expand=1E-3f,bool dynamic=false) {
    if(hasEmpty(1))
      compact();
    sizeType nrN=(sizeType)_bvh.size();
    bool updated=false;
    for(sizeType i=0; i<nrN; i++) {
      Node<T,BBOX>& n=const_cast<Node<T,BBOX>&>(_bvh[i]);
      if(n._cell != _verbose)
        n._bb.enlarged(expand,_dim);
      else if(n._l >= 0 && n._r >= 0) {
        n._bb=_bvh[n._l]._bb;
        n._bb.setUnion(_bvh[n._r]._bb);
        //dynamic SAH-based adjustment
        if(!dynamic)
          continue;
        const Node<T,BBOX>& l=_bvh[n._l];
        const Node<T,BBOX>& r=_bvh[n._r];
        if(l._cell != _verbose && r._cell == _verbose)
          updated=tryRotate(i,n._r)||updated;
        else if(l._cell == _verbose && r._cell != _verbose)
          updated=tryRotate(i,n._l)||updated;
        else if(l._cell == _verbose && r._cell == _verbose)
          updated=tryRecombine(i)||updated;
      }
    }
    return updated;
  }
  void updateBVH(sizeType i) {
    TV& bvh=const_cast<TV&>(_bvh);
    while(i>=0) {
      bvh[i]._bb=bvh[bvh[i]._l]._bb.getUnion(bvh[bvh[i]._r]._bb);
      i=bvh[i]._parent;
    }
  }
  void removeLeaf(sizeType i) {
    //make sure that first node is empty
    if(!hasEmpty(1)) {
      reserveEmpty(1);
      i+=1;
    }
    TV& bvh=const_cast<TV&>(_bvh);
    Node<T,BBOX>& c=bvh[i];
    ASSERT(c._cell != _verbose)
    //find parent
    if(c._parent == -1)
      bvh.clear();
    else {
      //use parent as sibling
      Node<T,BBOX>& p=bvh[c._parent];
      sizeType sid=p._l == i ? p._r : p._l;
      Node<T,BBOX>& s=bvh[sid];
      p._bb=s._bb;
      p._cell=s._cell;
      p._nrCell=s._nrCell;
      recomputeNrCell(c._parent);
      p._l=s._l;
      if(p._l >= 0)
        bvh[p._l]._parent=c._parent;
      p._r=s._r;
      if(p._r >= 0)
        bvh[p._r]._parent=c._parent;
      //remove sibling
      s._l=s._r=s._parent=-1;
      s._cell=_verbose;
      addEmpty(sid);
      //remove current
      c._l=c._r=c._parent=-1;
      c._cell=_verbose;
      addEmpty(i);
    }
  }
  void insertLeaf(T val,sizeType i,sizeType reserveBatch=256)
  {
    if(_bvh.empty() || !hasEmpty(2)) {
      reserveEmpty(reserveBatch);
      i+=reserveBatch;
    }
    TV& bvh=const_cast<TV&>(_bvh);
    if(i == (sizeType)bvh.size()) {
      Node<T,BBOX>& c=bvh[getEmpty()];
      c._cell=val;
      c._nrCell=1;
      c._parent=-1;
      return;
    }

    Node<T,BBOX>& c=bvh[i];
    ASSERT(c._cell != _verbose && val != _verbose)
    //find parent
    c._l=getEmpty();
    c._r=getEmpty();
    Node<T,BBOX>& cl=bvh[c._l];
    cl._l=cl._r=-1;
    cl._parent=i;
    cl._cell=val;
    cl._nrCell=1;
    Node<T,BBOX>& cr=bvh[c._r];
    cr._l=cr._r=-1;
    cr._parent=i;
    cr._cell=c._cell;
    cr._nrCell=c._nrCell;
    //recompute
    c._cell=_verbose;
    c._nrCell=cl._nrCell+cr._nrCell;
    recomputeNrCell(i);
    updateBVH(i);
  }
  sizeType findLeaf(const BBOX& bb) const {
    if(_bvh.empty())
      return 0;
    sizeType curr=(sizeType)_bvh.size()-1;
    while(_bvh[curr]._cell == _verbose) {
#define SURFACE(BB) (_dim==2) ? SurfaceArea<2>::area(BB) : SurfaceArea<3>::area(BB)
      BBOX bbl=_bvh[_bvh[curr]._l]._bb;
      BBOX bbr=_bvh[_bvh[curr]._r]._bb;
      scalar costL=SURFACE(bbr);
      scalar costR=SURFACE(bbl);
      bbl.setUnion(bb);
      bbr.setUnion(bb);
      costL+=SURFACE(bbl);
      costR+=SURFACE(bbr);
      if(costL < costR)
        curr=_bvh[curr]._l;
      else curr=_bvh[curr]._r;
#undef SURFACE
    }
    return curr;
  }
  //check
  sizeType parityCheck() const {
    if(_bvh.empty())
      return 0;
    ASSERT(_bvh.back()._parent == -1)
    sizeType nrLeaf=parityCheck((sizeType)_bvh.size()-1);
    ASSERT(_bvh.back()._nrCell*2-1+nrEmpty()==(sizeType)_bvh.size())
    return nrLeaf;
  }
  sizeType parityCheck(sizeType i) const {
    sizeType ret=(_bvh[i]._cell != _verbose) ? 1 : 0;
    if(_bvh[i]._l >= 0) {
      ASSERT(_bvh[_bvh[i]._l]._parent == i && _bvh[i]._cell == _verbose)
      ASSERT(_bvh[i]._nrCell == _bvh[_bvh[i]._l]._nrCell+_bvh[_bvh[i]._r]._nrCell)
      ret+=parityCheck(_bvh[i]._l);
    } else {
      ASSERT(_bvh[i]._cell != _verbose)
    }
    if(_bvh[i]._r >= 0) {
      ASSERT(_bvh[_bvh[i]._r]._parent == i && _bvh[i]._cell == _verbose)
      ASSERT(_bvh[i]._nrCell == _bvh[_bvh[i]._l]._nrCell+_bvh[_bvh[i]._r]._nrCell)
      ret+=parityCheck(_bvh[i]._r);
    } else {
      ASSERT(_bvh[i]._cell != _verbose)
    }
    return ret;
  }
  //query
  sizeType depth() const {
    return _bvh.empty()?0:depth((sizeType)_bvh.size()-1);
  }
  sizeType depth(sizeType id) const {
    if(_bvh[id]._cell != _verbose)
      return 1;
    else return std::max(depth(_bvh[id]._l),depth(_bvh[id]._r))+1;
  }
  template <typename GCALLBACK,typename PT=Vec3>
  void pointQuery(const PT& pos,scalar eps,GCALLBACK& cb,scalar& sqrDist) const {
    if(_bvh.empty())
      return;
    std::stack<sizeType> ss;
    ss.push((sizeType)_bvh.size()-1);
    while(!ss.empty()) {
      const Node<T,BBOX>& node=_bvh[ss.top()];
      ss.pop();
      if(node._bb.distTo(pos,_dim) > eps)
        continue;
      if(node._cell == _verbose) {
        ss.push(node._l);
        ss.push(node._r);
      } else {
        cb.updateDist(node,sqrDist);
        if(sqrDist == 0.0f)return;
      }
    }
  }
  template <typename GCALLBACK>
  void pointQuery(GCALLBACK& cb) const {
    if(_bvh.empty())
      return;
    std::stack<sizeType> ss;
    ss.push((sizeType)_bvh.size()-1);
    while(!ss.empty()) {
      sizeType id=ss.top();
      const Node<T,BBOX>& node=_bvh[id];
      ss.pop();
      if(!cb.validNode(node))
        continue;
      if(node._cell == _verbose) {
        ss.push(node._l);
        ss.push(node._r);
      } else {
        cb.updateDist(node);
      }
    }
  }
  template <typename GCALLBACK>
  void pointQueryId(GCALLBACK& cb) const {
    if(_bvh.empty())
      return;
    std::stack<sizeType> ss;
    ss.push((sizeType)_bvh.size()-1);
    while(!ss.empty()) {
      sizeType id=ss.top();
      const Node<T,BBOX>& node=_bvh[id];
      ss.pop();
      if(!cb.validNode(node))
        continue;
      if(node._cell == _verbose) {
        ss.push(node._l);
        ss.push(node._r);
      } else {
        cb.updateDist(node,id);
      }
    }
  }
  template <typename GCALLBACK,typename PT=Vec3>
  void pointDistQuery(const PT& pt,GCALLBACK& cb,PT& cp,PT& n,scalar& dist,scalar* minDist=NULL) const {
    if(_bvh.empty())
      return;
    std::stack<sizeType> ss;
    ss.push((sizeType)_bvh.size()-1);
    while(!ss.empty()) {
      const Node<T,BBOX>& node=_bvh[ss.top()];
      ss.pop();
      if(node._bb.distTo(pt,_dim) > std::min<scalar>(cb.depth(),dist))
        continue;
      if(node._cell == _verbose) {
        ss.push(node._l);
        ss.push(node._r);
      } else {
        cb.updateDist(node,pt,cp,n,dist,minDist);
      }
    }
  }
  template <typename T2,typename GCALLBACK,typename TV2>
  void interBodyQuery(const BVHQuery<T2,BBOX,TV2>& query2,GCALLBACK& cb) const {
    if(_bvh.empty() || query2._bvh.empty())
      return;
    std::stack<Vec2i> ss;
    ss.push(Vec2i(_bvh.size()-1,query2._bvh.size()-1));
    while(!ss.empty()) {
      sizeType rA=ss.top()[0];
      sizeType rB=ss.top()[1];
      ss.pop();

      if(_active && query2._active && !(*_active)[rA] && !(*(query2._active))[rB])continue;
      if(!_bvh[rA]._bb.intersect(query2._bvh[rB]._bb,_dim))continue;
      if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell != query2._verbose) {
        cb.onCell(_bvh[rA],query2._bvh[rB]);
      } else if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell == query2._verbose) {
        ss.push(Vec2i(rA,query2._bvh[rB]._l));
        ss.push(Vec2i(rA,query2._bvh[rB]._r));
      } else if(_bvh[rA]._cell == _verbose && query2._bvh[rB]._cell != query2._verbose) {
        ss.push(Vec2i(_bvh[rA]._l,rB));
        ss.push(Vec2i(_bvh[rA]._r,rB));
      } else {
        ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._l));
        ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._r));
        ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._l));
        ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._r));
      }
    }
  }
  template <typename T2,typename GCALLBACK,typename TV2>
  void interBodyDepthQuery(const BVHQuery<T2,BBOX,TV2>& query2,GCALLBACK& cb,scalar depth) const {
    if(_bvh.empty() || query2._bvh.empty())
      return;
    std::stack<Vec2i> ss;
    ss.push(Vec2i(_bvh.size()-1,query2._bvh.size()-1));
    while(!ss.empty()) {
      sizeType rA=ss.top()[0];
      sizeType rB=ss.top()[1];
      ss.pop();

      if(_bvh[rA]._bb.distTo(query2._bvh[rB]._bb,_dim) > depth)continue;
      if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell != query2._verbose) {
        cb.onCell(_bvh[rA],query2._bvh[rB]);
      } else if(_bvh[rA]._cell != _verbose && query2._bvh[rB]._cell == query2._verbose) {
        ss.push(Vec2i(rA,query2._bvh[rB]._l));
        ss.push(Vec2i(rA,query2._bvh[rB]._r));
      } else if(_bvh[rA]._cell == _verbose && query2._bvh[rB]._cell != query2._verbose) {
        ss.push(Vec2i(_bvh[rA]._l,rB));
        ss.push(Vec2i(_bvh[rA]._r,rB));
      } else {
        ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._l));
        ss.push(Vec2i(_bvh[rA]._l,query2._bvh[rB]._r));
        ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._l));
        ss.push(Vec2i(_bvh[rA]._r,query2._bvh[rB]._r));
      }
    }
  }
  template <typename T2,typename CACHE,typename TV2>
  void broadphaseQuery(const BVHQuery<T2,BBOX,TV2>& query2,std::vector<CACHE>& cache) const {
    AddCache<T2,CACHE> cb(cache);
    interBodyQuery(query2,cb);
  }
  const std::vector<char>* _active;
  const TV& _bvh;
  sizeType _dim;
  T _verbose;
private:
  //dynamic update
  void reconnect(sizeType id) {
    TV& bvh=const_cast<TV&>(_bvh);
    bvh[bvh[id]._l]._parent=id;
    bvh[bvh[id]._r]._parent=id;
    bvh[id]._bb=bvh[bvh[id]._l]._bb.getUnion(bvh[bvh[id]._r]._bb);
    bvh[id]._nrCell=bvh[bvh[id]._l]._nrCell+bvh[bvh[id]._r]._nrCell;
  }
  scalar surfaceArea(sizeType bid) const {
    if(_dim==2)
      return SurfaceArea<2>::area(_bvh[bid]._bb);
    else return SurfaceArea<3>::area(_bvh[bid]._bb);
  }
  scalar surfaceArea(sizeType l,sizeType r) const {
    if(_dim==2)
      return SurfaceArea<2>::area(_bvh[l]._bb.getUnion(_bvh[r]._bb));
    else return SurfaceArea<3>::area(_bvh[l]._bb.getUnion(_bvh[r]._bb));
  }
  bool tryRotate(sizeType pid,sizeType nlid)
  {
    TV& bvh=const_cast<TV&>(_bvh);
    sizeType& otherId=bvh[pid]._l==nlid?bvh[pid]._r:bvh[pid]._l;
    //compute cost
    scalar cost=surfaceArea(nlid)+surfaceArea(otherId);
    scalar cost1=surfaceArea(bvh[nlid]._r)+surfaceArea(bvh[nlid]._l,otherId);
    scalar cost2=surfaceArea(bvh[nlid]._l)+surfaceArea(bvh[nlid]._r,otherId);
    //rotate
    if(cost1<cost && cost1<cost2)
      std::swap(bvh[nlid]._r,otherId);
    else if(cost2<cost && cost2<cost1)
      std::swap(bvh[nlid]._l,otherId);
    else return false;
    //reconnect
    reconnect(nlid);
    reconnect(pid);
    //scalar costFinal=surfaceArea(bvh[pid]._l)+surfaceArea(bvh[pid]._r);
    return true;
  }
  bool tryRecombine(sizeType pid)
  {
    TV& bvh=const_cast<TV&>(_bvh);
    sizeType lid=bvh[pid]._l;
    sizeType rid=bvh[pid]._r;
    //compute cost
    scalar cost=surfaceArea(lid)+surfaceArea(rid);
    scalar cost1=surfaceArea(bvh[lid]._l,bvh[rid]._l)+surfaceArea(bvh[lid]._r,bvh[rid]._r);
    scalar cost2=surfaceArea(bvh[lid]._l,bvh[rid]._r)+surfaceArea(bvh[lid]._r,bvh[rid]._l);
    if(cost1<cost && cost1<cost2)
      std::swap(bvh[lid]._r,bvh[rid]._l);
    else if(cost2<cost && cost2<cost1)
      std::swap(bvh[lid]._r,bvh[rid]._r);
    else return false;
    //reconnect
    reconnect(lid);
    reconnect(rid);
    //scalar costFinal=surfaceArea(bvh[pid]._l)+surfaceArea(bvh[pid]._r);
    return true;
  }
  void recomputeNrCell(sizeType pid) {
    TV& bvh=const_cast<TV&>(_bvh);
    while((pid=bvh[pid]._parent)>=0)
      bvh[pid]._nrCell=bvh[bvh[pid]._l]._nrCell+bvh[bvh[pid]._r]._nrCell;
  }
  //garbage collection
  void reserveEmpty(sizeType nr)
  {
    TV& bvh=const_cast<TV&>(_bvh);
    sizeType sz=(sizeType)bvh.size();
    bvh.resize(sz+nr);
    //reserve
    std::copy_backward(bvh.begin(),bvh.begin()+sz,bvh.end());
    for(sizeType i=nr; i<nr+sz; i++) {
      if(bvh[i]._l >= 0)
        bvh[i]._l+=nr;
      if(bvh[i]._r >= 0)
        bvh[i]._r+=nr;
      if(bvh[i]._parent >= 0)
        bvh[i]._parent+=nr;
    }
    //reconnect
    for(sizeType i=0; i<nr; i++) {
      bvh[i]._l=bvh[i]._r=bvh[i]._parent=-1;
      bvh[i]._cell=_verbose;
      bvh[i]._parent=i+1;
    }
    if((sizeType)bvh.size()>nr && !isEmpty(nr))
      bvh[nr-1]._parent=-1;
  }
  sizeType getEmpty() {
    ASSERT(isEmpty(0))
    TV& bvh=const_cast<TV&>(_bvh);
    if(isEmpty((sizeType)bvh.size()-1)) {
      bvh[bvh.size()-2]._parent=-1;
      return (sizeType)bvh.size()-1;
    } else if(bvh[0]._parent>=0) {
      sizeType ret=bvh[0]._parent;
      bvh[0]._parent=bvh[ret]._parent;
      return ret;
    } else return 0;
  }
  sizeType nrEmpty() const
  {
    if(_bvh.empty() || !isEmpty(0))
      return 0;
    sizeType curr=0,nr=1;
    while(_bvh[curr]._parent>=0) {
      curr=_bvh[curr]._parent;
      nr++;
    }
    return nr;
  }
  bool isEmpty(sizeType id) const
  {
    return _bvh[id]._cell ==_verbose && _bvh[id]._l==-1 && _bvh[id]._r==-1;
  }
  bool hasEmpty(sizeType nr0) const {
    if(_bvh.empty() || !isEmpty(0))
      return 0>=nr0;
    sizeType curr=0,nr=1;
    if(nr>=nr0)
      return true;
    while(_bvh[curr]._parent>=0) {
      curr=_bvh[curr]._parent;
      nr++;
      if(nr>=nr0)
        return true;
    }
    return false;
  }
  void addEmpty(sizeType id)
  {
    TV& bvh=const_cast<TV&>(_bvh);
    sizeType tmp=bvh[0]._parent;
    bvh[0]._parent=id;
    bvh[id]._parent=tmp;
  }
};

PRJ_END
#endif
