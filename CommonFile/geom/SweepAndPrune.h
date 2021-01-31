#ifndef SWEEP_AND_PRUNE_H
#define SWEEP_AND_PRUNE_H

#include "BVHNode.h"
#include <CommonFile/Hash.h>
#include <unordered_map>
#include <unordered_set>

PRJ_BEGIN

template <typename T,int DIM>
struct SweepAndPrune
{
  struct BOXPointer {
    sizeType _boxId;
    bool _isMin;
  };
  SweepAndPrune(bool check=false):_check(check) {}
  void update(const BBox<T,DIM>& bb,sizeType i) {
    ASSERT_MSGV(i>=0 && i<(sizeType)_bbss.size(),"Invalid box index 0<=%d<%d",i,_bbss.size())
    for(sizeType d=0; d<DIM; d++) {
      //min
      _bbss[i]._minC[d]=bb._minC[d];
      adjust(d,_minIndex[d][i]);
      //max
      _bbss[i]._maxC[d]=bb._maxC[d];
      adjust(d,_maxIndex[d][i]);
    }
    parityCheck();
  }
  void add(const BBox<T,DIM>& bb) {
    sizeType i=_bbss.size();
    _bbss.push_back(bb);
    for(sizeType d=0; d<DIM; d++) {
      //min
      _minIndex[d].push_back(_minmax[d].size());
      _minmax[d].push_back({i,true});
      adjust(d,_minIndex[d][i]);
      //max
      _maxIndex[d].push_back(_minmax[d].size());
      _minmax[d].push_back({i,false});
      adjust(d,_maxIndex[d][i]);
    }
    parityCheck();
  }
  virtual void remove(sizeType i) {
    ASSERT_MSGV(i>=0 && i<(sizeType)_bbss.size(),"Invalid box index 0<=%d<%d",i,_bbss.size())
    _bbss.erase(_bbss.begin()+i);
    for(sizeType d=0; d<DIM; d++) {
      _minIndex[d].erase(_minIndex[d].begin()+i);
      _maxIndex[d].erase(_maxIndex[d].begin()+i);
      //reconnect
      sizeType k=0;
      for(sizeType j=0; j<(sizeType)_minmax[d].size(); j++)
        if(_minmax[d][j]._boxId!=i) {
          if(_minmax[d][j]._boxId>i)
            _minmax[d][j]._boxId--;
          _minmax[d][k]=_minmax[d][j];
          index(d,_minmax[d][j])=k;
          k++;
        }
      _minmax[d].resize(k);
    }
    parityCheck();
  }
  virtual void intersect(std::function<void(sizeType,sizeType)> pairs) const {
    std::unordered_map<Vec2i,sizeType,Hash> imap;
    std::unordered_set<sizeType> met;
    for(sizeType d=0; d<DIM; d++)
      for(sizeType off=0; off<(sizeType)_minmax[d].size(); off++)
        if(_minmax[d][off]._isMin) {
          for(sizeType i:met) {
            Vec2i pair(i,_minmax[d][off]._boxId);
            if(pair[0]>pair[1])
              std::swap(pair[0],pair[1]);
            imap[pair]++;
          }
          met.insert(_minmax[d][off]._boxId);
        } else {
          met.erase(_minmax[d][off]._boxId);
        }
    for(const std::pair<Vec2i,sizeType>& p:imap)
      if(p.second==DIM)
        pairs(p.first[0],p.first[1]);
  }
  void intersectBF(std::function<void(sizeType,sizeType)> pairs) const {
    for(sizeType i=0; i<(sizeType)_bbss.size(); i++)
      for(sizeType j=i+1; j<(sizeType)_bbss.size(); j++)
        if(_bbss[i].intersect(_bbss[j]))
          pairs(i,j);
  }
  void debugCompare() const {
    std::unordered_set<Vec2i,Hash> s,sBF;
    intersect([&](sizeType i,sizeType j) {
      s.insert(Vec2i(i,j));
    });
    intersectBF([&](sizeType i,sizeType j) {
      sBF.insert(Vec2i(i,j));
    });
    ASSERT(s==sBF)
  }
protected:
  virtual void adjust(sizeType d,sizeType id) {
    while(id>0 && val(d,_minmax[d][id-1])>val(d,_minmax[d][id])) {
      std::swap(index(d,_minmax[d][id-1]),index(d,_minmax[d][id]));
      std::swap(_minmax[d][id-1],_minmax[d][id]);
      id--;
    }
    while(id<(sizeType)_minmax[d].size()-1 && val(d,_minmax[d][id])>val(d,_minmax[d][id+1])) {
      std::swap(index(d,_minmax[d][id]),index(d,_minmax[d][id+1]));
      std::swap(_minmax[d][id],_minmax[d][id+1]);
      id++;
    }
  }
  sizeType& index(sizeType d,const BOXPointer& ptr) {
    if(ptr._isMin)
      return _minIndex[d][ptr._boxId];
    else return _maxIndex[d][ptr._boxId];
  }
  const sizeType& index(sizeType d,const BOXPointer& ptr) const {
    if(ptr._isMin)
      return _minIndex[d][ptr._boxId];
    else return _maxIndex[d][ptr._boxId];
  }
  T val(sizeType d,const BOXPointer& ptr) const {
    if(ptr._isMin)
      return _bbss[ptr._boxId]._minC[d];
    else return _bbss[ptr._boxId]._maxC[d];
  }
  void parityCheck() const {
    if(!_check)
      return;
    for(sizeType d=0; d<DIM; d++)
      for(sizeType i=0; i<(sizeType)_minmax[d].size(); i++) {
        ASSERT(index(d,_minmax[d][i])==i)
      }
    for(sizeType d=0; d<DIM; d++)
      for(sizeType i=0; i<(sizeType)_minmax[d].size()-1; i++) {
        ASSERT(val(d,_minmax[d][i])<=val(d,_minmax[d][i+1]))
      }
    for(sizeType d=0; d<DIM; d++)
      for(sizeType i=0; i<(sizeType)_bbss.size(); i++) {
        ASSERT(_minmax[d][_minIndex[d][i]]._boxId==i)
        ASSERT( _minmax[d][_minIndex[d][i]]._isMin)
        ASSERT(_minmax[d][_maxIndex[d][i]]._boxId==i)
        ASSERT(!_minmax[d][_maxIndex[d][i]]._isMin)
      }
  }
  //data
  std::vector<BOXPointer> _minmax[DIM];
  std::vector<sizeType> _minIndex[DIM];
  std::vector<sizeType> _maxIndex[DIM];
  std::vector<BBox<T,DIM>> _bbss;
  bool _check;
};
template <typename T,int DIM>
struct SweepAndPruneIncremental : public SweepAndPrune<T,DIM>
{
  using SweepAndPrune<T,DIM>::_minmax;
  using SweepAndPrune<T,DIM>::index;
  using SweepAndPrune<T,DIM>::val;
  SweepAndPruneIncremental(bool check=false):SweepAndPrune<T,DIM>(check) {}
  virtual void remove(sizeType i) override {
    SweepAndPrune<T,DIM>::remove(i);
    std::unordered_map<Vec2i,sizeType,Hash> imap;
    for(const std::pair<Vec2i,sizeType>& p:_imap)
      if(p.first[0]!=i && p.first[1]!=i)
        imap[Vec2i(p.first[0]>i?p.first[0]-1:p.first[0],
                   p.first[1]>i?p.first[1]-1:p.first[1])]=p.second;
    _imap.swap(imap);
  }
  virtual void intersect(std::function<void(sizeType,sizeType)> pairs) const override {
    for(const std::pair<Vec2i,sizeType>& p:_imap)
      if(p.second==DIM)
        pairs(p.first[0],p.first[1]);
  }
protected:
  virtual void adjust(sizeType d,sizeType id) override {
    while(id>0 && val(d,_minmax[d][id-1])>val(d,_minmax[d][id])) {
      if(_minmax[d][id-1]._isMin && !_minmax[d][id]._isMin) {
        Vec2i p(_minmax[d][id-1]._boxId,_minmax[d][id]._boxId);
        if(p[0]>p[1])
          std::swap(p[0],p[1]);
        std::unordered_map<Vec2i,sizeType,Hash>::iterator it=_imap.find(p);
        if(it!=_imap.end()) {
          it->second--;
          if(it->second==0)
            _imap.erase(it);
        }
      }
      if(!_minmax[d][id-1]._isMin && _minmax[d][id]._isMin) {
        Vec2i p(_minmax[d][id-1]._boxId,_minmax[d][id]._boxId);
        if(p[0]>p[1])
          std::swap(p[0],p[1]);
        _imap[p]++;
      }
      std::swap(index(d,_minmax[d][id-1]),index(d,_minmax[d][id]));
      std::swap(_minmax[d][id-1],_minmax[d][id]);
      id--;
    }
    while(id<(sizeType)_minmax[d].size()-1 && val(d,_minmax[d][id])>val(d,_minmax[d][id+1])) {
      if(_minmax[d][id]._isMin && !_minmax[d][id+1]._isMin) {
        Vec2i p(_minmax[d][id]._boxId,_minmax[d][id+1]._boxId);
        if(p[0]>p[1])
          std::swap(p[0],p[1]);
        std::unordered_map<Vec2i,sizeType,Hash>::iterator it=_imap.find(p);
        if(it!=_imap.end()) {
          it->second--;
          if(it->second==0)
            _imap.erase(it);
        }
      }
      if(!_minmax[d][id]._isMin && _minmax[d][id+1]._isMin) {
        Vec2i p(_minmax[d][id]._boxId,_minmax[d][id+1]._boxId);
        if(p[0]>p[1])
          std::swap(p[0],p[1]);
        _imap[p]++;
      }
      std::swap(index(d,_minmax[d][id]),index(d,_minmax[d][id+1]));
      std::swap(_minmax[d][id],_minmax[d][id+1]);
      id++;
    }
  }
  //data
  std::unordered_map<Vec2i,sizeType,Hash> _imap;
};

PRJ_END

#endif
