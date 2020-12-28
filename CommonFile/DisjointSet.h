#ifndef DISJOINT_SET_H
#define DISJOINT_SET_H

#include "MathBasic.h"
#include "IOBasic.h"

PRJ_BEGIN

template <typename WEIGHT>
struct DisjointSetElem : public SerializableBase {
  virtual bool read(std::istream& is,IOData* dat) {
    readBinaryData(_rank,is);
    readBinaryData(_p,is);
    readBinaryData(_size,is);
    readBinaryData(_Int,is,dat);
    return is.good();
  }
  virtual bool write(std::ostream& os,IOData* dat) const {
    writeBinaryData(_rank,os);
    writeBinaryData(_p,os);
    writeBinaryData(_size,os);
    writeBinaryData(_Int,os,dat);
    return os.good();
  }
  virtual std::shared_ptr<SerializableBase> copy() const {
    return std::shared_ptr<SerializableBase>(new DisjointSetElem<WEIGHT>);
  }
  virtual std::string type() const {
    return typeid(DisjointSetElem<WEIGHT>).name();
  }
  sizeType _rank;	//rank based merging
  sizeType _p;	//pointer to parent
  sizeType _size;	//size of current subtree, size=0 is a tag to avoid joining
  WEIGHT _Int;	//another criterion for merging
};
template <typename WEIGHT,typename ELEM=DisjointSetElem<WEIGHT> >
class DisjointSet : public SerializableBase
{
public:
  DisjointSet():_nrSet(0),_avoidId(-1) {}
  DisjointSet(sizeType elements):_avoidId(-1) {
    resize(elements);
  }
  void resize(sizeType elements) {
    _elts.clear();
    _elts.resize(elements);
    _nrSet=elements;
    for(sizeType i=0; i < elements; i++) {
      _elts[i]._rank=0;
      _elts[i]._size=1;
      _elts[i]._p=i;
      _elts[i]._Int=EigenTraits<WEIGHT>::value();
    }
  }
  virtual ~DisjointSet() {}
  virtual bool read(std::istream& is,IOData* dat) {
    readBinaryData(_elts,is,dat);
    readBinaryData(_nrSet,is,dat);
    readBinaryData(_avoidId,is,dat);
    return is.good();
  }
  virtual bool write(std::ostream& os,IOData* dat) const {
    writeBinaryData(_elts,os,dat);
    writeBinaryData(_nrSet,os,dat);
    writeBinaryData(_avoidId,os,dat);
    return os.good();
  }
  virtual std::shared_ptr<SerializableBase> copy() const {
    return std::shared_ptr<SerializableBase>(new DisjointSet<WEIGHT,ELEM>);
  }
  virtual std::string type() const {
    return typeid(DisjointSet<WEIGHT,ELEM>).name();
  }
  sizeType find(sizeType x) {
    sizeType y=x;
    bool avoidMerge=(y == _avoidId);
    while(y != _elts[y]._p) {
      y=_elts[y]._p;
      avoidMerge=avoidMerge || (y == _avoidId);
    }
    if(!avoidMerge)
      _elts[x]._p=y;
    return y;
  }
  bool avoid(sizeType x) const {
    sizeType y=x;
    bool avoidMerge=(y == _avoidId);
    while(y != _elts[y]._p) {
      y=_elts[y]._p;
      avoidMerge=avoidMerge || (y == _avoidId);
    }
    return avoidMerge;
  }
  void join(sizeType x,sizeType y) {
    if(x == y)
      return;
    if(y == _avoidId || (x != _avoidId && _elts[x]._rank > _elts[y]._rank)) {
      _elts[y]._p=x;
      _elts[x]._size+=_elts[y]._size;
    } else {
      _elts[x]._p=y;
      _elts[y]._size+=_elts[x]._size;
      if(_elts[x]._rank == _elts[y]._rank)
        _elts[y]._rank++;
    }
    _nrSet--;
  }
  void joinSafe(sizeType x,sizeType y) {
    if(_elts[x]._size == 0 || _elts[y]._size == 0)
      return;
    join(find(x),find(y));
  }
  void avoidMerge(sizeType avoidId) {
    _avoidId=avoidId;
  }
  sizeType size(sizeType x) const {
    return _elts[x]._size;
  }
  sizeType sizeR(sizeType x) {
    return _elts[find(x)]._size;
  }
  const WEIGHT& weight(sizeType x) const {
    return _elts[x]._Int;
  }
  WEIGHT& weightR(sizeType x) {
    return _elts[find(x)]._Int;
  }
  sizeType numSets() const {
    return _nrSet;
  }
  //std::vector<ELEM,std::fast_pool_allocator<ELEM> > _elts;
  std::vector<ELEM> _elts;
  sizeType _nrSet,_avoidId;
};

PRJ_END

#endif
