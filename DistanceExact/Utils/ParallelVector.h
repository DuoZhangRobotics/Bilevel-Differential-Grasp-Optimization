#ifndef PARALLEL_VECTOR_H
#define PARALLEL_VECTOR_H

#include <CommonFile/MathBasic.h>

PRJ_BEGIN

template <typename T>
class ParallelVector
{
public:
  typedef std::vector<T,Eigen::aligned_allocator<T> > vector_type;
  typedef typename std::vector<T,Eigen::aligned_allocator<T> >::const_iterator const_iterator;
  typedef typename std::vector<T,Eigen::aligned_allocator<T> >::iterator iterator;
  ParallelVector() {
    clear();
  }
  void clear() {
    clear(OmpSettings::getOmpSettings().nrThreads());
  }
  void clear(sizeType nr) {
    _blocks.assign(nr,vector_type());
  }
  void push_back(const T& newVal) {
    _blocks[id()].push_back(newVal);
  }
  template <typename IT>
  void insert(IT beg,IT end) {
    vector_type& v=_blocks[id()];
    v.insert(v.end(),beg,end);
  }
  const_iterator begin() const {
    const_cast<ParallelVector<T>&>(*this).join();
    return _blocks[0].begin();
  }
  const_iterator end() const {
    const_cast<ParallelVector<T>&>(*this).join();
    return _blocks[0].end();
  }
  iterator begin() {
    join();
    return _blocks[0].begin();
  }
  iterator end() {
    join();
    return _blocks[0].end();
  }
  const vector_type& getVector() const {
    const_cast<ParallelVector<T>*>(this)->join();
    return _blocks[0];
  }
  vector_type& getVector() {
    join();
    return _blocks[0];
  }
protected:
  sizeType id() const {
    return OmpSettings::getOmpSettings().threadId()%(sizeType)_blocks.size();
  }
  void join() {
    for(sizeType i=1; i<(sizeType)_blocks.size(); i++) {
      _blocks[0].insert(_blocks[0].end(),_blocks[i].begin(),_blocks[i].end());
      _blocks[i].clear();
    }
  }
  std::vector<vector_type> _blocks;
};
template <typename T>
class ParallelMatrix : public ParallelVector<T>
{
public:
  typedef std::vector<T,Eigen::aligned_allocator<T> > vector_type;
  typedef typename std::vector<T,Eigen::aligned_allocator<T> >::const_iterator const_iterator;
  typedef typename std::vector<T,Eigen::aligned_allocator<T> >::iterator iterator;
  ParallelMatrix() {}
  ParallelMatrix(T example) {
    assign(OmpSettings::getOmpSettings().nrThreads(),example);
  }
  ParallelMatrix(sizeType nr,T example) {
    assign(nr,example);
  }
  void assign(T example) {
    assign(OmpSettings::getOmpSettings().nrThreads(),example);
  }
  void assign(sizeType nr,T example) {
    _blocks.assign(nr,example);
  }
  void clear() {
    _blocks[0].setZero();
    assign((sizeType)_blocks.size(),_blocks[0]);
  }
  template <typename TOTHER>
  ParallelMatrix<T>& operator+=(const TOTHER& other) {
    _blocks[id()]+=other;
    return *this;
  }
  const T& getMatrixI() const {
    return _blocks[id()];
  }
  T& getMatrixI() {
    return _blocks[id()];
  }
  const T& getMatrix() const {
    const_cast<ParallelVector<T>*>(this)->join();
    return _blocks[0];
  }
  T& getMatrix() {
    join();
    return _blocks[0];
  }
protected:
  sizeType id() const {
    return OmpSettings::getOmpSettings().threadId()%(sizeType)_blocks.size();
  }
  void join() {
    for(sizeType i=1; i<(sizeType)_blocks.size(); i++) {
      _blocks[0]+=_blocks[i];
      _blocks[i].setZero();
    }
  }
  std::vector<T> _blocks;
};

PRJ_END

#endif
