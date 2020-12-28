#include "ParticleCD.h"

PRJ_BEGIN

//CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::CollisionInterfaceTpl():Serializable(typeid(CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>).name()) {}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::~CollisionInterfaceTpl() {}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::read(std::istream& is,IOData* dat) {
  readBinaryData(_pid,is);
  readBinaryData(_pidBK,is);
  readBinaryData(_index,is);
  readBinaryData(_indexBK,is);
  _gridStart.read(is);
  _gridEnd.read(is);
  readBinaryData(_maxIndex,is);
  readBinaryData(_gridLen,is);
  return is.good();
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::write(std::ostream& os,IOData* dat) const {
  writeBinaryData(_pid,os);
  writeBinaryData(_pidBK,os);
  writeBinaryData(_index,os);
  writeBinaryData(_indexBK,os);
  _gridStart.write(os);
  _gridEnd.write(os);
  writeBinaryData(_maxIndex,os);
  writeBinaryData(_gridLen,os);
  return os.good();
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
std::shared_ptr<SerializableBase> CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::copy() const {
  return std::shared_ptr<SerializableBase>(new CollisionInterfaceTpl);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::reset(const BBox<T>& bb,const Vec3i& nrCell) {
  _gridStart.reset(nrCell,bb,0,false,16);
  _gridEnd=_gridStart;
  _maxIndex=_gridStart.getNrCell()-Vec3i::Ones();
  if(_gridStart.getDim() == 2)
    _gridLen=std::min<T>(_gridStart.getCellSize().x(),_gridStart.getCellSize().y());
  else
    _gridLen=_gridStart.getCellSize().minCoeff();
  _warpMode=Vec3i::Zero();
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::reset(const BBox<T>& bb,const Vec3Type& cellSz) {
  Vec3i nrCell=ceilV((Vec3Type)(bb.getExtent().array()/cellSz.array()).matrix());
  reset(bb,nrCell);
  _warpMode=Vec3i::Zero();
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::resize(const sizeType& nrParticle) {
  _pid.resize(nrParticle);
  _pidBK.resize(nrParticle);
  _index.resize(nrParticle);
  _indexBK.resize(nrParticle);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::prepare(const PS_TYPE& pSet) {
  const sizeType n=pSet.size();
  if(n == 0)
    return;
  if((sizeType)_pid.size() != n)
    resize(pSet.size());

  calcHash(n,pSet);
  radixSort(n,&(_pid[0]),&(_pidBK[0]),&(_index[0]),&(_indexBK[0]));
  calcStartEnd(n);
}
//fast version
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const {
  Vec3Type radVec=Vec3Type::Zero();
  sizeType DIM=_gridStart.getDim();
  radVec.segment(0,DIM).setConstant(std::sqrt(radSqr));
  Vec3Type pos0=pos-radVec;
  Vec3Type pos1=pos+radVec;
  //warp compensate
  const Vec3i& nrP=_gridStart.getNrPoint();
  const Vec3i& stride=_gridStart.getStride();
  Vec3i warpCountMax=Vec3i::Zero();
  Vec3Type warpRange;
  for(sizeType d=0; d<DIM; d++)
    if(_warpMode[d] > 0) {
      //adjust bound
      T pos0Old=pos0[d];
      warpRange[d]=handleScalarWarp<T>(pos0[d],_gridStart.getBB()._minC[d],_gridStart.getBB()._maxC[d],NULL);
      pos1[d]+=pos0[d]-pos0Old;
      //adjust maxC
      handleScalarWarp<T>(pos1[d],_gridStart.getBB()._minC[d],_gridStart.getBB()._maxC[d],&warpCountMax[d]);
      ASSERT(warpCountMax[d] <= 0)
      warpCountMax[d]*=-1;
    }
  //increase one last cell
  Vec3i id0=floorV(_gridStart.getIndexFracSafe(pos0));
  Vec3i id1=ceilV(_gridStart.getIndexFracSafe(pos1));
  for(sizeType d=0; d<DIM; d++) {
    if(id1[d] > nrP[d]) {
      if(_warpMode[d] > 0) {
        id1[d]-=nrP[d];
        warpCountMax[d]++;
      } else id1[d]=std::min<sizeType>(id1[d],nrP[d]);
    }
    if(warpCountMax[d] == 0 && id1[d] <= id0[d])
      return false;
  }
  //setup iterator
  //beg
  Vec3i idBeg=id0;
  Vec3i warpCountBeg=Vec3i::Zero();
  //end
  Vec3i idEnd=id0;
  idEnd[0]=id1[0];
  Vec3i warpCountEnd=Vec3i::Zero();
  warpCountEnd[0]=warpCountMax[0];
  sizeType off=idBeg.dot(stride);
  //run iterator loop
  Vec3Type delta,delta2;
  while(!(idBeg == idEnd && warpCountBeg == warpCountEnd)) {
    //add neigh
    std::pair<sizeType,sizeType> range(_gridStart.get(off),_gridEnd.get(off));
    for(sizeType k=range.first; k < range.second; k++) {
      delta=EXTRACT_POS::extract(pSet[_index[k]])-pos;
      warp(delta2=delta,DIM);
      if(delta2.squaredNorm() < radSqr)
        if(f(EXTRACT_POS::add(pSet[_index[k]],delta2-delta),_index[k]))
          return true;
    }
    //advance iterator
    for(sizeType d=DIM-1; d>=0; d--) {
      idBeg[d]++;
      off+=stride[d];
      if(idBeg[d] == nrP[d] && warpCountBeg[d] < warpCountMax[d]) {
        idBeg[d]-=nrP[d];
        off-=nrP[d]*stride[d];
        warpCountBeg[d]++;
      }
      if(d == 0)
        break;
      if(idBeg[d] == id1[d] && warpCountBeg[d] == warpCountMax[d]) {
        idBeg[d]=id0[d];
        off+=(id0[d]-id1[d])*stride[d];
        warpCountBeg[d]=0;
      } else break;
    }
  }
  return false;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const {
  return fill(pSet,pos,radSqr,f);  //deprecated
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const {
  return fill(pSet,pos,radSqr,f);  //deprecated
}
//brute force version
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fillBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const {
  const sizeType n=pSet.size();
  sizeType DIM=_gridStart.getDim();
  Vec3Type delta,delta2;
  for(sizeType k=0; k < n; k++) {
    delta=EXTRACT_POS::extract(pSet[k])-pos;
    warp(delta2=delta,DIM);
    if(delta2.squaredNorm() < radSqr)
      if(f(EXTRACT_POS::assign(pSet[_index[k]],delta2-delta),k))
        return true;
  }
  return false;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill3DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const {
  return fillBF(pSet,pos,radSqr,f);  //deprecated
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill2DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const {
  return fillBF(pSet,pos,radSqr,f);  //deprecated
}
//neigh check
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::hasNeighBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
  CollisionFuncAlwaysTrue<typename PS_TYPE::value_type> f;
  return fillBF(pSet,pos,radSqr,f);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::hasNeigh(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
  CollisionFuncAlwaysTrue<typename PS_TYPE::value_type> f;
  return fill(pSet,pos,radSqr,f);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::hasNeigh3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
  return hasNeigh(pSet,pos,radSqr);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
bool CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::hasNeigh2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
  return hasNeigh(pSet,pos,radSqr);
}
//version that does not return
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
  CollisionFuncWrapper<typename PS_TYPE::value_type> wrapper(f);
  fill(pSet,pos,radSqr,wrapper);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
  fill(pSet,pos,radSqr,f);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
  fill(pSet,pos,radSqr,f);
}
//brute force version
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fillBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
  CollisionFuncWrapper<typename PS_TYPE::value_type> wrapper(f);
  fillBF(pSet,pos,radSqr,wrapper);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill3DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
  fillBF(pSet,pos,radSqr,f);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::fill2DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const {
  fillBF(pSet,pos,radSqr,f);
}
//get neigh
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
std::set<sizeType> CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getNeighBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
  CollisionFuncRecord<typename PS_TYPE::value_type> f;
  fillBF(pSet,pos,radSqr,f);
  return f._neigh;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
std::set<sizeType> CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getNeigh(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const {
  CollisionFuncRecord<typename PS_TYPE::value_type> f;
  fill(pSet,pos,radSqr,f);
  return f._neigh;
}
//misc
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
sizeType CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getDim() const {
  return _gridStart.getDim();
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
BBox<T> CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getBB() const {
  return _gridStart.getBB();
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
Vec3i CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getNrCell() const {
  return _gridStart.getNrCell();
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
const Grid<sizeType,T>& CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getStartGrid() const {
  return _gridStart;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
const Grid<sizeType,T>& CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getEndGrid() const {
  return _gridEnd;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
const std::vector<sizeType,Eigen::aligned_allocator<sizeType> >& CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getIndex() const {
  return _index;
}
//warpMode
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::setWarp(unsigned char D) {
  _warpMode[D]=1;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::clearWarp(unsigned char D) {
  _warpMode[D]=0;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
const Vec3i& CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::getWarp() const {
  return _warpMode;
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::warp(Vec3Type& delta,sizeType DIM) const {
  for(sizeType d=0; d<DIM; d++)
    if(_warpMode[d] > 0)
      handleScalarWarpDiff(delta[d],_gridStart.getBB()._minC[d],_gridStart.getBB()._maxC[d],NULL);
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::radixSort(const sizeType& n,sizeType* pid,sizeType* pidBK,sizeType* index,sizeType* indexBK) {
  sizeType bucket[sizeof(sizeType)][256];
  memset(bucket,0,sizeof(sizeType)*256*sizeof(sizeType));

  //count number of values in a bucket
  for(sizeType i=0; i<n; i++) {
    for(size_t p=0; p<sizeof(sizeType); p++)
      bucket[p][(pid[i]>>(p*8))&255]++;
    index[i]=i;
  }

  //accumulate buckets
  for(size_t p=0; p<sizeof(sizeType); p++)
    for(sizeType i=1; i<256; i++)
      bucket[p][i]+=bucket[p][i-1];
  for(size_t p=0; p<sizeof(sizeType); p++)
    for(sizeType i=255; i>=1; i--)
      bucket[p][i]=bucket[p][i-1];
  for(size_t j=0; j<sizeof(sizeType); j++)
    bucket[j][0]=0;

  //redistribute
  for(size_t p=0; p<sizeof(sizeType); p++) {
    for(sizeType i=0; i<n; i++) {
      sizeType& pos=bucket[p][(pid[i]>>(p*8))&255];
      indexBK[pos]=index[i];
      pidBK[pos]=pid[i];
      pos++;
    }
    std::swap(pid,pidBK);
    std::swap(index,indexBK);
  }
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::calcHash(const sizeType& n,const PS_TYPE& pSet) {
  const sizeType maxVal=_gridStart.getSzLinear();
  sizeType DIM=_gridStart.getDim();
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<n; i++) {
    Vec3Type pos=EXTRACT_POS::extract(pSet[i]);
    for(sizeType d=0; d<DIM; d++)
      if(_warpMode[d] > 0)
        handleScalarWarp<T>(pos[d],_gridStart.getBB()._minC[d],_gridStart.getBB()._maxC[d],NULL);
    Vec3i base=floorV(_gridStart.getIndexFracSafe(pos));
    _pid[i]=_gridStart.getIndex(base);
    ASSERT(_pid[i] >= 0 && _pid[i] <= maxVal)
  }
}
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
void CollisionInterfaceTpl<T,PS_TYPE,EXTRACT_POS>::calcStartEnd(const sizeType& n) {
  const sizeType maxVal=_gridStart.getSzLinear();
  _gridStart.init(0);
  _gridEnd.init(0);
  OMP_PARALLEL_FOR_
  for(sizeType i=1; i<n; i++) {
    if(_pid[i] != _pid[i-1]) {
      if(_pid[i] != maxVal)
        _gridStart.get(_pid[i])=i;
      _gridEnd.get(_pid[i-1])=i;
    }
  }

  if(_pid[0] != maxVal)
    _gridStart.get(_pid[0])=0;
  if(_pid[n-1] != maxVal)
    _gridEnd.get(_pid[n-1])=n;
}

//debug
template <int DIM>
void debugParticleCollisionInterface(sizeType N,const Vec3i& warpMode)
{
  BBox<scalar> bb;
  CollisionInterfaceDirect cd;
  std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
  for(sizeType i=0; i<N; i++) {
    vss.push_back(Vec3::Zero());
    vss.back().segment<DIM>(0).setRandom();
    bb.setUnion(vss.back());
  }
  bb.enlarged(0.01f,DIM);
  scalar rad=0.1f;
  cd.reset(bb,(Vec3)Vec3::Constant(rad));
  cd.prepare(vss);
  for(sizeType d=0; d<DIM; d++)
    if(warpMode[d] > 0)
      cd.setWarp((unsigned char)d);
  OMP_PARALLEL_FOR_
  for(sizeType i=0; i<N; i++) {
    std::set<sizeType> neigh=cd.getNeigh(vss,vss[i],rad*rad);
    std::set<sizeType> neighBF=cd.getNeighBF(vss,vss[i],rad*rad);
    ASSERT(neigh == neighBF && cd.hasNeigh(vss,vss[i],rad*rad) == cd.hasNeighBF(vss,vss[i],rad*rad))
  }
}
template <int DIM>
void debugParticleCollisionInterfaceAllWarp(sizeType N)
{
  INFO("Debugging ParticleCD all warp!")
  debugParticleCollisionInterface<DIM>(N,Vec3i::Zero());
  for(sizeType d=0; d<DIM; d++)
    debugParticleCollisionInterface<DIM>(N,Vec3i::Unit(d));
}
void debugParticleCollisionInterfaceAllDIM(sizeType N)
{
  INFO("Debugging ParticleCD all DIM!")
  debugParticleCollisionInterfaceAllWarp<2>(N);
  debugParticleCollisionInterfaceAllWarp<3>(N);
}

//instance
template class CollisionInterfaceTpl<scalarF,ParticleSetF,ExtractPosParticle<ParticleSetF::ParticleType> >;
template class CollisionInterfaceTpl<scalarD,ParticleSetD,ExtractPosParticle<ParticleSetD::ParticleType> >;
template class CollisionInterfaceTpl<scalarF,std::vector<Vec3f,Eigen::aligned_allocator<Vec3f> >,ExtractPosDirect<Vec3f> >;
template class CollisionInterfaceTpl<scalarD,std::vector<Vec3d,Eigen::aligned_allocator<Vec3d> >,ExtractPosDirect<Vec3d> >;

PRJ_END
