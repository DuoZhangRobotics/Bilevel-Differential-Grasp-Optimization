#ifndef PARTICLE_CD_H
#define PARTICLE_CD_H

#include "GridBasic.h"
#include "CollisionFunc.h"

PRJ_BEGIN

//collision detection
template <typename T,typename PS_TYPE,typename EXTRACT_POS>
class CollisionInterfaceTpl : public Serializable
{
public:
  typedef typename ScalarUtil<T>::ScalarVec3 Vec3Type;
  EIGEN_DEVICE_FUNC CollisionInterfaceTpl();
  EIGEN_DEVICE_FUNC virtual ~CollisionInterfaceTpl();
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  void reset(const BBox<T>& bb,const Vec3i& nrCell);
  void reset(const BBox<T>& bb,const Vec3Type& cellSz);
  void resize(const sizeType& nrParticle);
  virtual void prepare(const PS_TYPE& pSet);
  //fast version
  bool fill(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const;
  bool fill3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const;
  bool fill2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const;
  //brute force version
  bool fillBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const;
  bool fill3DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const;
  bool fill2DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFuncEarlyStop<typename PS_TYPE::value_type>& f) const;
  //neigh check
  bool hasNeighBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const;
  bool hasNeigh(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const;
  bool hasNeigh3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const;
  bool hasNeigh2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const;
  //version that does not return
  void fill(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const;
  void fill3D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const;
  void fill2D(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const;
  //brute force version
  void fillBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const;
  void fill3DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const;
  void fill2DBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr,CollisionFunc<typename PS_TYPE::value_type>& f) const;
  //get neigh
  std::set<sizeType> getNeighBF(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const;
  std::set<sizeType> getNeigh(const PS_TYPE& pSet,const Vec3Type& pos,const T& radSqr) const;
  //misc
  sizeType getDim() const;
  BBox<T> getBB() const;
  Vec3i getNrCell() const;
  const Grid<sizeType,T>& getStartGrid() const;
  const Grid<sizeType,T>& getEndGrid() const;
  const std::vector<sizeType,Eigen::aligned_allocator<sizeType> >& getIndex() const;
  //warpMode
  void setWarp(unsigned char D);
  void clearWarp(unsigned char D);
  const Vec3i& getWarp() const;
public:
  void warp(Vec3Type& delta,sizeType DIM) const;
  static void radixSort(const sizeType& n,sizeType* pid,sizeType* pidBK,sizeType* index,sizeType* indexBK);
  void calcHash(const sizeType& n,const PS_TYPE& pSet);
  void calcStartEnd(const sizeType& n);
  std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _pid;
  std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _pidBK;
  std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _index;
  std::vector<sizeType,Eigen::aligned_allocator<sizeType> > _indexBK;
  Grid<sizeType,T> _gridStart,_gridEnd;
  Vec3i _maxIndex;
  T _gridLen;
  //warp
  Vec3i _warpMode;
};

//debug
template <int DIM>
void debugParticleCollisionInterface(sizeType N,const Vec3i& warpMode);
template <int DIM>
void debugParticleCollisionInterfaceAllWarp(sizeType N);
void debugParticleCollisionInterfaceAllDIM(sizeType N);

template <typename T,typename PS_TYPE>
class CollisionInterface : public CollisionInterfaceTpl<T,PS_TYPE,ExtractPosParticle<typename PS_TYPE::ParticleType> > {};
typedef CollisionInterfaceTpl<scalarF,std::vector<Vec3f,Eigen::aligned_allocator<Vec3f> >,ExtractPosDirect<Vec3f> > CollisionInterfaceDirectF;
typedef CollisionInterfaceTpl<scalarD,std::vector<Vec3d,Eigen::aligned_allocator<Vec3d> >,ExtractPosDirect<Vec3d> > CollisionInterfaceDirectD;
typedef CollisionInterfaceTpl<scalar,std::vector<Vec3,Eigen::aligned_allocator<Vec3> >,ExtractPosDirect<Vec3> > CollisionInterfaceDirect;

PRJ_END

#endif
