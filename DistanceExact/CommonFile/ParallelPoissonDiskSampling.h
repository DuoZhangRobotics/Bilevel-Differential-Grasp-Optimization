#ifndef PARALLEL_POISSON_DISK_SAMPLING_H
#define PARALLEL_POISSON_DISK_SAMPLING_H

#include "ParticleSet.h"
#include "GridBasic.h"
#include "ObjMesh.h"

PRJ_BEGIN

template <typename T> class ObjMeshTpl;

class ParallelPoissonDiskSampling
{
public:
  struct HashEntry {
    HashEntry():_key(-1),_from(0),_to(0),_sample(-1) {}
    sizeType _key;
    sizeType _from;
    sizeType _to;
    sizeType _sample;
  };
  struct PhaseInfo {
    PhaseInfo() {}
    PhaseInfo(const Vec3i& id):_id(id),_from(0),_to(0) {}
    Vec3i _id;
    sizeType _from;
    sizeType _to;
  };
  ParallelPoissonDiskSampling(sizeType dim);
  virtual ~ParallelPoissonDiskSampling() {}
  void sample();
  void sample(ObjMesh& mesh);
  void sample(const Vec3& ctr,scalar rad);
  void sample(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss,
              const std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> >& iss,sizeType* nrP=NULL);
  void sample(const ScalarField& field,bool vol=false);
  void poissonFilter(std::vector<Vec3d,Eigen::aligned_allocator<Vec3d> >& obj_p,std::vector<Vec3d,Eigen::aligned_allocator<Vec3d> >& obj_n);
  scalar getDensity() const {
    return _density;
  }
  scalar& getDensity() {
    return _density;
  }
  void setDensity(scalar d) {
    _density = d;
  }
  scalar getRadius() const {
    return _radius;
  }
  scalar& getRadius() {
    return _radius;
  }
  void setRadius(scalar r) {
    _radius = r;
  }
  ParticleSetN& getPSet() {
    return _PSet;
  }
  const ParticleSetN& getPSet() const {
    return _PSet;
  }
  const ParticleSetN& getRawPSet() const {
    return _rawPSet;
  }
  void buildPhase(sizeType phaseSz);
  void generateRawSample(ObjMesh& mesh);
  void generateRawSample(const std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& vss,
                         const std::vector<Vec4i,Eigen::aligned_allocator<Vec4i> >& iss,sizeType* nrP);
  void generateRawSample2D(const ScalarField& phi,bool vol);
  void generateRawSample3D(const ScalarField& phi,bool vol);
protected:
  void generateHashTable();
  void generatePoissonDisk();
  void parityCheck();
  bool generateRawPSet(Vec4& bary,sizeType& id,const std::vector<scalar>& poss,const sizeType& nrT,sizeType& seed) const;
  bool generateRawPSet(Vec3& bary,sizeType& id,const std::vector<scalar>& poss,const sizeType& nrT,sizeType& seed) const;
  sizeType calcKey(const Vec3& pos) const;
  sizeType calcPhase(const Vec3& pos) const;
  sizeType findEntry(const sizeType& key) const;
  sizeType insertEntry(const sizeType& key);
  scalar dist(const ParticleN<scalar>& a,const ParticleN<scalar>& b) const;
protected:
  //samples
  ParticleSetN _rawPSet;
  ParticleSetN _PSet;
  //spatial data structure
  std::shared_ptr<sizeType> _key;
  std::shared_ptr<sizeType> _index;
  std::vector<HashEntry> _hashTable;
  sizeType _szTable;
  scalar _cellSz;
  BBox<scalar> _bb;
  Vec3i _stride;
  //param
  scalar _density;
  scalar _radius;
  sizeType _dim;
  sizeType _ratioSz;
  sizeType _bucketSz;
  sizeType _k;
  bool _euclid;
  //phase
  sizeType _phaseSz;
  Grid<sizeType,scalar> _phase;
  std::vector<PhaseInfo> _phaseMap;
  std::vector<sizeType> _phaseGroup;
};

PRJ_END

#endif
