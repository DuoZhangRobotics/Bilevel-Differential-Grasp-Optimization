#ifndef PARTICLE_SET_H
#define PARTICLE_SET_H

#include "Config.h"
#include "MathBasic.h"
#include "IO.h"
#include <type_traits>

PRJ_BEGIN

template <typename P_TYPE>
struct ParticleSetTpl;
template <typename T>
struct HasPos : std::false_type {};
template <typename T>
struct HasVel : std::false_type {
  static void addData(VTKWriter<typename T::scalarType>& writer,const ParticleSetTpl<T>& pset) {}
};
template <typename T>
struct HasNormal : std::false_type {};
template <typename T>
struct HasDensity : std::false_type {
  static void addData(VTKWriter<typename T::scalarType>& writer,const ParticleSetTpl<T>& pset) {}
};
template <typename T>
struct HasPressure : std::false_type {
  static void addData(VTKWriter<typename T::scalarType>& writer,const ParticleSetTpl<T>& pset) {}
};
template <typename T>
struct HasWeight : std::false_type {
  static void addData(VTKWriter<typename T::scalarType>& writer,const ParticleSetTpl<T>& pset) {}
};
template <typename T>
struct HasTargetVolume : std::false_type {
  static void addData(VTKWriter<typename T::scalarType>& writer,const ParticleSetTpl<T>& pset) {}
};
template <typename T>
struct HasType : std::false_type {
  static void addData(VTKWriter<typename T::scalarType>& writer,const ParticleSetTpl<T>& pset) {}
};

template <typename A,typename B,bool has>
struct CopyPos {
  static FORCE_INLINE void copy(A& a,const B& b) {
    a._pos.x()=(typename A::scalarType)b._pos.x();
    a._pos.y()=(typename A::scalarType)b._pos.y();
    a._pos.z()=(typename A::scalarType)b._pos.z();
  }
};
template <typename A,typename B>
struct CopyPos<A,B,false> {
  static FORCE_INLINE void copy(A& a,const B& b) {}
};

template <typename A,typename B,bool has>
struct CopyVel {
  static FORCE_INLINE void copy(A& a,const B& b) {
    a._vel.x()=(typename A::scalarType)b._vel.x();
    a._vel.y()=(typename A::scalarType)b._vel.y();
    a._vel.z()=(typename A::scalarType)b._vel.z();
  }
};
template <typename A,typename B>
struct CopyVel<A,B,false> {
  static FORCE_INLINE void copy(A& a,const B& b) {}
};

template <typename A,typename B,bool has>
struct CopyNormal {
  static FORCE_INLINE void copy(A& a,const B& b) {
    a._normal.x()=(typename A::scalarType)b._normal.x();
    a._normal.y()=(typename A::scalarType)b._normal.y();
    a._normal.z()=(typename A::scalarType)b._normal.z();
  }
};
template <typename A,typename B>
struct CopyNormal<A,B,false> {
  static FORCE_INLINE void copy(A& a,const B& b) {}
};

//============================================================================
// Flip is one of the most efficient surface tracker
template <typename P_TYPE>
struct ParticleSetTpl : public Serializable {
  typedef typename P_TYPE::scalarType scalarType;
  typedef typename ScalarUtil<scalarType>::ScalarVec3 Vec3Type;
  typedef typename ScalarUtil<scalarType>::ScalarMat3 Mat3Type;
  typedef P_TYPE ParticleType;
  typedef P_TYPE value_type;
public:
  EIGEN_DEVICE_FUNC ParticleSetTpl():Serializable(typeid(ParticleSetTpl<P_TYPE>).name()) {}
  EIGEN_DEVICE_FUNC virtual ~ParticleSetTpl() {}
  virtual bool read(std::istream& is,IOData* dat) {
    readBinaryData(_pSet,is);
    return is.good();
  }
  virtual bool write(std::ostream& os,IOData* dat) const {
    writeBinaryData(_pSet,os);
    return os.good();
  }
  virtual std::shared_ptr<SerializableBase> copy() const {
    return std::shared_ptr<SerializableBase>(new ParticleSetTpl<P_TYPE>());
  }
  bool readFrame(std::istream& is) {
    sizeType nr;
    char comma;
    is >> nr;
    resize(nr);
    for(sizeType i=0; i<nr; i++) {
      Vec3Type& pos=get(i)._pos;
      is >> comma >> pos.x();
      ASSERT(comma == ',');
      is >> comma >> pos.y();
      ASSERT(comma == ',');
      is >> comma >> pos.z();
      ASSERT(comma == ',');
    }
    /*for(sizeType i=0;i<nr;i++)
    {
    	vec3Type& vel=get(i)._vel;
    	is >> comma >> vel.x();ASSERT(comma == ',');
    	is >> comma >> vel.y();ASSERT(comma == ',');
    	is >> comma >> vel.z();ASSERT(comma == ',');
    }*/
    return is.good();
  }
  bool writeFrame(std::ostream& os) const {
    os << size() << ",";
    for(sizeType i=0; i<(sizeType)_pSet.size(); i++) {
      const Vec3Type& pos=get(i)._pos;
      os << pos.x() << "," << pos.y() << "," << pos.z() << ",";
    }
    /*for(sizeType i=0;i<(sizeType)_pSet.size();i++)
    {
    	const vec3Type& vel=get(i)._vel;
    	os << vel.x() << "," << vel.y() << "," << vel.z();
    	if(i < (sizeType)_pSet.size()-1)
    		os << ",";
    }*/
    return os.good();
  }
  void writeVTK(const std::string& path) const {
    VTKWriter<scalarType> writer("Particles",path,true);
    writeVTK(writer);
  }
  void writeVTK(VTKWriter<scalarType>& writer) const;
  void writeVTK(const Mat3Type& R,const Vec3Type& C,const std::string& path) const {
    VTKWriter<scalarType> writer("Particles",path,true);
    writeVTK(R,C,writer);
  }
  void writeVTK(const Mat3Type& R,const Vec3Type& C,VTKWriter<scalarType>& writer) const;
  void writeNormalVTK(const std::string& path,scalar len) const {
    std::vector<Vec3Type,Eigen::aligned_allocator<Vec3Type> > pos;
    const sizeType n=size();
    for(sizeType i=0; i<n; i++) {
      pos.push_back(_pSet[i]._pos);
      pos.push_back(_pSet[i]._pos+_pSet[i]._normal*len);
    }

    VTKWriter<scalarType> writer("Particle Normals",path,true);
    writer.appendPoints(pos.begin(),pos.end());
    typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itBeg(0,2,0);
    typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itEnd(_pSet.size(),2,0);
    writer.appendCells(itBeg,itEnd,
                       VTKWriter<scalarType>::LINE);
  }
  void writeVelVTK(const std::string& path,const scalarType len) const {
    std::vector<Vec3Type,Eigen::aligned_allocator<Vec3Type> > pos;
    std::vector<scalarType> css;
    const sizeType n=size();
    for(sizeType i=0; i<n; i++) {
      pos.push_back(_pSet[i]._pos);
      pos.push_back(_pSet[i]._pos+_pSet[i]._vel*len);
      css.push_back(0.0f);
      css.push_back(1.0f);
    }

    VTKWriter<scalarType> writer("Particle Normals",path,true);
    writer.appendPoints(pos.begin(),pos.end());
    writer.appendCustomPointData("Color",css.begin(),css.end());
    typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itBeg(0,2,0);
    typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itEnd(_pSet.size(),2,0);
    writer.appendCells(itBeg,itEnd,VTKWriter<scalarType>::LINE);
  }
  sizeType size() const {
    return _pSet.size();
  }
  FORCE_INLINE P_TYPE& operator[](const sizeType& i) {
    return get(i);
  }
  FORCE_INLINE const P_TYPE& operator[](const sizeType& i) const {
    return get(i);
  }
  FORCE_INLINE P_TYPE& get(const sizeType& i) {
    return _pSet[i];
  }
  FORCE_INLINE const P_TYPE& get(const sizeType& i) const {
    return _pSet[i];
  }
  //P_TYPE* getPtr(){return &(_pSet[0]);}
  //const P_TYPE* getPtr() const{return &(_pSet[0]);}
  void clear() {
    _pSet.clear();
  }
  void resize(const sizeType& n) {
    _pSet.resize(n);
  }
  void addParticle(const P_TYPE& p) {
    _pSet.push_back(p);
  }
  void swap(ParticleSetTpl& other) {
    _pSet.swap(other._pSet);
  }
  void append(const ParticleSetTpl& other) {
    _pSet.insert(_pSet.end(),other.begin(),other.end());
  }
  typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::iterator begin() {
    return _pSet.begin();
  }
  typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::iterator end() {
    return _pSet.end();
  }
  typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::const_iterator begin() const {
    return _pSet.begin();
  }
  typename std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> >::const_iterator end() const {
    return _pSet.end();
  }
  template<typename P_TYPE_FROM>
  ParticleSetTpl& copy(const ParticleSetTpl<P_TYPE_FROM>& other) {
    resize(other.size());
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<other.size(); i++)
      get(i).template copy<P_TYPE_FROM>(other.get(i));
    return *this;
  }
  template<typename COMP>
  void sort(const COMP& cmp) {
    std::sort(_pSet.begin(),_pSet.end(),cmp);
  }
protected:
  std::vector<P_TYPE,Eigen::aligned_allocator<P_TYPE> > _pSet;
};
template <typename T>
struct ParticleBase : public Serializable {
  typedef T scalarType;
  typedef sizeType SIZE;
  typedef char CHAR;
  typedef typename ScalarUtil<T>::ScalarVec3 PT3;
  typedef typename ScalarUtil<T>::ScalarMat3 MAT3;
  EIGEN_DEVICE_FUNC ParticleBase():Serializable(typeid(ParticleBase<T>).name()),_pos(0.0f,0.0f,0.0f) {}
  template<typename OTHER>
  ParticleBase& copy(const OTHER& p) {
    CopyPos<ParticleBase,OTHER,HasPos<OTHER>::value>::copy(*this,p);
    return *this;
  }
  virtual bool read(std::istream& is) {
    return readBinaryData(_pos,is).good();
  }
  virtual bool write(std::ostream& os) const {
    return writeBinaryData(_pos,os).good();
  }
  virtual bool operator==(const ParticleBase& other) const {
    return other._pos == _pos;
  }
  PT3 _pos;
protected:
  ParticleBase(const std::string& name):Serializable(name),_pos(0.0f,0.0f,0.0f) {}
};
template <typename T>
struct Particle : public ParticleBase<T> {
  using typename ParticleBase<T>::PT3;
  EIGEN_DEVICE_FUNC Particle():ParticleBase<T>(typeid(Particle<T>).name()),_vel(0.0f,0.0f,0.0f) {}
  template<typename OTHER>
  Particle& copy(const OTHER& p) {
    ParticleBase<T>::copy(p);
    CopyVel<Particle,OTHER,HasVel<OTHER>::value>::copy(*this,p);
    return *this;
  }
  virtual bool read(std::istream& is) {
    return ParticleBase<T>::read(is) && readBinaryData(_vel,is).good();
  }
  virtual bool write(std::ostream& os) const {
    return ParticleBase<T>::write(os) && writeBinaryData(_vel,os).good();
  }
  virtual bool operator==(const Particle& other) const {
    return ParticleBase<T>::operator==(other) && other._vel == _vel;
  }
  PT3 _vel;
protected:
  Particle(const std::string& name):ParticleBase<T>(name),_vel(0.0f,0.0f,0.0f) {}
};
template <typename T>
struct ParticleN : public Particle<T> {
  using typename Particle<T>::PT3;
  EIGEN_DEVICE_FUNC ParticleN():Particle<T>(typeid(ParticleN<T>).name()),_normal(0.0f,0.0f,0.0f) {}
  template<typename OTHER>
  ParticleN& copy(const OTHER& p) {
    Particle<T>::copy(p);
    CopyNormal<ParticleN,OTHER,HasNormal<OTHER>::value>::copy(*this,p);
    return *this;
  }
  virtual bool read(std::istream& is) {
    return Particle<T>::read(is) && readBinaryData(_normal,is).good();
  }
  virtual bool write(std::ostream& os) const {
    return Particle<T>::write(os) && writeBinaryData(_normal,os).good();
  }
  virtual bool operator==(const ParticleN& other) const {
    return Particle<T>::operator==(other) && other._normal == _normal;
  }
  PT3 _normal;
protected:
  ParticleN(const std::string& name):Particle<T>(name),_normal(0.0f,0.0f,0.0f) {}
};
template <typename T>
struct ParticleSPH : public ParticleN<T> {
  using typename Particle<T>::PT3;
  using typename Particle<T>::MAT3;
  EIGEN_DEVICE_FUNC ParticleSPH():ParticleN<T>(typeid(ParticleSPH<T>).name()) {
    _oldPos.setZero();
    _oldVel.setZero();
    _density=0.0f;
    _pressure=0.0f;
    _weight=0.0f;
    _targetVolume=0.0f;
    _type=0;
  }
  template<typename OTHER>
  ParticleSPH& copy(const OTHER& p) {
    ParticleN<T>::copy(p);
    const ParticleSPH<T>* pSPH=static_cast<const ParticleSPH<T>*>(&p);
    if(pSPH) {
      _density=(T)(pSPH->_density);
      _pressure=(T)(pSPH->_pressure);
      _weight=(T)(pSPH->_weight);
      _targetVolume=(T)(pSPH->_targetVolume);
      _type=pSPH->_type;
    }
    return *this;
  }
  virtual bool read(std::istream& is) {
    return ParticleN<T>::read(is) &&
           readBinaryData(_density,is).good() &&
           readBinaryData(_pressure,is).good() &&
           readBinaryData(_weight,is).good() &&
           readBinaryData(_targetVolume,is).good() &&
           readBinaryData(_type,is).good();
  }
  virtual bool write(std::ostream& os) const {
    return ParticleN<T>::write(os) &&
           writeBinaryData(_density,os).good() &&
           writeBinaryData(_pressure,os).good() &&
           writeBinaryData(_weight,os).good() &&
           writeBinaryData(_targetVolume,os).good() &&
           writeBinaryData(_type,os).good();
  }
  virtual bool operator==(const ParticleSPH& other) const {
    return ParticleN<T>::operator==(other) &&
           other._density == _density &&
           other._pressure == _pressure &&
           other._weight == _weight &&
           other._targetVolume == _targetVolume &&
           other._type == _type;
  }
  PT3 _oldPos,_oldVel;
  T _density,_pressure;
  T _weight,_targetVolume;
  char _type;
protected:
  ParticleSPH(const std::string& name):ParticleN<T>(name) {
    _oldPos.setZero();
    _oldVel.setZero();
    _density=0.0f;
    _pressure=0.0f;
    _weight=0.0f;
    _targetVolume=0.0f;
    _type=0;
  }
};

//============================================================================
// For VTK Writer
template <typename PS_TYPE,typename VALUE_TYPE>
struct PSetIter {
  PSetIter(const PS_TYPE& ps,sizeType id):_ps(ps),_id(id) {}
  void operator++() {
    _id++;
  }
  bool operator!=(const PSetIter& other) const {
    return _id != other._id;
  }
  virtual const VALUE_TYPE& operator*() const=0;
  virtual VALUE_TYPE& operator*()=0;
  const PS_TYPE& _ps;
  sizeType _id;
};
#define PSET_ITER(TYPE,VAL,FIELD)	\
template <typename PS_TYPE>	\
struct PSetIter##VAL : public PSetIter<PS_TYPE,typename PS_TYPE::ParticleType::TYPE>	\
{	\
    typedef PSetIter<PS_TYPE,typename PS_TYPE::ParticleType::TYPE> Parent; \
    typedef typename PS_TYPE::ParticleType::TYPE value_type;	\
    using Parent::_ps; \
    using Parent::_id; \
    PSetIter##VAL(const PS_TYPE& ps,sizeType id):Parent(ps,id){}	\
    virtual const value_type& operator*() const{const typename PS_TYPE::ParticleType& p=_ps[_id];return FIELD;}	\
    virtual value_type& operator*(){typename PS_TYPE::ParticleType& p=const_cast<typename PS_TYPE::ParticleType&>(_ps[_id]);return FIELD;}	\
};
PSET_ITER(PT3,Pos,p._pos)
PSET_ITER(PT3,Vel,p._vel)
PSET_ITER(PT3,Normal,p._normal)
PSET_ITER(PT3,OldPos,p._oldPos)
PSET_ITER(PT3,OldVel,p._oldVel)
PSET_ITER(scalarType,Density,p._density)
PSET_ITER(scalarType,Pressure,p._pressure)
PSET_ITER(scalarType,Weight,p._weight)
PSET_ITER(scalarType,TargetVolume,p._targetVolume)
PSET_ITER(CHAR,Type,p._type)
PSET_ITER(PT3,PosAddVel,p._pos+p._vel)
PSET_ITER(scalarType,NormalX,p._normal[0])
PSET_ITER(scalarType,NormalY,p._normal[1])
PSET_ITER(scalarType,NormalZ,p._normal[2])
template <typename PS_TYPE>
struct PSetIterTransPos : public PSetIter<PS_TYPE,typename PS_TYPE::ParticleType::PT3> {
  typedef PSetIter<PS_TYPE,typename PS_TYPE::ParticleType::PT3> Parent;
  typedef typename PS_TYPE::Vec3Type value_type;
  typedef typename PS_TYPE::Mat3Type MAT3;
  using Parent::_ps;
  using Parent::_id;
  PSetIterTransPos(const PS_TYPE& ps,sizeType id):Parent(ps,id) {}
  virtual const value_type& operator*() const {
    const typename PS_TYPE::ParticleType& p=_ps[_id];
    const_cast<value_type&>(_tmp)=_R*p._pos+_C;
    return _tmp;
  }
  virtual value_type& operator*() {
    typename PS_TYPE::ParticleType& p=const_cast<typename PS_TYPE::ParticleType&>(_ps[_id]);
    _tmp=_R*p._pos+_C;
    return _tmp;
  }
  value_type _C,_tmp;
  MAT3 _R;
};

template <typename P_TYPE>
void ParticleSetTpl<P_TYPE>::writeVTK(VTKWriter<scalarType>& writer) const
{
  writer.setRelativeIndex();
  writer.appendPoints(PSetIterPos<ParticleSetTpl<P_TYPE> >(*this,0),PSetIterPos<ParticleSetTpl<P_TYPE> >(*this,size()));
  typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itBeg(0,1,0);
  typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itEnd(_pSet.size(),1,0);
  writer.appendCells(itBeg,itEnd,VTKWriter<scalarType>::POINT,true);
  HasVel<P_TYPE>::addData(writer,*this);
  HasDensity<P_TYPE>::addData(writer,*this);
  HasPressure<P_TYPE>::addData(writer,*this);
  HasWeight<P_TYPE>::addData(writer,*this);
  HasTargetVolume<P_TYPE>::addData(writer,*this);
  HasType<P_TYPE>::addData(writer,*this);
}

template <typename P_TYPE>
void ParticleSetTpl<P_TYPE>::writeVTK(const Mat3Type& R,const Vec3Type& C,VTKWriter<scalarType>& writer) const
{
  writer.setRelativeIndex();
  PSetIterTransPos<ParticleSetTpl<P_TYPE> > beg(*this,0);
  PSetIterTransPos<ParticleSetTpl<P_TYPE> > end(*this,size());
  beg._R=end._R=R;
  beg._C=end._C=C;
  writer.appendPoints(beg,end);
  typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itBeg(0,1,0);
  typename VTKWriter<scalarType>::template IteratorIndex<Vec3i> itEnd(_pSet.size(),1,0);
  writer.appendCells(itBeg,itEnd,VTKWriter<scalarType>::POINT,true);
  HasVel<P_TYPE>::addData(writer,*this);
  HasDensity<P_TYPE>::addData(writer,*this);
  HasPressure<P_TYPE>::addData(writer,*this);
  HasWeight<P_TYPE>::addData(writer,*this);
  HasTargetVolume<P_TYPE>::addData(writer,*this);
  HasType<P_TYPE>::addData(writer,*this);
}

typedef ParticleSetTpl<ParticleBase<scalarF> > ParticleSetBF;
typedef ParticleSetTpl<Particle<scalarF> > ParticleSetF;
typedef ParticleSetTpl<ParticleN<scalarF> > ParticleSetNF;
typedef ParticleSetTpl<ParticleSPH<scalarF> > ParticleSPHF;

typedef ParticleSetTpl<ParticleBase<scalarD> > ParticleSetBD;
typedef ParticleSetTpl<Particle<scalarD> > ParticleSetD;
typedef ParticleSetTpl<ParticleN<scalarD> > ParticleSetND;
typedef ParticleSetTpl<ParticleSPH<scalarD> > ParticleSetSPHD;

typedef ParticleSetTpl<ParticleBase<scalar> > ParticleSetB;
typedef ParticleSetTpl<Particle<scalar> > ParticleSet;
typedef ParticleSetTpl<ParticleN<scalar> > ParticleSetN;
typedef ParticleSetTpl<ParticleSPH<scalar> > ParticleSetSPH;

template <typename T>
struct HasPos<ParticleBase<T> > : public std::true_type {};

template <typename T>
struct HasPos<Particle<T> > : public std::true_type {};
template <typename T>
struct HasVel<Particle<T> > : public std::true_type {
  typedef ParticleSetTpl<Particle<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointVectorData("velocity",PSetIterVel<PS_TYPE>(pset,0),PSetIterVel<PS_TYPE>(pset,pset.size()));
  }
};

template <typename T>
struct HasPos<ParticleN<T> > : public std::true_type {};
template <typename T>
struct HasVel<ParticleN<T> > : public std::true_type {
  typedef ParticleSetTpl<ParticleN<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointVectorData("velocity",PSetIterVel<PS_TYPE>(pset,0),PSetIterVel<PS_TYPE>(pset,pset.size()));
  }
};
template <typename T>
struct HasNormal<ParticleN<T> > : public std::true_type {};

template <typename T>
struct HasVel<ParticleSPH<T> > : public std::true_type {
  typedef ParticleSetTpl<ParticleSPH<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointVectorData("velocity",PSetIterVel<PS_TYPE>(pset,0),PSetIterVel<PS_TYPE>(pset,pset.size()));
  }
};
template <typename T>
struct HasDensity<ParticleSPH<T> > : std::true_type {
  typedef ParticleSetTpl<ParticleSPH<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointData("density",PSetIterDensity<PS_TYPE>(pset,0),PSetIterDensity<PS_TYPE>(pset,pset.size()));
  }
};
template <typename T>
struct HasPressure<ParticleSPH<T> > : std::false_type {
  typedef ParticleSetTpl<ParticleSPH<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointData("pressure",PSetIterPressure<PS_TYPE>(pset,0),PSetIterPressure<PS_TYPE>(pset,pset.size()));
  }
};
template <typename T>
struct HasWeight<ParticleSPH<T> > : std::true_type {
  typedef ParticleSetTpl<ParticleSPH<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointData("weight",PSetIterWeight<PS_TYPE>(pset,0),PSetIterWeight<PS_TYPE>(pset,pset.size()));
  }
};
template <typename T>
struct HasTargetVolume<ParticleSPH<T> > : std::false_type {
  typedef ParticleSetTpl<ParticleSPH<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointData("targetVolume",PSetIterTargetVolume<PS_TYPE>(pset,0),PSetIterTargetVolume<PS_TYPE>(pset,pset.size()));
  }
};
template <typename T>
struct HasType<ParticleSPH<T> > : std::false_type {
  typedef ParticleSetTpl<ParticleSPH<T> > PS_TYPE;
  static void addData(VTKWriter<T>& writer,const PS_TYPE& pset) {
    writer.appendCustomPointData("type",PSetIterType<PS_TYPE>(pset,0),PSetIterType<PS_TYPE>(pset,pset.size()));
  }
};

PRJ_END

#endif
