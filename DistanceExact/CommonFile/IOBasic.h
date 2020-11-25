#ifndef IO_BASIC_H
#define IO_BASIC_H

#include "IOFwd.h"
#include <iostream>
#include "MathBasic.h"
#include <Eigen/Sparse>
#include <unordered_map>
#include "CollisionDetection.h"

PRJ_BEGIN

//minimal serializable system
std::shared_ptr<IOData> getIOData();
std::ostream& writeSerializableData(std::shared_ptr<SerializableBase> val,std::ostream& os,IOData* dat=NULL);
std::istream& readSerializableData(std::shared_ptr<SerializableBase>& val,std::istream& is,IOData* dat=NULL);
void registerType(IOData* dat,std::shared_ptr<SerializableBase> type);
template <typename T>FORCE_INLINE void registerType(IOData* dat)
{
  registerType(dat,std::shared_ptr<SerializableBase>(new T));
}

//io for shared_ptr
template <typename T>FORCE_INLINE std::ostream& writeBinaryData(std::shared_ptr<T> val,std::ostream& os,IOData* dat=NULL)
{
  std::shared_ptr<SerializableBase> valS=std::dynamic_pointer_cast<SerializableBase>(val);
  return writeSerializableData(valS,os,dat);
}
template <typename T>FORCE_INLINE std::istream& readBinaryData(std::shared_ptr<T>& val,std::istream& is,IOData* dat=NULL)
{
  std::shared_ptr<SerializableBase> valS;
  readSerializableData(valS,is,dat);
  val=std::dynamic_pointer_cast<T>(valS);
  return is;
}
template <typename T>FORCE_INLINE std::ostream& writeBinaryData(std::weak_ptr<T> val,std::ostream& os,IOData* dat=NULL)
{
  std::shared_ptr<T> val_shared=val.lock();
  std::shared_ptr<SerializableBase> valS=std::dynamic_pointer_cast<SerializableBase>(val_shared);
  return writeSerializableData(valS,os,dat);
}
template <typename T>FORCE_INLINE std::istream& readBinaryData(std::weak_ptr<T>& val,std::istream& is,IOData* dat=NULL)
{
  std::shared_ptr<SerializableBase> valS;
  readSerializableData(valS,is,dat);
  val=std::dynamic_pointer_cast<T>(valS);
  return is;
}
FORCE_INLINE std::ostream& writeBinaryData(const SerializableBase& val,std::ostream& os,IOData* dat=NULL)
{
  val.write(os,dat);
  return os;
}
FORCE_INLINE std::istream& readBinaryData(SerializableBase& val,std::istream& is,IOData* dat=NULL)
{
  val.read(is,dat);
  return is;
}

//io for basic type
struct IOData;
typedef Eigen::Triplet<scalarF,sizeType> TripF;
typedef Eigen::Triplet<scalarD,sizeType> TripD;
#define IO_BASIC_DECL(T)	\
std::ostream& writeBinaryData(const T& val,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(T& val,std::istream& is,IOData* dat=NULL);
IO_BASIC_DECL(char)
IO_BASIC_DECL(unsigned char)
IO_BASIC_DECL(short)
IO_BASIC_DECL(unsigned short)
IO_BASIC_DECL(int)
IO_BASIC_DECL(unsigned int)
//IO_BASIC_DECL(scalarD)
IO_BASIC_DECL(bool)
#ifdef USE_QUAD_SIZE
IO_BASIC_DECL(sizeType)
#endif
IO_BASIC_DECL(TripF)
IO_BASIC_DECL(TripD)
#undef IO_BASIC_DECL

//io for string
std::istream& readBinaryData(std::string& str,std::istream& is,IOData* dat=NULL);
std::ostream& writeBinaryData(const std::string& str,std::ostream& os,IOData* dat=NULL);

//io for float is double
std::ostream& writeBinaryData(scalarF val,std::ostream& os,IOData* dat=NULL);
std::istream& readBinaryData(scalarF& val,std::istream& is,IOData* dat=NULL);
std::ostream& writeBinaryData(scalarD val,std::ostream& os,IOData* dat=NULL);
std::istream& readBinaryData(scalarD& val,std::istream& is,IOData* dat=NULL);

//io for fixed matrix
#define NO_CONFLICT
#define FIXED_ONLY
#include "BeginAllEigen.h"
//redefine atomic operation
#undef NAME_EIGEN
#define NAME_EIGEN(type,NAME,size1,size2) \
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat=NULL);
//realize
NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE()
NAME_EIGEN_MAT_ALLTYPES_SPECIALSIZE()
#include "EndAllEigen.h"

//io for non-fixed matrix
#define NO_CONFLICT
#define NON_FIXED_ONLY
#include "BeginAllEigen.h"
//redefine atomic operation
#undef NAME_EIGEN
#define NAME_EIGEN(type,NAME,size1,size2) \
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat=NULL);
//realize
NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE()
#include "EndAllEigen.h"

//io for quaternion
#define IO_FIXED_QUAT_DECL(NAMEQ,NAMET,NAMEA)	\
std::ostream& writeBinaryData(const NAMEQ& v,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(NAMEQ& v,std::istream& is,IOData* dat=NULL);         \
std::ostream& writeBinaryData(const NAMET& v,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(NAMET& v,std::istream& is,IOData* dat=NULL);         \
std::ostream& writeBinaryData(const NAMEA& v,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(NAMEA& v,std::istream& is,IOData* dat=NULL);
IO_FIXED_QUAT_DECL(Quatd,Transd,Affined)
IO_FIXED_QUAT_DECL(Quatf,Transf,Affinef)
#undef IO_FIXED_QUAT_DECL

//io for shape
#define IO_SHAPE_DECL(TYPE)	\
std::ostream& writeBinaryData(const TYPE& b,std::ostream& os,IOData* dat=NULL);	\
std::istream& readBinaryData(TYPE& b,std::istream& is,IOData* dat=NULL);
#define IO_SHAPE_ALL_DECL(NAME) \
IO_SHAPE_DECL(BBox<NAME>) \
IO_SHAPE_DECL(LineSegTpl<NAME>) \
IO_SHAPE_DECL(PlaneTpl<NAME>) \
IO_SHAPE_DECL(TriangleTpl<NAME>) \
IO_SHAPE_DECL(TetrahedronTpl<NAME>) \
IO_SHAPE_DECL(OBBTpl<NAME>) \
IO_SHAPE_DECL(KDOP18<NAME>) \
IO_SHAPE_DECL(Sphere<NAME>)
IO_SHAPE_ALL_DECL(scalarD)
IO_SHAPE_ALL_DECL(scalarF)
IO_SHAPE_ALL_DECL(sizeType)
IO_SHAPE_ALL_DECL(char)
IO_SHAPE_ALL_DECL(unsigned char)
#undef IO_SHAPE_ALL_DECL
#undef IO_SHAPE_DECL

//io for unordered_map
template <typename K,typename V,typename H>
FORCE_INLINE std::ostream& writeBinaryData(const std::unordered_map<K,V,H>& val,std::ostream& os,IOData* dat=NULL)
{
  sizeType nr=(sizeType)val.size();
  writeBinaryData(nr,os);
  for(typename std::unordered_map<K,V,H>::const_iterator beg=val.begin(),end=val.end(); beg!=end; beg++) {
    writeBinaryData(beg->first,os,dat);
    writeBinaryData(beg->second,os,dat);
  }
  return os;
}
template <typename K,typename V,typename H>
FORCE_INLINE std::istream& readBinaryData(std::unordered_map<K,V,H>& val,std::istream& is,IOData* dat=NULL)
{
  K k;
  V v;
  sizeType nr;
  val.clear();
  readBinaryData(nr,is);
  for(sizeType i=0; i<nr; i++) {
    readBinaryData(k,is,dat);
    readBinaryData(v,is,dat);
    val[k]=v;
  }
  return is;
}

PRJ_END

#endif
