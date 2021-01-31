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
void registerTypeAs(IOData* dat,std::shared_ptr<SerializableBase> alias,std::shared_ptr<SerializableBase> type);
template <typename T>FORCE_INLINE void registerType(IOData* dat)
{
  registerType(dat,std::shared_ptr<SerializableBase>(new T));
}
template <typename A,typename T>FORCE_INLINE void registerTypeAs(IOData* dat)
{
  registerTypeAs(dat,std::shared_ptr<SerializableBase>(new A),std::shared_ptr<SerializableBase>(new T));
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

//io for matrix
template <typename type,int size1,int size2>
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat=NULL)
{
  sizeType d0=v.rows();
  sizeType d1=v.cols();
  os.write((char*)&d0,sizeof(sizeType));
  os.write((char*)&d1,sizeof(sizeType));
  for(sizeType r=0; r<d0; r++)
    for(sizeType c=0; c<d1; c++)
      writeBinaryData(v(r,c),os);
  return os;
}
template <typename type,int size1,int size2>
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat=NULL)
{
  sizeType d0;
  sizeType d1;
  is.read((char*)&d0,sizeof(sizeType));
  is.read((char*)&d1,sizeof(sizeType));
  v.resize(d0,d1);
  for(sizeType r=0; r<d0; r++)
    for(sizeType c=0; c<d1; c++)
      readBinaryData(v(r,c),is);
  return is;
}

//io for quaternion
template <typename type>
std::ostream& writeBinaryData(const Eigen::Quaternion<type>& v,std::ostream& os,IOData* dat=NULL)
{
  writeBinaryData(v.w(),os);
  writeBinaryData(v.x(),os);
  writeBinaryData(v.y(),os);
  writeBinaryData(v.z(),os);
  return os;
}
template <typename type>
std::istream& readBinaryData(Eigen::Quaternion<type>& v,std::istream& is,IOData* dat=NULL)
{
  readBinaryData(v.w(),is);
  readBinaryData(v.x(),is);
  readBinaryData(v.y(),is);
  readBinaryData(v.z(),is);
  return is;
}
template <typename type>
std::ostream& writeBinaryData(const Eigen::Translation<type,3>& v,std::ostream& os,IOData* dat=NULL)
{
  writeBinaryData(v.x(),os);
  writeBinaryData(v.y(),os);
  writeBinaryData(v.z(),os);
  return os;
}
template <typename type>
std::istream& readBinaryData(Eigen::Translation<type,3>& v,std::istream& is,IOData* dat=NULL)
{
  readBinaryData(v.x(),is);
  readBinaryData(v.y(),is);
  readBinaryData(v.z(),is);
  return is;
}
template <typename type>
std::ostream& writeBinaryData(const Eigen::Transform<type,3,Eigen::Affine>& v,std::ostream& os,IOData* dat=NULL)
{
  Eigen::Matrix<type,3,3> r=v.translation();
  Eigen::Matrix<type,3,1> t=v.linear();
  writeBinaryData(r,os);
  writeBinaryData(t,os);
  return os;
}
template <typename type>
std::istream& readBinaryData(Eigen::Transform<type,3,Eigen::Affine>& v,std::istream& is,IOData* dat=NULL)
{
  Eigen::Matrix<type,3,3> r;
  Eigen::Matrix<type,3,1> t;
  readBinaryData(r,is);
  readBinaryData(t,is);
  v.translation()=r;
  v.linear()=t;
  return is;
}

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
