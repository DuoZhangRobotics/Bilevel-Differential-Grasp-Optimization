#include "IOBasic.h"
#include "CollisionDetection.h"
#include <unordered_map>
#include <fstream>

PRJ_BEGIN

//SerializableBase
bool SerializableBase::read(std::istream& is) {
  FUNCTION_NOT_IMPLEMENTED
  return false;
}
bool SerializableBase::write(std::ostream& os) const {
  FUNCTION_NOT_IMPLEMENTED
  return false;
}
bool SerializableBase::read(std::istream& is,IOData* dat) {
  return read(is);
}
bool SerializableBase::write(std::ostream& os,IOData* dat) const {
  return write(os);
}
bool SerializableBase::read(const std::string& path) {
  std::ifstream is(path,std::ios::binary);
  std::shared_ptr<IOData> dat=getIOData();
  return read(is,dat.get());
}
bool SerializableBase::write(const std::string& path) const {
  std::ofstream os(path,std::ios::binary);
  std::shared_ptr<IOData> dat=getIOData();
  return write(os,dat.get());
}
std::shared_ptr<SerializableBase> SerializableBase::copy() const {
  FUNCTION_NOT_IMPLEMENTED
  return std::shared_ptr<SerializableBase>();
}

//Serializable
std::string Serializable::type() const {
  return _type;
}
void Serializable::setType(const std::string& type) {
  _type=type;
}
void Serializable::setSerializableType(const std::string& type) {
  _type=type;
}

//IOData
struct IOData {
public:
  typedef std::shared_ptr<SerializableBase> SPTR;
  typedef std::allocator<std::pair<const SPTR,sizeType> > ALLOC_MAP_WR;
  typedef std::allocator<std::pair<const sizeType,SPTR> > ALLOC_MAP_RD;
  typedef std::allocator<std::pair<const std::string,SPTR> > ALLOC_TYPE_SET;
  typedef std::unordered_map<SPTR,sizeType,std::hash<SPTR>,std::equal_to<SPTR>,ALLOC_MAP_WR> MAP_WR;
  typedef std::unordered_map<sizeType,SPTR,std::hash<sizeType>,std::equal_to<sizeType>,ALLOC_MAP_RD> MAP_RD;
  typedef std::unordered_map<std::string,SPTR,std::hash<std::string>,std::equal_to<std::string>,ALLOC_TYPE_SET> TYPE_SET;
public:
  IOData():_index(0) {}
  sizeType getIndex() {
    return _index++;
  }
  void registerType(std::shared_ptr<SerializableBase> type) {
    ASSERT_MSG(!type->type().empty(),"Given type doesn't support shared_ptr serialization!");
    if(_tSet.find(type->type()) != _tSet.end())
      ASSERT_MSGV(typeid(*_tSet[type->type()]) == typeid(*type),"Conflicit type id: %s",type->type().c_str())
      _tSet[type->type()]=type;
  }
  template <typename T>void createNew(std::istream& is,std::shared_ptr<T>& val) const {
    std::string type;
    readBinaryData(type,is);
    ASSERT_MSG(type != "","Type not found!")
    for(TYPE_SET::const_iterator beg=_tSet.begin(),end=_tSet.end(); beg!=end; beg++)
      if(beg->first == type) {
        val=std::dynamic_pointer_cast<T>(beg->second->copy());
        return;
      }
    ASSERT_MSG(false,"Cannot find compatible type!")
  }
  MAP_WR _ptrMapWr;
  MAP_RD _ptrMapRd;
private:
  sizeType _index;
  TYPE_SET _tSet;
};
//minimal serializable system
std::shared_ptr<IOData> getIOData()
{
  return std::shared_ptr<IOData>(new IOData);
}
std::ostream& writeSerializableData(std::shared_ptr<SerializableBase> val,std::ostream& os,IOData* dat)
{
  ASSERT_MSG(dat,"You must provide pool for serialize shared_ptr!")
  if(val.get() == NULL) {
    writeBinaryData((sizeType)-1,os);
    return os;
  }
  std::shared_ptr<SerializableBase> ptrS=std::dynamic_pointer_cast<SerializableBase>(val);
  ASSERT_MSG(ptrS,"Not serializable type!")
  IOData::MAP_WR::const_iterator it=dat->_ptrMapWr.find(ptrS);
  if(it == dat->_ptrMapWr.end()) {
    sizeType id=dat->getIndex();
    writeBinaryData(id,os,dat);
    writeBinaryData(val->type(),os,dat);
    dat->_ptrMapWr[val]=id;
    writeBinaryData(*val,os,dat);
  } else {
    writeBinaryData(it->second,os,dat);
  }
  return os;
}
std::istream& readSerializableData(std::shared_ptr<SerializableBase>& val,std::istream& is,IOData* dat)
{
  ASSERT_MSG(dat,"You must provide pool for serialize shared_ptr!")
  sizeType id;
  readBinaryData(id,is,dat);
  if(id == (sizeType)-1) {
    val.reset((Serializable*)NULL);
    return is;
  }
  IOData::MAP_RD::const_iterator it=dat->_ptrMapRd.find(id);
  if(it == dat->_ptrMapRd.end()) {
    dat->createNew(is,val);
    dat->_ptrMapRd[id]=val;
    readBinaryData(*val,is,dat);
  } else {
    val=std::dynamic_pointer_cast<SerializableBase>(dat->_ptrMapRd[id]);
  }
  return is;
}
void registerType(IOData* dat,std::shared_ptr<SerializableBase> type)
{
  dat->registerType(type);
}

//io for basic type
#define IO_BASIC(T)                                                                                             \
std::ostream& writeBinaryData(const T& val,std::ostream& os,IOData* dat){os.write((char*)&val,sizeof(T));return os;}	\
std::istream& readBinaryData(T& val,std::istream& is,IOData* dat){is.read((char*)&val,sizeof(T));return is;}
IO_BASIC(char)
IO_BASIC(unsigned char)
IO_BASIC(short)
IO_BASIC(unsigned short)
IO_BASIC(int)
IO_BASIC(unsigned int)
//IO_BASIC(scalarD)
IO_BASIC(bool)
#ifdef USE_QUAD_SIZE
IO_BASIC(sizeType)
#endif
IO_BASIC(TripF)
IO_BASIC(TripD)
#undef IO_BASIC

//io for std::string
std::istream& readBinaryData(std::string& str,std::istream& is,IOData* dat)
{
  sizeType len;
  readBinaryData(len,is,dat);
  str.assign(len,' ');
  is.read(&(str[0]),len);
  return is;
}
std::ostream& writeBinaryData(const std::string& str,std::ostream& os,IOData* dat)
{
  const sizeType len=(sizeType)str.length();
  writeBinaryData(len,os,dat);
  return os.write(str.c_str(),len);
}

//type invariant IO
#if defined(FOUND_QUADMATH) && defined(__GNUC__)
#include "quadmath.h"
typedef __float128 FLOAT_TO_TYPE;
#else
typedef double FLOAT_TO_TYPE;
#endif
template <typename T>
struct IOType {
  typedef T TO_TYPE;
};
template <>
struct IOType<scalarF> {
  typedef FLOAT_TO_TYPE TO_TYPE;
};
template <>
struct IOType<scalarD> {
  typedef FLOAT_TO_TYPE TO_TYPE;
};

//io for float is double
std::ostream& writeBinaryData(scalarF val,std::ostream& os,IOData* dat)
{
  IOType<scalarF>::TO_TYPE valD=val;
  os.write((char*)&valD,sizeof(IOType<scalarF>::TO_TYPE));
  return os;
}
std::istream& readBinaryData(scalarF& val,std::istream& is,IOData* dat)
{
  IOType<scalarF>::TO_TYPE valD;
  is.read((char*)&valD,sizeof(IOType<scalarF>::TO_TYPE));
  val=(scalarF)valD;
  return is;
}
std::ostream& writeBinaryData(scalarD val,std::ostream& os,IOData* dat)
{
  IOType<scalarD>::TO_TYPE valD=std::convert<IOType<scalarD>::TO_TYPE>()(val);
  os.write((char*)&valD,sizeof(IOType<scalarD>::TO_TYPE));
  return os;
}
std::istream& readBinaryData(scalarD& val,std::istream& is,IOData* dat)
{
  IOType<scalarD>::TO_TYPE valD;
  is.read((char*)&valD,sizeof(IOType<scalarD>::TO_TYPE));
  val=scalarD(valD);
  return is;
}

//io for fixed matrix
#define NO_CONFLICT
#define FIXED_ONLY
#include "BeginAllEigen.h"
//redefine atomic operation
#undef NAME_EIGEN
#define NAME_EIGEN(type,NAME,size1,size2) \
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat) \
{ \
  typedef Eigen::Matrix<type,size1,size2> TYPE; \
  typedef IOType<type>::TO_TYPE TO_TYPE;  \
  sizeType d0=TYPE::RowsAtCompileTime;os.write((char*)&d0,sizeof(sizeType)); \
  sizeType d1=TYPE::ColsAtCompileTime;os.write((char*)&d1,sizeof(sizeType)); \
  for(sizeType r=0;r<TYPE::RowsAtCompileTime;r++) \
  for(sizeType c=0;c<TYPE::ColsAtCompileTime;c++) \
  { \
    TO_TYPE val=std::convert<TO_TYPE>()(v(r,c)); \
    os.write((char*)&val,sizeof(TO_TYPE)); \
  } \
  return os; \
} \
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat) \
{ \
  typedef Eigen::Matrix<type,size1,size2> TYPE; \
  typedef IOType<type>::TO_TYPE TO_TYPE;  \
  sizeType d0;is.read((char*)&d0,sizeof(sizeType));ASSERT(d0 == TYPE::RowsAtCompileTime) \
  sizeType d1;is.read((char*)&d1,sizeof(sizeType));ASSERT(d1 == TYPE::ColsAtCompileTime) \
  for(sizeType r=0;r<TYPE::RowsAtCompileTime;r++) \
  for(sizeType c=0;c<TYPE::ColsAtCompileTime;c++) \
  { \
    TO_TYPE val; \
    is.read((char*)&val,sizeof(TO_TYPE)); \
    v(r,c)=TYPE::Scalar(val); \
  } \
  return is; \
}
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
std::ostream& writeBinaryData(const Eigen::Matrix<type,size1,size2>& v,std::ostream& os,IOData* dat) \
{ \
  typedef IOType<type>::TO_TYPE TO_TYPE;  \
  sizeType d0=v.rows(); \
  sizeType d1=v.cols(); \
  os.write((char*)&d0,sizeof(sizeType)); \
  os.write((char*)&d1,sizeof(sizeType)); \
  if(d0 <= 1 || d1 <= 1) { \
    if(d0*d1 > 0)os.write((char*)(v.cast<TO_TYPE>().eval().data()),sizeof(TO_TYPE)*d0*d1); \
    return os; \
  } \
  Eigen::Matrix<TO_TYPE,-1,1> sub; \
  sub.setZero(d1); \
  for(sizeType r=0; r<d0; r++) { \
    sub=v.row(r).cast<TO_TYPE>(); \
    if(d1 > 0)os.write((char*)(sub.data()),sizeof(TO_TYPE)*d1); \
  } \
  return os; \
} \
std::istream& readBinaryData(Eigen::Matrix<type,size1,size2>& v,std::istream& is,IOData* dat) \
{ \
  typedef Eigen::Matrix<type,size1,size2> TYPE; \
  typedef IOType<type>::TO_TYPE TO_TYPE;  \
  typedef Eigen::Matrix<TO_TYPE,TYPE::RowsAtCompileTime,TYPE::ColsAtCompileTime> TONAME; \
  sizeType d0; \
  is.read((char*)&d0,sizeof(sizeType)); \
  sizeType d1; \
  is.read((char*)&d1,sizeof(sizeType)); \
  v.resize(d0,d1); \
  if(d0 <= 1 || d1 <= 1) { \
    TONAME vcast(d0,d1); \
    if(d0*d1 > 0)is.read((char*)(vcast.data()),sizeof(TO_TYPE)*d0*d1); \
    v=vcast.cast<TYPE::Scalar>(); \
    return is; \
  } \
  Eigen::Matrix<TO_TYPE,-1,1> sub; \
  sub.setZero(d1); \
  for(sizeType r=0; r<d0; r++) { \
    if(d1 > 0)is.read((char*)(sub.data()),sizeof(TO_TYPE)*d1); \
    v.row(r)=sub.cast<TYPE::Scalar>(); \
  } \
  return is; \
}
//realize
NAME_EIGEN_ROWCOL_ALLTYPES_SPECIALSIZE()
#include "EndAllEigen.h"

//io for quaternion
#define IO_FIXED_QUAT(NAMEQ,NAMET,NAMEA,TO_TYPE)			\
std::ostream& writeBinaryData(const NAMEQ& v,std::ostream& os,IOData* dat)	\
{									\
    TO_TYPE val=std::convert<TO_TYPE>()(v.w());			\
	os.write((char*)&val,sizeof(TO_TYPE));				\
    val=std::convert<TO_TYPE>()(v.x());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
    val=std::convert<TO_TYPE>()(v.y());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
    val=std::convert<TO_TYPE>()(v.z());					\
	os.write((char*)&val,sizeof(TO_TYPE));				\
    return os;                                          \
}                                                       \
std::istream& readBinaryData(NAMEQ& v,std::istream& is,IOData* dat)		\
{									\
    TO_TYPE val;                                        \
	is.read((char*)&val,sizeof(TO_TYPE));				\
    v.w()=NAMEQ::Scalar(val);                           \
	is.read((char*)&val,sizeof(TO_TYPE));				\
    v.x()=NAMEQ::Scalar(val);                           \
	is.read((char*)&val,sizeof(TO_TYPE));				\
    v.y()=NAMEQ::Scalar(val);                           \
	is.read((char*)&val,sizeof(TO_TYPE));				\
    v.z()=NAMEQ::Scalar(val);                           \
    return is;                                          \
}                                                       \
std::ostream& writeBinaryData(const NAMET& v,std::ostream& os,IOData* dat)	\
{									\
    TO_TYPE val=std::convert<TO_TYPE>()(v.x());         \
	os.write((char*)&val,sizeof(TO_TYPE));				\
    val=std::convert<TO_TYPE>()(v.y());                 \
	os.write((char*)&val,sizeof(TO_TYPE));				\
    val=std::convert<TO_TYPE>()(v.z());                 \
	os.write((char*)&val,sizeof(TO_TYPE));				\
    return os;                                          \
}                                                       \
std::istream& readBinaryData(NAMET& v,std::istream& is,IOData* dat)		\
{									\
    TO_TYPE val;                                        \
	is.read((char*)&val,sizeof(TO_TYPE));				\
    v.x()=NAMEQ::Scalar(val);                           \
	is.read((char*)&val,sizeof(TO_TYPE));				\
    v.y()=NAMEQ::Scalar(val);                           \
	is.read((char*)&val,sizeof(TO_TYPE));				\
    v.z()=NAMEQ::Scalar(val);                           \
    return is;                                          \
}                                                       \
std::ostream& writeBinaryData(const NAMEA& v,std::ostream& os,IOData* dat)	\
{                                   \
  Eigen::Matrix<NAMEA::Scalar,3,3> l=v.linear();                        \
  Eigen::Matrix<NAMEA::Scalar,3,1> t=v.translation();                   \
  writeBinaryData(l,os,dat);                                            \
  writeBinaryData(t,os,dat);                                            \
  return os;                                                            \
}                                                                       \
std::istream& readBinaryData(NAMEA& v,std::istream& is,IOData* dat)		\
{									\
  Eigen::Matrix<NAMEA::Scalar,3,3> l;                                   \
  Eigen::Matrix<NAMEA::Scalar,3,1> t;                                   \
  readBinaryData(l,is,dat);                                             \
  readBinaryData(t,is,dat);                                             \
  v.linear()=l;                                                         \
  v.translation()=t;                                                    \
  return is;                                                            \
}
IO_FIXED_QUAT(Quatd,Transd,Affined,IOType<scalarD>::TO_TYPE)
IO_FIXED_QUAT(Quatf,Transf,Affinef,IOType<scalarF>::TO_TYPE)
#undef IO_FIXED_QUAT

//io for shape
#define IO_SHAPE(TYPE) \
std::ostream& writeBinaryData(const TYPE& b,std::ostream& os,IOData* dat)   \
{   \
    b.write(os);    \
    return os;  \
}   \
std::istream& readBinaryData(TYPE& b,std::istream& is,IOData* dat)  \
{   \
    b.read(is);    \
    return is;  \
}
#define IO_SHAPE_ALL(NAME) \
IO_SHAPE(BBox<NAME>) \
IO_SHAPE(LineSegTpl<NAME>) \
IO_SHAPE(PlaneTpl<NAME>) \
IO_SHAPE(TriangleTpl<NAME>) \
IO_SHAPE(TetrahedronTpl<NAME>) \
IO_SHAPE(OBBTpl<NAME>) \
IO_SHAPE(KDOP18<NAME>) \
IO_SHAPE(Sphere<NAME>)
IO_SHAPE_ALL(scalarD)
IO_SHAPE_ALL(scalarF)
IO_SHAPE_ALL(sizeType)
IO_SHAPE_ALL(char)
IO_SHAPE_ALL(unsigned char)
#undef IO_SHAPE_ALL
#undef IO_SHAPE

PRJ_END
