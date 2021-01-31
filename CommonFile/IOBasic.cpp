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
    if(_tSet.find(type->type()) != _tSet.end()) {
      ASSERT_MSGV(typeid(*_tSet[type->type()]) == typeid(*type),"Conflicit type id: %s",type->type().c_str())
    }
    _tSet[type->type()]=type;
  }
  void registerTypeAs(std::shared_ptr<SerializableBase> alias,std::shared_ptr<SerializableBase> type) {
    ASSERT_MSG(!alias->type().empty(),"Given type doesn't support shared_ptr serialization!");
    if(_tSet.find(alias->type()) != _tSet.end()) {
      ASSERT_MSGV(typeid(*_tSet[alias->type()]) == typeid(*type),"Conflicit type id: %s",alias->type().c_str())
    }
    _tSet[alias->type()]=type;
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
void registerTypeAs(IOData* dat,std::shared_ptr<SerializableBase> alias,std::shared_ptr<SerializableBase> type)
{
  dat->registerTypeAs(alias,type);
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

//io for float is double
std::ostream& writeBinaryData(scalarF val,std::ostream& os,IOData*)
{
  double valD=val;
  os.write((char*)&valD,sizeof(double));
  return os;
}
std::istream& readBinaryData(scalarF& val,std::istream& is,IOData*)
{
  double valD;
  is.read((char*)&valD,sizeof(double));
  val=(scalarF)valD;
  return is;
}
std::ostream& writeBinaryData(scalarD val,std::ostream& os,IOData*)
{
  double valD=std::convert<double>()(val);
  os.write((char*)&valD,sizeof(double));
  return os;
}
std::istream& readBinaryData(scalarD& val,std::istream& is,IOData*)
{
  double valD;
  is.read((char*)&valD,sizeof(double));
  val=scalarD(valD);
  return is;
}

PRJ_END
