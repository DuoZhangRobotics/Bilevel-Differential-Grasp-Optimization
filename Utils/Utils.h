#ifndef UTILS_H
#define UTILS_H

#include <experimental/filesystem>
#include <tinyxml2.h>
#include "Scalar.h"

PRJ_BEGIN

extern std::vector<std::string> split(const std::string& l,const std::string& sep=" ");
extern bool beginsWith(const std::string& l,const std::string& s);
extern bool endsWith(const std::string& l,const std::string& s);
extern std::string toUpper(const std::string& l);

//filesystem
bool exists(const std::experimental::filesystem::v1::path& path);
void removeDir(const std::experimental::filesystem::v1::path& path);
void create(const std::experimental::filesystem::v1::path& path);
void recreate(const std::experimental::filesystem::v1::path& path);
std::vector<std::experimental::filesystem::v1::path> files(const std::experimental::filesystem::v1::path& path);
std::vector<std::experimental::filesystem::v1::path> directories(const std::experimental::filesystem::v1::path& path);
void sortFilesByNumber(std::vector<std::experimental::filesystem::v1::path>& files);
bool isDir(const std::experimental::filesystem::v1::path& path);
size_t fileSize(const std::experimental::filesystem::v1::path& path);

//property tree
template <typename SCALAR>
struct PtreeSafeType {
  typedef SCALAR SAFE_SCALAR;
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const SCALAR& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(SAFE_SCALAR(val));
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),SAFE_SCALAR(val));
  }
  static SCALAR get(const tinyxml2::XMLElement& pt,const std::string& name,const SCALAR& val)
  {
    SAFE_SCALAR def=val;
    if(name.empty())
      pt.QueryFloatText(&def);
    else pt.QueryFloatAttribute(name.c_str(),&def);
    return def;
  }
  static SCALAR get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      ASSERT(pt.QueryFloatText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryFloatAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<sizeType> {
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const sizeType& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(val);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),(int)val);
  }
  static sizeType get(const tinyxml2::XMLElement& pt,const std::string& name,const sizeType& val)
  {
    int def=val;
    if(name.empty())
      pt.QueryIntText(&def);
    else pt.QueryIntAttribute(name.c_str(),&def);
    return def;
  }
  static sizeType get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    int def=0;
    if(name.empty()) {
      ASSERT(pt.QueryIntText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryIntAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<int> {
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const int& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(val);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),val);
  }
  static int get(const tinyxml2::XMLElement& pt,const std::string& name,const int& val)
  {
    int def=val;
    if(name.empty())
      pt.QueryIntText(&def);
    else pt.QueryIntAttribute(name.c_str(),&def);
    return def;
  }
  static int get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    int def=0;
    if(name.empty()) {
      ASSERT(pt.QueryIntText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryIntAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<bool> {
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const bool& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(val);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),val);
  }
  static bool get(const tinyxml2::XMLElement& pt,const std::string& name,const bool& val)
  {
    bool def=val;
    if(name.empty())
      pt.QueryBoolText(&def);
    else pt.QueryBoolAttribute(name.c_str(),&def);
    return def;
  }
  static bool get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    bool def=true;
    if(name.empty()) {
      ASSERT(pt.QueryBoolText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryBoolAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<char> {
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const char& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(val);
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),val);
  }
  static char get(const tinyxml2::XMLElement& pt,const std::string& name,const char& val)
  {
    int def=val;
    if(name.empty())
      pt.QueryIntText(&def);
    else pt.QueryIntAttribute(name.c_str(),&def);
    return def;
  }
  static char get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    int def=0;
    if(name.empty()) {
      ASSERT(pt.QueryIntText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryIntAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<scalarD> {
  typedef scalarF SAFE_SCALAR;
  //text
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const scalarD& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(std::to_double(val));
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),std::to_double(val));
  }
  static scalarD get(const tinyxml2::XMLElement& pt,const std::string& name,const scalarD& val)
  {
    SAFE_SCALAR def=val;
    if(name.empty())
      pt.QueryFloatText(&def);
    else pt.QueryFloatAttribute(name.c_str(),&def);
    return def;
  }
  static scalarD get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      ASSERT(pt.QueryFloatText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryFloatAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<__float128> {
  typedef scalarF SAFE_SCALAR;
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const __float128& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(std::to_double(val));
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),std::to_double(val));
  }
  static __float128 get(const tinyxml2::XMLElement& pt,const std::string& name,const __float128& val)
  {
    SAFE_SCALAR def=val;
    if(name.empty())
      pt.QueryFloatText(&def);
    else pt.QueryFloatAttribute(name.c_str(),&def);
    return def;
  }
  static __float128 get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      ASSERT(pt.QueryFloatText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryFloatAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<mpfr::mpreal> {
  typedef scalarF SAFE_SCALAR;
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const mpfr::mpreal& val)
  {
    if(name.empty()) {
      std::string text=std::to_string(std::to_double(val));
      pt.SetText(text.c_str());
    } else pt.SetAttribute(name.c_str(),SAFE_SCALAR(std::to_double(val)));
  }
  static mpfr::mpreal get(const tinyxml2::XMLElement& pt,const std::string& name,const mpfr::mpreal& val)
  {
    SAFE_SCALAR def=std::to_double(val);
    if(name.empty())
      pt.QueryFloatText(&def);
    else pt.QueryFloatAttribute(name.c_str(),&def);
    return def;
  }
  static mpfr::mpreal get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    SAFE_SCALAR def=0;
    if(name.empty()) {
      ASSERT(pt.QueryFloatText(&def)==tinyxml2::XML_SUCCESS)
    } else {
      ASSERT(pt.QueryFloatAttribute(name.c_str(),&def)==tinyxml2::XML_SUCCESS)
    }
    return def;
  }
};
template <>
struct PtreeSafeType<std::string> {
  static void put(tinyxml2::XMLElement& pt,const std::string& name,const std::string& val)
  {
    if(name.empty()) {
      pt.SetText(val.c_str());
    } else pt.SetAttribute(name.c_str(),val.c_str());
  }
  static std::string get(const tinyxml2::XMLElement& pt,const std::string& name,const std::string& val)
  {
    if(name.empty())
      return pt.GetText();
    else {
      const char* ret=pt.Attribute(name.c_str());
      if(!ret)
        return val;
      else return ret;
    }
  }
  static std::string get(const tinyxml2::XMLElement& pt,const std::string& name)
  {
    const char* ret=NULL;
    if(name.empty())
      ret=pt.GetText();
    else ret=pt.Attribute(name.c_str());
    ASSERT(ret!=0)
    return ret;
  }
};

//put/get
const tinyxml2::XMLElement* getChild(const tinyxml2::XMLElement& pt,const std::string& name);
tinyxml2::XMLElement* getChild(tinyxml2::XMLElement& pt,const std::string& name);
tinyxml2::XMLElement* addChild(tinyxml2::XMLElement& pt,const std::string& name);
const tinyxml2::XMLElement* getAttributeInfo(const tinyxml2::XMLElement& pt,std::string& name);
tinyxml2::XMLElement* getAttributeInfoPut(tinyxml2::XMLElement& pt,std::string& name);
bool hasAttribute(const tinyxml2::XMLElement& pt,const std::string& name);
template <typename T>
void put(tinyxml2::XMLElement& pt,const std::string& name,const T& val)
{
  std::string nameProcessed=name;
  tinyxml2::XMLElement* e=getAttributeInfoPut(pt,nameProcessed);
  ASSERT(e)
  PtreeSafeType<T>::put(*e,nameProcessed,val);
}
template <typename T>
void put(tinyxml2::XMLDocument& pt,const std::string& name,const T& val)
{
  put<T>(*(pt.RootElement()),name,val);
}
template <typename T>
void putCond(tinyxml2::XMLElement& pt,const std::string& path,T val)
{
  if(!hasAttribute(pt,path))
    put<T>(pt,path,val);
}
template <typename T>
void putCond(tinyxml2::XMLDocument& pt,const std::string& path,T val)
{
  putCond<T>(*(pt.RootElement()),path,val);
}
template <typename T>
T get(const tinyxml2::XMLElement& pt,const std::string& name,const T& val)
{
  std::string nameProcessed=name;
  const tinyxml2::XMLElement* e=getAttributeInfo(pt,nameProcessed);
  if(!e)
    return val;
  else return PtreeSafeType<T>::get(*e,nameProcessed,val);
}
template <typename T>
T get(const tinyxml2::XMLDocument& pt,const std::string& name,const T& val)
{
  return get<T>(*(pt.RootElement()),name,val);
}
template <typename T>
T get(const tinyxml2::XMLElement& pt,const std::string& name)
{
  std::string nameProcessed=name;
  const tinyxml2::XMLElement* e=getAttributeInfo(pt,nameProcessed);
  ASSERT(e)
  return PtreeSafeType<T>::get(*e,nameProcessed);
}
template <typename T>
T get(const tinyxml2::XMLDocument& pt,const std::string& name)
{
  return get<T>(*(pt.RootElement()),name);
}

//parsePtree
std::vector<std::string> toParams(sizeType argc,char** argv);
std::string parseProps(sizeType argc,char** argv,tinyxml2::XMLElement& pt);
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLElement& pt);
std::string parseProps(sizeType argc,char** argv,tinyxml2::XMLDocument& pt);
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLDocument& pt);
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLElement& pt,const std::string& path,const std::string& def,sizeType r=-1,sizeType c=-1)
{
  EIGEN_VEC ret;
  if(r<0&&c<0)
    ret.setZero();
  else if(r>=0&&c<0)
    ret.setZero(r);
  else ret.setZero(r,c);
  //read
  std::string str=get<std::string>(pt,path.c_str(),def.c_str());
  std::vector<std::string> strs=split(str," _,");
  sizeType newSz=0;
  for(sizeType j=0; j<(sizeType)strs.size(); j++)
    if(!strs[j].empty())
      strs[newSz++]=strs[j];
  strs.resize(newSz);
  if(str.empty()) //fixing boost bug
    strs.clear();
  ASSERT_MSGV(ret.size()==0 || (sizeType)strs.size() == ret.size(),
              "Try reading: %s, Vec size mismatch, size(strs)=%ld size(ret)=%ld!",
              path.c_str(),strs.size(),ret.size())
  if(ret.size() == 0)
    ret.resize((sizeType)strs.size());
  for(sizeType i=0; i<ret.size(); i++) {
    typename PtreeSafeType<typename EIGEN_VEC::Scalar>::SAFE_SCALAR val;
    std::istringstream(strs[i]) >> val;
    ret[i]=val;
  }
  return ret;
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLDocument& pt,const std::string& path,const std::string& def,sizeType r=-1,sizeType c=-1)
{
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,def,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLElement& pt,const std::string& path,const EIGEN_VEC& ret,sizeType r=-1,sizeType c=-1)
{
  std::string val;
  for(sizeType i=0; i<ret.size(); i++)
    if(i == 0)
      val+=std::to_string((typename PtreeSafeType<typename EIGEN_VEC::Scalar>::SAFE_SCALAR)ret[i]);
    else val+=" "+std::to_string((typename PtreeSafeType<typename EIGEN_VEC::Scalar>::SAFE_SCALAR)ret[i]);
  return parsePtreeDef<EIGEN_VEC>(pt,path,val,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLDocument& pt,const std::string& path,const EIGEN_VEC& ret,sizeType r=-1,sizeType c=-1)
{
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,ret,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLElement& pt,const std::string& path,typename EIGEN_VEC::Scalar val,sizeType r=-1,sizeType c=-1)
{
  EIGEN_VEC ret;
  if(r<0&&c<0)
    ret.setZero();
  else if(r>=0&&c<0)
    ret.setZero(r);
  else ret.setZero(r,c);
  ret.setConstant(val);
  return parsePtreeDef<EIGEN_VEC>(pt,path,ret,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtreeDef(const tinyxml2::XMLDocument& pt,const std::string& path,typename EIGEN_VEC::Scalar val,sizeType r=-1,sizeType c=-1)
{
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,val,r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtree(const tinyxml2::XMLElement& pt,const std::string& path,sizeType r=-1,sizeType c=-1)
{
  return parsePtreeDef<EIGEN_VEC>(pt,path,"",r,c);
}
template <typename EIGEN_VEC>
EIGEN_VEC parsePtree(const tinyxml2::XMLDocument& pt,const std::string& path,sizeType r=-1,sizeType c=-1)
{
  return parsePtreeDef<EIGEN_VEC>(*(pt.RootElement()),path,"",r,c);
}
template <typename EIGEN_VEC>
void putPtree(tinyxml2::XMLElement& pt,const std::string& path,const EIGEN_VEC& ret)
{
  std::string val;
  for(sizeType i=0; i<ret.size(); i++)
    if(i == 0)
      val+=std::to_string((typename PtreeSafeType<typename EIGEN_VEC::Scalar>::SAFE_SCALAR)ret[i]);
    else val+=" "+std::to_string((typename PtreeSafeType<typename EIGEN_VEC::Scalar>::SAFE_SCALAR)ret[i]);
  put<std::string>(pt,path.c_str(),val.c_str());
}
template <typename EIGEN_VEC>
void putPtree(tinyxml2::XMLDocument& pt,const std::string& path,const EIGEN_VEC& ret)
{
  putPtree<EIGEN_VEC>(*(pt.RootElement()),path,ret);
}
template <typename EIGEN_VEC>
std::string toStringPtree(const EIGEN_VEC& v)
{
  std::ostringstream oss;
  for(sizeType i=0; i<v.size(); i++) {
    oss << v[i];
    if(i<v.size()-1)
      oss << ",";
  }
  return oss.str();
}
template <typename T>
std::string toStringPtree(const std::vector<T>& v)
{
  std::ostringstream oss;
  for(sizeType i=0; i<(sizeType)v.size(); i++) {
    oss << v[i];
    if(i<(sizeType)v.size()-1)
      oss << ",";
  }
  return oss.str();
}

PRJ_END

#endif
