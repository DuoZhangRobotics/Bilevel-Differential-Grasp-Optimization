#include "Utils.h"
#include <iostream>
#include <iomanip>
#include <CommonFile/IOBasic.h>

namespace std
{
double to_double(char a) {
  return a;
}
double to_double(unsigned char a) {
  return a;
}
double to_double(float a) {
  return a;
}
double to_double(double a) {
  return a;
}
double to_double(sizeType a) {
  return double(a);
}
}

PRJ_BEGIN

std::vector<std::string> split(const std::string& l,const std::string& sep)
{
  sizeType i,last=-1;
  std::vector<std::string> ret;
  for(i=0; i<(sizeType)l.size(); i++)
    if(std::find(sep.begin(),sep.end(),l[i])!=sep.end()) {
      if(last<i-1)
        ret.push_back(l.substr(last+1,i-last-1));
      last=i;
    }
  if(last<i-1)
    ret.push_back(l.substr(last+1,i-last-1));
  return ret;
}
bool beginsWith(const std::string& l,const std::string& s)
{
  return l.length()>=s.length() && l.substr(0,s.length())==s;
}
bool endsWith(const std::string& l,const std::string& s)
{
  return l.length()>=s.length() && l.substr(l.size()-s.length(),s.length())==s;
}
std::string toUpper(const std::string& l) {
  std::string ret=l;
  std::for_each(ret.begin(),ret.end(),[&](char c) {
    return std::toupper(c);
  });
  return ret;
}

//filesystem
bool notDigit(char c)
{
  return !std::isdigit(c);
}
bool lessDirByNumber(std::experimental::filesystem::v1::path A,std::experimental::filesystem::v1::path B)
{
  std::string strA=A.string();
  std::string::iterator itA=std::remove_if(strA.begin(),strA.end(),notDigit);
  strA.erase(itA,strA.end());

  std::string strB=B.string();
  std::string::iterator itB=std::remove_if(strB.begin(),strB.end(),notDigit);
  strB.erase(itB,strB.end());

  sizeType idA,idB;
  std::istringstream(strA) >> idA;
  std::istringstream(strB) >> idB;
  return idA < idB;
}
bool exists(const std::experimental::filesystem::v1::path& path)
{
  return std::experimental::filesystem::v1::exists(path);
}
void removeDir(const std::experimental::filesystem::v1::path& path)
{
  if(std::experimental::filesystem::v1::exists(path))
    try {
      std::experimental::filesystem::v1::remove_all(path);
    } catch(...) {}
}
void create(const std::experimental::filesystem::v1::path& path)
{
  if(!std::experimental::filesystem::v1::exists(path))
    try {
      std::experimental::filesystem::v1::create_directory(path);
    } catch(...) {}
}
void recreate(const std::experimental::filesystem::v1::path& path)
{
  removeDir(path);
  try {
    std::experimental::filesystem::v1::create_directory(path);
  } catch(...) {}
}
std::vector<std::experimental::filesystem::v1::path> files(const std::experimental::filesystem::v1::path& path)
{
  std::vector<std::experimental::filesystem::v1::path> ret;
  for(std::experimental::filesystem::v1::directory_iterator beg(path),end; beg!=end; beg++)
    if(std::experimental::filesystem::v1::is_regular_file(*beg))
      ret.push_back(*beg);
  return ret;
}
std::vector<std::experimental::filesystem::v1::path> directories(const std::experimental::filesystem::v1::path& path)
{
  std::vector<std::experimental::filesystem::v1::path> ret;
  for(std::experimental::filesystem::v1::directory_iterator beg(path),end; beg!=end; beg++)
    if(std::experimental::filesystem::v1::is_directory(*beg))
      ret.push_back(*beg);
  return ret;
}
void sortFilesByNumber(std::vector<std::experimental::filesystem::v1::path>& files)
{
  sort(files.begin(),files.end(),lessDirByNumber);
}
bool isDir(const std::experimental::filesystem::v1::path& path)
{
  return std::experimental::filesystem::v1::is_directory(path);
}
size_t fileSize(const std::experimental::filesystem::v1::path& path)
{
  return (size_t)std::experimental::filesystem::v1::file_size(path);
}

//put/get
const tinyxml2::XMLElement* getChild(const tinyxml2::XMLElement& pt,const std::string& name)
{
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement())
    if(v->Name()==name)
      return v;
  return NULL;
}
tinyxml2::XMLElement* getChild(tinyxml2::XMLElement& pt,const std::string& name)
{
  for(tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement())
    if(v->Name()==name)
      return v;
  return NULL;
}
tinyxml2::XMLElement* addChild(tinyxml2::XMLElement& pt,const std::string& name)
{
  //for(tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement())
  //  if(v->Name()==name)
  //    return v;
  tinyxml2::XMLElement* node=pt.GetDocument()->NewElement(name.c_str());
  pt.InsertEndChild(node);
  return node;
}
const tinyxml2::XMLElement* getAttributeInfo(const tinyxml2::XMLElement& pt,std::string& name)
{
  const tinyxml2::XMLElement* curr=&pt;
  std::vector<std::string> paths=split(name,".");
  for(sizeType i=0; i<(sizeType)paths.size(); i++)
    if(paths[i]=="<xmlattr>") {
      ASSERT(i==(sizeType)paths.size()-2)
      name=paths[i+1];
      return curr;
    } else {
      curr=getChild(*curr,paths[i]);
      if(!curr)
        return curr;
    }
  name="";
  return curr;
}
tinyxml2::XMLElement* getAttributeInfoPut(tinyxml2::XMLElement& pt,std::string& name)
{
  tinyxml2::XMLElement* curr=&pt;
  std::vector<std::string> paths=split(name,".");
  for(sizeType i=0; i<(sizeType)paths.size(); i++)
    if(paths[i]=="<xmlattr>") {
      ASSERT(i==(sizeType)paths.size()-2)
      name=paths[i+1];
      return curr;
    } else {
      tinyxml2::XMLElement* tmp=getChild(*curr,paths[i]);
      if(!tmp) {
        tinyxml2::XMLElement* node=pt.GetDocument()->NewElement(paths[i].c_str());
        curr->InsertEndChild(node);
        curr=node;
      } else {
        curr=tmp;
      }
    }
  name="";
  return curr;
}
bool hasAttribute(const tinyxml2::XMLElement& pt,const std::string& name)
{
  std::string nameProcessed=name;
  const tinyxml2::XMLElement* e=getAttributeInfo(pt,nameProcessed);
  if(!e)
    return false;
  else if(nameProcessed.empty())
    return true;
  else return e->Attribute(nameProcessed.c_str())!=NULL;
}

//parsePtree
std::vector<std::string> toParams(sizeType argc,char** argv)
{
  std::vector<std::string> params;
  for(sizeType i=0; i<argc; i++)
    params.push_back(argv[i]);
  return params;
}
std::string parseProps(sizeType argc,char** argv,tinyxml2::XMLElement& pt)
{
  return parseProps(toParams(argc,argv),pt);
}
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLElement& pt)
{
  std::string addParam;
  for(sizeType i=0; i<(sizeType)params.size(); i++) {
    const std::string& str=params[i];
    size_t pos=str.find("=");
    if(pos != std::string::npos) {
      std::string LHS=str.substr(0,pos);
      std::string RHS=str.substr(pos+1);
      put<std::string>(pt,LHS.c_str(),RHS.c_str());
      addParam+="_"+str;
    }
  }
  return addParam+"_";
}
std::string parseProps(sizeType argc,char** argv,tinyxml2::XMLDocument& pt)
{
  return parseProps(argc,argv,*(pt.RootElement()));
}
std::string parseProps(const std::vector<std::string>& params,tinyxml2::XMLDocument& pt)
{
  return parseProps(params,*(pt.RootElement()));
}

PRJ_END
