#ifndef OPTIONS_H
#define OPTIONS_H

#include <Utils/Scalar.h>
#include <functional>
#include <cxxabi.h>
#include <iomanip>
#include <sstream>
#include <map>
#include <set>

PRJ_BEGIN

template <typename TVar,typename T>
struct FunctionWrapper : std::function<TVar(void*,TVar*)>
{
  FunctionWrapper(std::function<TVar(T&,TVar*)> inner):_inner(inner) {}
  TVar operator()(void* ptr,TVar* val) {
    return _inner(*(T*)ptr,val);
  }
  std::function<TVar(T&,TVar*)> _inner;
};
class Options
{
public:
  template <typename T>
  void registerIntOption(const std::string& name,std::function<int(T&,int*)> f) {
    ASSERT(_allNames.find(name)==_allNames.end())
    _intOption[demangleTypeName<T>()][name]=FunctionWrapper<int,T>(f);
    _allNames.insert(std::string(demangleTypeName<T>())+"."+name);
    _allTypes.insert(demangleTypeName<T>());
  }
  template <typename T>
  void registerBoolOption(const std::string& name,std::function<bool(T&,bool*)> f) {
    ASSERT(_allNames.find(name)==_allNames.end())
    _boolOption[demangleTypeName<T>()][name]=FunctionWrapper<bool,T>(f);
    _allNames.insert(std::string(demangleTypeName<T>())+"."+name);
    _allTypes.insert(demangleTypeName<T>());
  }
  template <typename T>
  void registerFloatOption(const std::string& name,std::function<double(T&,double*)> f) {
    ASSERT(_allNames.find(name)==_allNames.end())
    _floatOption[demangleTypeName<T>()][name]=FunctionWrapper<double,T>(f);
    _allNames.insert(std::string(demangleTypeName<T>())+"."+name);
    _allTypes.insert(demangleTypeName<T>());
  }
  template <typename T>
  void registerFloatQOption(const std::string& name,std::function<__float128(T&,__float128*)> f) {
    ASSERT(_allNames.find(name)==_allNames.end())
    _floatQOption[demangleTypeName<T>()][name]=FunctionWrapper<__float128,T>(f);
    _allNames.insert(std::string(demangleTypeName<T>())+"."+name);
    _allTypes.insert(demangleTypeName<T>());
  }
  template <typename T>
  void registerMPFROption(const std::string& name,std::function<mpfr::mpreal(T&,mpfr::mpreal*)> f) {
    ASSERT(_allNames.find(name)==_allNames.end())
    _MPFROption[demangleTypeName<T>()][name]=FunctionWrapper<mpfr::mpreal,T>(f);
    _allNames.insert(std::string(demangleTypeName<T>())+"."+name);
    _allTypes.insert(demangleTypeName<T>());
  }
  template <typename TYPE,typename T>
  void setOptions(const std::string& name,T val) {
    std::string fullname=std::string(demangleTypeName<TYPE>())+"."+name;
    ASSERT_MSGV(_allNames.find(fullname)!=_allNames.end(),"Unknown parameter: %s",name.c_str())
    _vss[demangleTypeName<TYPE>()][name]=std::to_string(val);
  }
  template <typename T>
  bool hasType() {
    return _allTypes.find(demangleTypeName<T>())!=_allTypes.end();
  }
  template <typename T>
  void setOptions(T* t,bool fetchDefault=true) {
    if(_intOption.find(demangleTypeName<T>())!=_intOption.end())
      setOptions<T,int>(t,_intOption.find(demangleTypeName<T>())->second,fetchDefault);
    if(_boolOption.find(demangleTypeName<T>())!=_boolOption.end())
      setOptions<T,bool>(t,_boolOption.find(demangleTypeName<T>())->second,fetchDefault);
    if(_floatOption.find(demangleTypeName<T>())!=_floatOption.end())
      setOptions<T,double>(t,_floatOption.find(demangleTypeName<T>())->second,fetchDefault);
    if(_floatQOption.find(demangleTypeName<T>())!=_floatQOption.end())
      setOptions<T,__float128>(t,_floatQOption.find(demangleTypeName<T>())->second,fetchDefault);
    if(_MPFROption.find(demangleTypeName<T>())!=_MPFROption.end())
      setOptions<T,mpfr::mpreal>(t,_MPFROption.find(demangleTypeName<T>())->second,fetchDefault);
  }
  template <typename T,typename T2>
  void setOptions(T* t,const std::map<std::string,std::function<T2(void*,T2*)>>& ops,bool fetchDefault) {
    for(typename std::map<std::string,std::function<T2(void*,T2*)>>::const_iterator beg=ops.begin(),end=ops.end(); beg!=end; beg++) {
      std::map<std::string,std::map<std::string,std::string>>::const_iterator it=_vss.find(demangleTypeName<T>());
      if(it!=_vss.end()) {
        std::map<std::string,std::string>::const_iterator it2=it->second.find(beg->first);
        if(it2!=it->second.end()) {
          T2 t2;
          std::istringstream(it2->second) >> t2;
          beg->second(t,&t2);
        } else if(fetchDefault) {
          std::ostringstream iss;
          iss << beg->second(t,NULL);
          _vss[demangleTypeName<T>()][beg->first]=iss.str();
        }
      } else if(fetchDefault) {
        std::ostringstream iss;
        iss << beg->second(t,NULL);
        _vss[demangleTypeName<T>()][beg->first]=iss.str();
      }
    }
  }
  void print() const
  {
    sizeType maxWidthL=0;
    sizeType maxWidthR=0;
    for(std::map<std::string,std::map<std::string,std::string>>::const_iterator beg=_vss.begin(),end=_vss.end(); beg!=end; beg++) {
      for(std::map<std::string,std::string>::const_iterator beg2=beg->second.begin(),end2=beg->second.end(); beg2!=end2; beg2++) {
        maxWidthL=std::max<sizeType>(maxWidthL,beg->first.size()+beg2->first.size()+1);
        maxWidthR=std::max<sizeType>(maxWidthR,beg2->second.size());
      }
    }
    INFO((std::string("[Options]")+std::string(maxWidthL+maxWidthR+2,'-')).c_str());
    for(std::map<std::string,std::map<std::string,std::string>>::const_iterator beg=_vss.begin(),end=_vss.end(); beg!=end; beg++) {
      for(std::map<std::string,std::string>::const_iterator beg2=beg->second.begin(),end2=beg->second.end(); beg2!=end2; beg2++) {
        std::ostringstream oss;
        oss << "[Options]" << std::left << std::setw(maxWidthL) << (beg->first+"."+beg2->first) << ": " << std::left << std::setw(maxWidthR) << beg2->second;
        INFO(oss.str().c_str())
      }
      INFO((std::string("[Options]")+std::string(maxWidthL+maxWidthR+2,'-')).c_str());
    }
  }
private:
  std::string demangle(const char* mangled)
  {
    int status;
    std::unique_ptr<char[],void(*)(void*)> result(abi::__cxa_demangle(mangled,0,0,&status),std::free);
    return result.get() ? std::string(result.get()) : "error occurred";
  }
  template <typename T>
  std::string demangleTypeName()
  {
    return demangle(typeid(T).name());
  }
  std::map<std::string,std::map<std::string,std::string>> _vss;
  std::set<std::string> _allNames,_allTypes;
  //type->name->function
  std::map<std::string,std::map<std::string,std::function<int(void*,int*)>>> _intOption;
  std::map<std::string,std::map<std::string,std::function<bool(void*,bool*)>>> _boolOption;
  std::map<std::string,std::map<std::string,std::function<double(void*,double*)>>> _floatOption;
  std::map<std::string,std::map<std::string,std::function<__float128(void*,__float128*)>>> _floatQOption;
  std::map<std::string,std::map<std::string,std::function<mpfr::mpreal(void*,mpfr::mpreal*)>>> _MPFROption;
};
#define REGISTER_INT_TYPE(NAME,T,TMEMBER,EXPR)    \
ops.registerIntOption<T>(NAME,[](T& t,int* val) { \
  if(val)   \
    EXPR=(TMEMBER)*val; \
  return EXPR;  \
});
#define REGISTER_BOOL_TYPE(NAME,T,TMEMBER,EXPR)    \
ops.registerBoolOption<T>(NAME,[](T& t,bool* val) { \
  if(val)   \
    EXPR=(TMEMBER)*val; \
  return EXPR;  \
});
#define REGISTER_FLOAT_TYPE(NAME,T,TMEMBER,EXPR)    \
ops.registerFloatOption<T>(NAME,[](T& t,double* val) { \
  if(val)   \
    EXPR=(TMEMBER)*val; \
  return EXPR;  \
});
#define REGISTER_FLOAT128_TYPE(NAME,T,TMEMBER,EXPR)    \
ops.registerFloatQOption<T>(NAME,[](T& t,__float128* val) { \
  if(val)   \
    EXPR=(TMEMBER)*val; \
  return EXPR;  \
});
#define REGISTER_MPFR_TYPE(NAME,T,TMEMBER,EXPR)    \
ops.registerMPFROption<T>(NAME,[](T& t,mpfr::mpreal* val) { \
  if(val)   \
    EXPR=(TMEMBER)*val; \
  return EXPR;  \
});

PRJ_END

#endif
