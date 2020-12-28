#include "MathBasic.h"
#include <iostream>
#include <random>
#include <ctime>
#include <map>

PRJ_BEGIN

#define DECL_FLOAT_TRAITS(T)  \
bool compL(T a,T b) {return a<b;} \
bool compLE(T a,T b) {return a<=b;} \
bool compG(T a,T b)  {return a>b;}  \
bool compGE(T a,T b) {return a>=b;} \
T compMin(T a,T b) {return std::min(a,b);}  \
T compMax(T a,T b) {return std::max(a,b);}  \
T floorV(T a) {return std::floor(a);} \
T ceilV(T a) {return std::ceil(a);}
DECL_FLOAT_TRAITS(char)
DECL_FLOAT_TRAITS(unsigned char)
DECL_FLOAT_TRAITS(scalarF)
DECL_FLOAT_TRAITS(scalarD)
DECL_FLOAT_TRAITS(sizeType)
#undef DECL_FLOAT_TRAITS

//RandEngine
//#define USE_BOOST_RANDOM
bool USEDeterministic=false;
std::random_device vgenDevice;
struct STDRandom {
  typedef int result_type;
  typedef std::map<sizeType,std::pair<sizeType,std::vector<int> > > record;
  STDRandom():_id(-1) {}
  static int min() {
    return 0;
  }
  static int max() {
    return RAND_MAX;
  }
  int operator()() {
    int ret=0;
    if(_id >= 0) {
      std::pair<sizeType,std::vector<int> >& pair=_record[_id];
      while(pair.first >= (sizeType)pair.second.size())
#ifdef USE_BOOST_RANDOM
        pair.second.push_back(vgen());
#else
        pair.second.push_back(rand());
#endif
      ret=pair.second[pair.first++];
    } else {
#ifdef USE_BOOST_RANDOM
      ret=vgen();
#else
      ret=rand();
#endif
    }
    return ret;
  }
  void seed(unsigned int s) {
#ifdef USE_BOOST_RANDOM
    vgen.seed(s);
#else
    srand(s);
#endif
  }
  std::mt19937 vgen;
  record _record;
  sizeType _id;
} vgenPesudo;
//real
scalar RandEngine::randR01()
{
  scalar ret=0;
  OMP_CRITICAL_
  if(USEDeterministic)
    ret=std::uniform_real_distribution<double>(0,1)(vgenPesudo);
  else ret=std::uniform_real_distribution<double>(0,1)(vgenDevice);
  return ret;
}
scalar RandEngine::randR(scalar l,scalar u)
{
  if(l == u)
    return l;
  else if(l > u)
    std::swap(l,u);
  scalar ret=0;
  OMP_CRITICAL_
  if(USEDeterministic)
    ret=std::uniform_real_distribution<double>(l,u)(vgenPesudo);
  else ret=std::uniform_real_distribution<double>(l,u)(vgenDevice);
  return ret;
}
scalar RandEngine::randSR(sizeType seed,scalar l,scalar u)
{
  if(l == u)
    return l;
  else if(l > u)
    std::swap(l,u);
  scalar ret=0;
  OMP_CRITICAL_ {
    vgenPesudo.seed((const uint32_t)seed);
    ret=std::uniform_real_distribution<double>(l,u)(vgenPesudo);
  }
  return ret;
}
//int
sizeType RandEngine::randSI(sizeType seed)
{
  return randSI(seed,0,std::numeric_limits<int>::max());
}
sizeType RandEngine::randI(sizeType l,sizeType u)
{
  if(l == u)
    return l;
  else if(l > u)
    std::swap(l,u);
  sizeType ret=0;
  OMP_CRITICAL_
  if(USEDeterministic)
    ret=std::uniform_int_distribution<sizeType>(l,u)(vgenPesudo);
  else ret=std::uniform_int_distribution<sizeType>(l,u)(vgenDevice);
  return ret;
}
sizeType RandEngine::randSI(sizeType seed,sizeType l,sizeType u)
{
  if(l == u)
    return l;
  else if(l > u)
    std::swap(l,u);
  sizeType ret=0;
  OMP_CRITICAL_ {
    vgenPesudo.seed((const uint32_t)seed);
    ret=std::uniform_int_distribution<sizeType>(l,u)(vgenPesudo);
  }
  return ret;
}
//normal
scalar RandEngine::normal(scalar mean,scalar cov)
{
  double ret=0;
  OMP_CRITICAL_
  if(USEDeterministic)
    ret=std::normal_distribution<double>(mean,cov)(vgenPesudo);
  else ret=std::normal_distribution<double>(mean,cov)(vgenDevice);
  return ret;
}
scalar RandEngine::normal()
{
  double ret=0;
  OMP_CRITICAL_
  if(USEDeterministic)
    ret=std::normal_distribution<double>()(vgenPesudo);
  else ret=std::normal_distribution<double>()(vgenDevice);
  return ret;
}
//settings
void RandEngine::seedTime()
{
  std::time_t result=std::time(nullptr);
  std::asctime(std::localtime(&result));
  seed(result);
}
void RandEngine::seed(sizeType i)
{
  vgenPesudo.seed((const uint32_t)i);
}
void RandEngine::useDeterministic()
{
  USEDeterministic=true;
}
void RandEngine::useNonDeterministic()
{
  USEDeterministic=false;
}
//common random number
const std::vector<int>& RandEngine::getRecordedRandom(sizeType id)
{
  ASSERT_MSGV(false,"Invalid random id: %ld!",id)
  STDRandom::record& rcd=vgenPesudo._record;
  return rcd[id].second;
}
void RandEngine::resetAllRecordRandom()
{
  STDRandom::record& rcd=vgenPesudo._record;
  for(STDRandom::record::iterator
      beg=rcd.begin(),end=rcd.end(); beg!=end; beg++)
    beg->second.first=0;
}
void RandEngine::resetRecordRandom(sizeType id)
{
  ASSERT_MSGV(id >= 0,"Invalid random id: %ld!",id)
  STDRandom::record& rcd=vgenPesudo._record;
  rcd[id].first=0;
}
void RandEngine::beginRecordRandom(sizeType id)
{
  ASSERT_MSGV(id >= 0,"Invalid random id: %ld!",id)
  STDRandom::record& rcd=vgenPesudo._record;
  if(rcd.find(id) == rcd.end())
    rcd[id]=make_pair((sizeType)0,std::vector<int>());
  vgenPesudo._id=id;
}
void RandEngine::endRecordRandom()
{
  vgenPesudo._id=-1;
}
void RandEngine::clearAllRecord()
{
  ASSERT_MSG(vgenPesudo._id == -1,"We must have vgenPesudo._id == -1 to modify cache!")
  STDRandom::record& rcd=vgenPesudo._record;
  rcd.clear();
}
void RandEngine::clearRecord(sizeType id)
{
  ASSERT_MSG(vgenPesudo._id == -1,"We must have vgenPesudo._id == -1 to modify cache!")
  STDRandom::record& rcd=vgenPesudo._record;
  if(rcd.find(id) != rcd.end())
    rcd.erase(id);
}
std::vector<sizeType> RandEngine::getRecordId()
{
  std::vector<sizeType> ret;
  STDRandom::record& rcd=vgenPesudo._record;
  for(STDRandom::record::const_iterator
      beg=rcd.begin(),end=rcd.end(); beg!=end; beg++)
    ret.push_back(beg->first);
  return ret;
}
void RandEngine::printRecordId()
{
  std::cout << "Random ids: ";
  STDRandom::record& rcd=vgenPesudo._record;
  for(STDRandom::record::const_iterator
      beg=rcd.begin(),end=rcd.end(); beg!=end; beg++)
    std::cout << beg->first << " ";
  std::cout << std::endl;
}

//OpenMP Settings
const OmpSettings& OmpSettings::getOmpSettings()
{
  return _ompSettings;
}
OmpSettings& OmpSettings::getOmpSettingsNonConst()
{
  return _ompSettings;
}
int OmpSettings::nrThreads() const
{
  return _nrThreads;
}
#ifdef NO_OPENMP
int OmpSettings::threadId() const
{
  return 0;
}
void OmpSettings::setNrThreads(int nr) {}
void OmpSettings::useAllThreads() {}
OmpSettings::OmpSettings():_nrThreads(1) {}
#else
int OmpSettings::threadId() const
{
  return omp_get_thread_num();
}
void OmpSettings::setNrThreads(int nr)
{
  nr=std::max<int>(nr,1);
  omp_set_num_threads(nr);
  _nrThreads=nr;
}
void OmpSettings::useAllThreads()
{
  _nrThreads=omp_get_num_procs();
}
OmpSettings::OmpSettings():_nrThreads(std::max<int>(omp_get_num_procs(),2)*3/4) {}
#endif
OmpSettings OmpSettings::_ompSettings;

PRJ_END
