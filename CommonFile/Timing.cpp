#include "Timing.h"
#include "MathBasic.h"
#include <vector>
#include <chrono>

PRJ_BEGIN

bool __timing=true;
std::vector<char> __timingSetting;
std::vector<std::pair<std::string,std::chrono::time_point<std::chrono::steady_clock>>> __timer;

std::string replace_all(std::string str,const std::string& from,const std::string& to) {
  size_t start_pos=0;
  while((start_pos=str.find(from,start_pos))!=std::string::npos) {
    str.replace(start_pos,from.length(),to);
    start_pos+=to.length(); // Handles case where 'to' is a substring of 'from'
  }
  return str;
}
void saveTimingSetting()
{
  __timingSetting.push_back(__timing);
}
void loadTimingSetting()
{
  ASSERT_MSG(!__timingSetting.empty(),"loading timing setting with no setting available!")
  __timing=(bool)__timingSetting.back();
}
void enableTiming()
{
  __timing=true;
}
void disableTiming()
{
  __timing=false;
}
std::string TGETNAME(const std::string& path,const std::string& pathFrm,tinyxml2::XMLElement& pt)
{
  int def=0;
  std::string pathReal=path;
  int fid=pt.QueryIntAttribute(pathFrm.c_str(),&def);
  replace_all(pathReal,"#",std::to_string(fid));
  return pathReal;
}
void TFRMRESET(const std::string& pathFrm,tinyxml2::XMLElement& pt)
{
  pt.SetAttribute(pathFrm.c_str(),0);
}
void TFRMADVANCE(const std::string& pathFrm,tinyxml2::XMLElement& pt)
{
  int def=0;
  int fid=pt.QueryIntAttribute(pathFrm.c_str(),&def);
  pt.SetAttribute(pathFrm.c_str(),fid+1);
}

void TBEG()
{
  TBEG("");
}
void TBEG(const std::string& name)
{
  OMP_CRITICAL_ {
    __timer.push_back(make_pair(name,std::chrono::steady_clock::now()));
  }
}
void TENDT(const std::string& path,tinyxml2::XMLElement& pt)
{
  OMP_CRITICAL_ {
    ASSERT_MSG(!__timer.empty(),"Incorrect calling of TEND without TBEG!")
    std::chrono::time_point<std::chrono::steady_clock> ptt=std::chrono::steady_clock::now();
    std::chrono::duration<double> diff=ptt-__timer.back().second;
    pt.SetAttribute(path.c_str(),diff.count());
    __timer.pop_back();
  }
}
void TENDT(const std::string& path,const std::string& pathFrm,tinyxml2::XMLElement& pt)
{
  TENDT(TGETNAME(path,pathFrm,pt),pt);
}
void TEND(std::ostream& os)
{
  OMP_CRITICAL_ {
    ASSERT_MSG(!__timer.empty(),"Incorrect calling of TEND without TBEG!")
    if(__timing) {
      std::chrono::time_point<std::chrono::steady_clock> pt=std::chrono::steady_clock::now();
      std::chrono::duration<double> diff=pt-__timer.back().second;
      os << __timer.back().first << ": " << diff.count();
    }
    __timer.pop_back();
  }
}
void TEND()
{
  TEND(std::cout);
  OMP_CRITICAL_ {
    if(__timing)
      std::cout << std::endl;
  }
}
scalarD TENDV()
{
  scalarD ret;
  OMP_CRITICAL_ {
    ASSERT_MSG(!__timer.empty(),"Incorrect calling of TEND without TBEG!")
    std::chrono::time_point<std::chrono::steady_clock> pt=std::chrono::steady_clock::now();
    std::chrono::duration<double> diff=pt-__timer.back().second;
    ret=diff.count();
    __timer.pop_back();
  }
  return ret;
}

PRJ_END
