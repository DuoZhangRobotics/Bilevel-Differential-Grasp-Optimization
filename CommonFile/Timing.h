#ifndef TIMING_H
#define TIMING_H

#include "Config.h"
#include <iostream>
#include <tinyxml2.h>

PRJ_BEGIN

//timing
void saveTimingSetting();
void loadTimingSetting();
void enableTiming();
void disableTiming();
std::string TGETNAME(const std::string& path,const std::string& pathFrm,tinyxml2::XMLElement& pt);
void TFRMRESET(const std::string& pathFrm,tinyxml2::XMLElement& pt);
void TFRMADVANCE(const std::string& pathFrm,tinyxml2::XMLElement& pt);

void TBEG();
void TBEG(const std::string& name);
void TENDT(const std::string& path,tinyxml2::XMLElement& pt);
void TENDT(const std::string& path,const std::string& pathFrm,tinyxml2::XMLElement& pt);
void TEND(std::ostream& os);
void TEND();
scalarD TENDV();

PRJ_END

#endif
