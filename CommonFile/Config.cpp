#include "Config.h"
#include <cstdio>
#include <string>
#include <cstdarg>
#include <stdexcept>
#include <iostream>
#include <assert.h>

//coutPrintf
bool isInt(char c)
{
  return c == 'd';
}
bool isUInt(char c)
{
  return c == 'u';
}
bool isFloat(char c)
{
  return c == 'f' || c == 'F';
}
bool isString(char c)
{
  return c == 's';
}
bool isType(char c)
{
  return c == 'd' || c == 'i' || c == 'u' || c == 'o' ||
         c == 'x' || c == 'X' || c == 'f' || c == 'F' ||
         c == 'e' || c == 'E' || c == 'g' || c == 'G' ||
         c == 'a' || c == 'A' || c == 'c' || c == 's' ||
         c == 'p' || c == 'n';
}
void coutPrintf(std::string format,...)
{
  va_list args;
  va_start(args,format);
  while(format[0] != '[')
    format.erase(format.begin());
  while(format.back() != '\n')
    format.pop_back();
  for(size_t i=0; i<format.size();)
    if(format[i] == '%') {
      size_t j=i+1;
      while(j<format.size())
        if(isType(format[j])) {
          if(isFloat(format[j])) {
            std::cout << va_arg(args,scalarD);
            break;
          } else if(isString(format[j])) {
            std::cout << va_arg(args,const char*);
            break;
          } else if(isInt(format[j])) {
            std::cout << va_arg(args,sizeType);
            break;
          } else if(isUInt(format[j])) {
            std::cout << va_arg(args,unsigned int);
            break;
          } else {
            assert(false);
          }
        } else j++;
      assert(j<format.size());
      i=j+1;
    } else {
      std::cout << format[i];
      i++;
    }
  std::cout << std::flush;
  va_end(args);
}
