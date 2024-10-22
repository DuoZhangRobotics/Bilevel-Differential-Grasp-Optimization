#ifndef CONFIG_H
#define CONFIG_H

#define NAMESPACE COMMON
#define PRJ_BEGIN namespace NAMESPACE {
#define PRJ_END	}
#define USE_PRJ_NAMESPACE using namespace NAMESPACE;

//#define _DEBUG
#ifdef _DEBUG
#define CUSTOM_DEBUG
#endif
#define PROFILE

#ifdef _MSC_VER
#define ALIGN_8 __declspec(align(8))
#define ALIGN_16 __declspec(align(16))
#define FORCE_INLINE __forceinline
#else
#define ALIGN_8 __attribute__((aligned (8)))
#define ALIGN_16 __attribute__((aligned (16)))
#define FORCE_INLINE inline
#endif

#include <string>
void coutPrintf(const std::string fmt,...);
//#define NDEBUG
#ifndef NDEBUG

#include <assert.h>
#define ASSERT(x) do{assert((x));}while(0);
#define ASSERT_MSG(x,msg) do{if(!(x)){coutPrintf("[ERROR] %s \n",msg);fflush(stdout);assert(false);}}while(0);
#define ASSERT_MSGV(x,fmt,...) do{if(!(x)){coutPrintf("[ERROR] " fmt " \n",__VA_ARGS__);fflush(stdout);assert(false);}}while(0);

#else

#ifdef _MSC_VER
#pragma warning(disable:4552)
#pragma warning(disable:4553)
#endif
#define ASSERT(x) do{if(!(x)){exit(EXIT_FAILURE);}}while(0);
#define ASSERT_MSG(x,msg) do{if(!(x)){WARNING(msg);exit(EXIT_FAILURE);}}while(0);
#define ASSERT_MSGV(x,fmt,...) do{if(!(x)){WARNINGV(fmt,__VA_ARGS__);exit(EXIT_FAILURE);}}while(0);

#endif

#define WARNING(msg) do{coutPrintf("[WARNING] \x1B[31m %s \x1B[0m\n",msg);fflush(stdout);}while(0);
#define WARNINGV(fmt,...) do{coutPrintf("[WARNING] \x1B[31m " fmt " \x1B[0m\n",__VA_ARGS__);fflush(stdout);}while(0);
#define INFO(msg) do{coutPrintf("[INFO] %s \n",msg);fflush(stdout);}while(0);
#define INFOV(fmt,...) do{coutPrintf("[INFO] " fmt " \n",__VA_ARGS__);fflush(stdout);}while(0);
#define NOTIFY_MSG(msg) do{coutPrintf("[NOTIFY] %s \n",msg);fflush(stdout);}while(0); getchar();
#define NOTIFY_MSGV(fmt,...) do{coutPrintf("[NOTIFY] " fmt " \n",__VA_ARGS__);fflush(stdout);}while(0); getchar();

//OpenMP only support signed variable as index
#include <stdint.h>
#define USE_QUAD_SIZE
#ifdef USE_QUAD_SIZE
typedef int64_t sizeType;
#else
typedef int sizeType;
#endif

typedef double scalarD;
typedef float scalarF;

typedef int vtkSizeType;
#define INVALID sizeType(-1)

#define FUNCTION_NOT_IMPLEMENTED ASSERT_MSGV(false,"Function \"%s\" not implemented!",__FUNCTION__)

#endif
