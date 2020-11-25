#include "ObjMesh.h"
#include "IO.h"
#include "DisjointSet.h"
#include "CollisionDetection.h"
#include <stack>
#include <set>

USE_PRJ_NAMESPACE

//helper
template<typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T getAngle3D(const typename ScalarUtil<T>::ScalarVec3& a,const typename ScalarUtil<T>::ScalarVec3& b)
{
  T denom=std::max<T>(a.norm()*b.norm(),ScalarUtil<T>::scalar_eps());
  return acos(std::max<T>(std::min<T>(a.dot(b)/denom,(T)(1.0f-ScalarUtil<T>::scalar_eps())),(T)(-1.0f+ScalarUtil<T>::scalar_eps())));
}
template<typename T>
FORCE_INLINE EIGEN_DEVICE_FUNC T getAngle2D(const typename ScalarUtil<T>::ScalarVec2& a,const typename ScalarUtil<T>::ScalarVec2& b)
{
  T denom=std::max<T>(a.norm()*b.norm(),ScalarUtil<T>::scalar_eps());
  return acos(std::max<T>(std::min<T>(a.dot(b)/denom,(T)(1.0f-ScalarUtil<T>::scalar_eps())),(T)(-1.0f+ScalarUtil<T>::scalar_eps())));
}
template <typename T>
class WriteObjVertex
{
public:
  static void write(char* buf,T& a,T& b,T& c) {
#ifdef _MSC_VER
    sscanf_s(buf,"%f %f %f",&a,&b,&c);
#else
    sscanf(buf,"%f %f %f",&a,&b,&c);
#endif
  }
};
template <>
class WriteObjVertex<scalarD>
{
public:
  static void write(char* buf,scalarD& a,scalarD& b,scalarD& c) {
    double ta,tb,tc;
#ifdef _MSC_VER
    sscanf_s(buf,"%lf %lf %lf",&ta,&tb,&tc);
#else
    sscanf(buf,"%lf %lf %lf",&ta,&tb,&tc);
#endif
    a=std::convert<double>()(ta);
    b=std::convert<double>()(tb);
    c=std::convert<double>()(tc);
  }
};
template <typename TC>
void findInsert(std::map<int,int>& Ring,const TC& tris,int currR)
{
  for(typename TC::const_iterator beg=tris.begin(),end=tris.end(); beg!=end; beg++)
    if(Ring.find(*beg) == Ring.end())
      Ring.insert(std::pair<int,int>(*beg,currR));
}

#define SCALAR_NAME scalarF
#define OBJMESH ObjMeshF
#include "ObjMeshImpl.h"
#undef OBJMESH
#undef SCALAR_NAME

#define SCALAR_NAME scalarD
#define OBJMESH ObjMeshD
#include "ObjMeshImpl.h"
#undef OBJMESH
#undef SCALAR_NAME
