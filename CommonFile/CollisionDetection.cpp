#include "CollisionDetectionImpl.h"

PRJ_BEGIN

#define DEF_BBOX(T,DIM) \
template struct BBox<T,DIM>; \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<sizeType>(const BBox<sizeType,DIM>& other);  \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<unsigned char>(const BBox<unsigned char,DIM>& other);  \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<scalarF>(const BBox<scalarF,DIM>& other);  \
template typename BBox<T,DIM>::BBox<T,DIM>& BBox<T,DIM>::copy<scalarD>(const BBox<scalarD,DIM>& other);

//this is a workaround for fast gauss transform used in GraspPlanner
#ifdef ALL_TYPES
DEF_BBOX(__float128,2)
DEF_BBOX(mpfr::mpreal,3)
template class Sphere<__float128>;
template class Sphere<mpfr::mpreal>;
#endif

DEF_BBOX(scalarD,2)
DEF_BBOX(scalarD,3)
template typename BBox<scalarF,2>::BBox<scalarF,2>& BBox<scalarF,2>::operator=<scalarD>(const BBox<scalarD,2>& other);
template typename BBox<scalarF,3>::BBox<scalarF,3>& BBox<scalarF,3>::operator=<scalarD>(const BBox<scalarD,3>& other);
template void TriangleTpl<scalarD>::calcPointDist(const Vec3d& pt,scalarD& sqrDistance,Vec3d& cp,Vec3d& b) const;
template void TriangleTpl<scalarD>::calcPointDist(const Vec3d& pt,scalarD& sqrDistance,Vec3d& cp,Vec4d& b) const;
template class LineSegTpl<scalarD>;
template class PlaneTpl<scalarD>;
template class TriangleTpl<scalarD>;
template class TetrahedronTpl<scalarD>;
template class OBBTpl<scalarD,2>;
template class OBBTpl<scalarD,3>;
template class KDOP18<scalarD>;
template class Sphere<scalarD>;

DEF_BBOX(scalarF,2)
DEF_BBOX(scalarF,3)
template typename BBox<scalarD,2>::BBox<scalarD,2>& BBox<scalarD,2>::operator=<scalarF>(const BBox<scalarF,2>& other);
template typename BBox<scalarD,3>::BBox<scalarD,3>& BBox<scalarD,3>::operator=<scalarF>(const BBox<scalarF,3>& other);
template void TriangleTpl<scalarF>::calcPointDist(const Vec3f& pt,scalarF& sqrDistance,Vec3f& cp,Vec3f& b) const;
template void TriangleTpl<scalarF>::calcPointDist(const Vec3f& pt,scalarF& sqrDistance,Vec3f& cp,Vec4f& b) const;
template class LineSegTpl<scalarF>;
template class PlaneTpl<scalarF>;
template class TriangleTpl<scalarF>;
template class TetrahedronTpl<scalarF>;
template class OBBTpl<scalarF,2>;
template class OBBTpl<scalarF,3>;
template class KDOP18<scalarF>;
template class Sphere<scalarF>;

DEF_BBOX(sizeType,2)
DEF_BBOX(sizeType,3)
template void TriangleTpl<sizeType>::calcPointDist(const Vec3i& pt,sizeType& sqrDistance,Vec3i& cp,Vec3i& b) const;
template void TriangleTpl<sizeType>::calcPointDist(const Vec3i& pt,sizeType& sqrDistance,Vec3i& cp,Vec4i& b) const;
template class LineSegTpl<sizeType>;
template class PlaneTpl<sizeType>;
template class TriangleTpl<sizeType>;
template class TetrahedronTpl<sizeType>;
template class OBBTpl<sizeType,2>;
template class OBBTpl<sizeType,3>;
template class KDOP18<sizeType>;
template class Sphere<sizeType>;

DEF_BBOX(char,2)
DEF_BBOX(char,3)
template void TriangleTpl<char>::calcPointDist(const Vec3c& pt,char& sqrDistance,Vec3c& cp,Vec3c& b) const;
template void TriangleTpl<char>::calcPointDist(const Vec3c& pt,char& sqrDistance,Vec3c& cp,Vec4c& b) const;
template class LineSegTpl<char>;
template class PlaneTpl<char>;
template class TriangleTpl<char>;
template class TetrahedronTpl<char>;
template class OBBTpl<char,2>;
template class OBBTpl<char,3>;
template class KDOP18<char>;
template class Sphere<char>;

DEF_BBOX(unsigned char,2)
DEF_BBOX(unsigned char,3)
template void TriangleTpl<unsigned char>::calcPointDist(const Vec3uc& pt,unsigned char& sqrDistance,Vec3uc& cp,Vec3uc& b) const;
template void TriangleTpl<unsigned char>::calcPointDist(const Vec3uc& pt,unsigned char& sqrDistance,Vec3uc& cp,Vec4uc& b) const;
template class LineSegTpl<unsigned char>;
template class PlaneTpl<unsigned char>;
template class TriangleTpl<unsigned char>;
template class TetrahedronTpl<unsigned char>;
template class OBBTpl<unsigned char,2>;
template class OBBTpl<unsigned char,3>;
template class KDOP18<unsigned char>;
template class Sphere<unsigned char>;

PRJ_END
