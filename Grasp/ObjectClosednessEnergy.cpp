#include "ObjectClosednessEnergy.h"
#include <CommonFile/geom/StaticGeom.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
ObjectClosednessEnergy<T>::ObjectClosednessEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T coef):ArticulatedObjective<T>(planner,obj),_coef(_planner.area()*coef) {}
template <typename T>
int ObjectClosednessEnergy<T>::operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h)
{
  const ObjMeshGeomCellExact& distCalc=_obj.dist();
  for(sizeType i=0; i<_planner.body().nrJ(); i++) {
    const Mat3XT& pss=_planner.pnss()[i].first;
    for(sizeType j=0; j<pss.cols(); j++) {
      Vec3T pL=ROTI(info._TM,i)*pss.col(j)+CTRI(info._TM,i);
      Vec3T n,normal;
      Mat3T hessian;
      Vec2i feat;
      T dist=distCalc.template closest<T>(pL,n,normal,hessian,feat);
      e+=dist*dist*_coef/2;
      if(g) {
        ROTI(*g,i)+=normal*pss.col(j).transpose()*dist*_coef;
        CTRI(*g,i)+=normal*dist*_coef;
      }
      if(h) {
        hessian=(hessian*dist+normal*normal.transpose())*_coef;
        Vec4T cH(pss.col(j)[0],pss.col(j)[1],pss.col(j)[2],1);
        Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->coeffRef(0,i*12)));
        for(sizeType r=0; r<4; r++)
          for(sizeType c=0; c<4; c++)
            hBlk.template block<3,3>(r*3,c*3)+=hessian*cH[r]*cH[c];
      }
    }
  }
  return 0;
}
//instance
PRJ_BEGIN
template class ObjectClosednessEnergy<double>;
#ifdef ALL_TYPES
template class ObjectClosednessEnergy<__float128>;
template class ObjectClosednessEnergy<mpfr::mpreal>;
#endif
PRJ_END
