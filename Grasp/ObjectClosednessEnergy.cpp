#include "ObjectClosednessEnergy.h"
#include "GraspPlanner.h"
#include <CommonFile/geom/StaticGeom.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
ObjectClosednessEnergy<T>::ObjectClosednessEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T coef):ArticulatedObjective<T>(planner,obj),_coef(_planner.area()*coef) {}
template <typename T>
int ObjectClosednessEnergy<T>::operator()(const Vec&,const PBDArticulatedGradientInfo<T>& info,sizeType,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec*,MatT*)
{
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> terms;
  for(sizeType i=0; i<_planner.body().nrJ(); i++)
    for(sizeType j=0; j<_planner.pnss()[i].first.cols(); j++)
      terms.push_back(Vec2i(i,j));

  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType termId=0; termId<(sizeType)terms.size(); termId++)
      addTerm(terms[termId],info,e,g,h);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType termId=0; termId<(sizeType)terms.size(); termId++)
      addTerm(terms[termId],info,e,g,h);
  }
  return 0;
}
template <typename T>
std::string ObjectClosednessEnergy<T>::name() const
{
  return "ObjectClosednessEnergy(coef="+std::to_string(_coef)+")";
}
template <typename T>
void ObjectClosednessEnergy<T>::addTerm(const Vec2i& termId,const PBDArticulatedGradientInfo<T>& info,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const
{
  sizeType i=termId[0];
  sizeType j=termId[1];
  const Mat3XT& pss=_planner.pnss()[i].first;
  const ObjMeshGeomCellExact& distCalc=_obj.dist();
  Vec3T pL=ROTI(info._TM,i)*pss.col(j)+CTRI(info._TM,i);
  Vec3T n,normal;
  Mat3T hessian;
  Vec2i feat;
  T dist=distCalc.template closest<T>(pL,n,normal,hessian,feat);
  T E=dist*dist*_coef/2;
  e+=E;
  if(g) {
    ROTI(g->getMatrixI(),i)+=normal*pss.col(j).transpose()*dist*_coef;
    CTRI(g->getMatrixI(),i)+=normal*dist*_coef;
  }
  if(h) {
    hessian=(hessian*dist+normal*normal.transpose())*_coef;
    Vec4T cH(pss.col(j)[0],pss.col(j)[1],pss.col(j)[2],1);
    Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->getMatrixI().coeffRef(0,i*12)));
    for(sizeType r=0; r<4; r++)
      for(sizeType c=0; c<4; c++)
        hBlk.template block<3,3>(r*3,c*3)+=hessian*cH[r]*cH[c];
  }
}
//instance
PRJ_BEGIN
template class ObjectClosednessEnergy<double>;
#ifdef ALL_TYPES
template class ObjectClosednessEnergy<__float128>;
template class ObjectClosednessEnergy<mpfr::mpreal>;
#endif
PRJ_END
