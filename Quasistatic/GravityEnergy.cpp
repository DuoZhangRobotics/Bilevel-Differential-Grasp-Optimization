#include "GravityEnergy.h"
#include "GraspPlanner.h"
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
GravityEnergy<T>::GravityEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const Vec3T& g,T coef)
  :ArticulatedObjective<T>(obj,"GravityEnergy(coef="+std::to_string(coef)+",g[0]="+std::to_string(g[0])+",g[1]="+std::to_string(g[1])+",g[2]="+std::to_string(g[2])+")",info,planner,object),_g(g),_coef(coef) {}
template <typename T>
int GravityEnergy<T>::operator()(const Vec&,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>*,Vec*,STrips*)
{
  for(sizeType i=0; i<_planner.body().nrJ(); i++) {
    const Joint& joint=_planner.body().joint(i);
    if(joint._M<=0)
      continue;
    Vec3T c=joint.getC().template cast<T>();
    Vec3T cG=ROTI(_info._TM,i)*c+CTRI(_info._TM,i);
    e+=-cG.dot(_g)*joint._M*_coef;
    if(g) {
      ROTI(g->getMatrixI(),i)-=_g*c.transpose()*joint._M*_coef;
      CTRI(g->getMatrixI(),i)-=_g*joint._M*_coef;
    }
  }
  return 0;
}
//instance
PRJ_BEGIN
template class GravityEnergy<double>;
#ifdef ALL_TYPES
template class GravityEnergy<__float128>;
template class GravityEnergy<mpfr::mpreal>;
#endif
PRJ_END
