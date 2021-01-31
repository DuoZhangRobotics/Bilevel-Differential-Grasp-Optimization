#include "CentroidClosednessEnergy.h"
#include "GraspPlanner.h"
#include <CommonFile/geom/StaticGeom.h>
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
CentroidClosednessEnergy<T>::CentroidClosednessEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T coef)
  :ArticulatedObjective<T>(obj,"CentroidClosednessEnergy(coef="+std::to_string(coef)+")",info,planner,object),_coef(coef)
{
  _centroid.resize(planner.body().nrJ());
  for(sizeType i=0; i<(sizeType)_centroid.size(); i++) {
    ObjMesh m;
    planner.body().getGeom().getG(i).getMesh(m);
    _centroid[i]=m.getVolumeCentroid().template cast<T>();
  }
}
template <typename T>
int CentroidClosednessEnergy<T>::operator()(const Vec&,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec*,STrips*)
{
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType i=0; i<(sizeType)_centroid.size(); i++)
      addTerm(i,e,g,h);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_centroid.size(); i++)
      addTerm(i,e,g,h);
  }
  return 0;
}
template <typename T>
void CentroidClosednessEnergy<T>::addTerm(sizeType i,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const
{
  if(std::isfinite(_centroid[i][0])) {
    Vec3T c=ROTI(_info._TM,i)*_centroid[i]+CTRI(_info._TM,i);
    T E=c.squaredNorm()*_coef/2;
    e+=E;
    if(g) {
      ROTI(g->getMatrixI(),i)+=c*_centroid[i].transpose()*_coef;
      CTRI(g->getMatrixI(),i)+=c*_coef;
    }
    if(h) {
      Vec4T cH(_centroid[i][0],_centroid[i][1],_centroid[i][2],1);
      Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->getMatrixI().coeffRef(0,i*12)));
      for(sizeType r=0; r<4; r++)
        for(sizeType c=0; c<4; c++)
          hBlk.template block<3,3>(r*3,c*3).diagonal().array()+=cH[r]*cH[c]*_coef;
    }
  }
}
//instance
PRJ_BEGIN
template class CentroidClosednessEnergy<double>;
#ifdef ALL_TYPES
template class CentroidClosednessEnergy<__float128>;
template class CentroidClosednessEnergy<mpfr::mpreal>;
#endif
PRJ_END
