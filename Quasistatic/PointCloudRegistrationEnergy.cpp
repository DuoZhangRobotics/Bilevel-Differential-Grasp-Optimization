#include "PointCloudRegistrationEnergy.h"
#include "GraspPlanner.h"
#include <Environment/Environment.h>

USE_PRJ_NAMESPACE

template <typename T>
bool PointCloudRegistrationEnergy<T>::PointCloudMap::operator<(const PointCloudMap& other) const {
  return _pid<other._pid;
}
template <typename T>
bool PointCloudRegistrationEnergy<T>::PointCloudMap::operator!=(const PointCloudMap& other) const {
  return *this<other || other<*this;
}
template <typename T>
bool PointCloudRegistrationEnergy<T>::PointCloudMap::operator==(const PointCloudMap& other) const {
  return !(*this<other) && !(other<*this);
}
template <typename T>
PointCloudRegistrationEnergy<T>::PointCloudRegistrationEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T coef)
  :ArticulatedObjective<T>(obj,"PointCloudRegistrationEnergy(coef="+std::to_string(coef)+")",info,planner,object),_coef(coef)
{
  for(sizeType i=0; i<object.idss().size(); i++)
    if(object.idss()[i]>=0)
      addPointCloudMap(i,object.idss()[i],Vec3T::Constant(ScalarUtil<T>::scalar_nanq()));
}
template <typename T>
void PointCloudRegistrationEnergy<T>::addPointCloudMap(sizeType pid,sizeType oid,const Vec3T& objPos)
{
  PointCloudMap pm;
  pm._pid=pid;
  pm._oid=oid;
  pm._objPos=objPos;

  ASSERT_MSGV(pm._pid>=0 && pm._pid<_object.pss().cols(),"Invalid point id: 0<=%d<%d",pm._pid,_object.pss().cols())
  ASSERT_MSGV(pm._oid>=0 && pm._oid<_planner.body().nrDOF()/6,"Invalid object id: 0<=%d<%d",pm._oid,_planner.body().nrDOF()/6)
  pm._oid=pm._oid*2+2;  //remap to linkid

  typename std::vector<PointCloudMap>::iterator it=std::lower_bound(_pmss.begin(),_pmss.end(),pm);
  if(it==_pmss.end())
    _pmss.insert(it,pm);
  else if(*it!=pm)
    _pmss.insert(it,pm);
  else *it=pm;
}
template <typename T>
bool PointCloudRegistrationEnergy<T>::existPointCloudMap(sizeType pid,sizeType oid) const
{
  PointCloudMap pm;
  pm._pid=pid;
  pm._oid=oid;

  ASSERT_MSGV(pm._pid>=0 && pm._pid<_object.pss().cols(),"Invalid point id: 0<=%d<%d",pm._pid,_object.pss().cols())
  ASSERT_MSGV(pm._oid>=0 && pm._oid<_planner.body().nrDOF()/6,"Invalid object id: 0<=%d<%d",pm._oid,_planner.body().nrDOF()/6)
  pm._oid=pm._oid*2+2;  //remap to linkid

  typename std::vector<PointCloudMap>::const_iterator it=std::lower_bound(_pmss.begin(),_pmss.end(),pm);
  return it!=_pmss.end() && *it==pm;
}
template <typename T>
void PointCloudRegistrationEnergy<T>::clearPointCloudMap()
{
  _pmss.clear();
}
//objective
template <typename T>
int PointCloudRegistrationEnergy<T>::operator()(const Vec&,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec*,STrips*)
{
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType i=0; i<(sizeType)_pmss.size(); i++)
      addTerm(_pmss[i],e,g,h);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)_pmss.size(); i++)
      addTerm(_pmss[i],e,g,h);
  }
  return 0;
}
template <typename T>
void PointCloudRegistrationEnergy<T>::addTerm(const PointCloudMap& pm,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const
{
  Mat3T hessian;
  Mat3X4T tLink=TRANSI(_info._TM,pm._oid);
  Vec3T p=_object.pss().col(pm._pid),normal;
  Vec3T pL=ROTI(_info._TM,pm._oid).transpose()*(p-CTR(tLink));
  T D=_coef,DD=_coef;

  //mode selection
  if(std::isfinite(pm._objPos[0])) {
    //this indicates mapping to a known local point
    normal=pL-pm._objPos;
    e+=normal.squaredNorm()*_coef/2;
    if(h)
      hessian=Mat3T::Identity()*DD;
  } else {
    //this indicates finding the closest point on the mesh
    const Environment<T>& env=_planner.env(pm._oid);
    T phi=env.phi(pL,g?&normal:NULL);
    e+=phi*phi*_coef/2;
    D*=phi;
    if(h) {
      env.phiGrad(pL,&hessian);
      hessian=(hessian*D+normal*normal.transpose()*DD).eval();
    }
  }

  //gradient calculation
  if(g) {
    p-=CTR(tLink);
    CTRI(g->getMatrixI(),pm._oid)+=-D*ROT(tLink)*normal;
    ROTI(g->getMatrixI(),pm._oid)+=D*p*normal.transpose();
  }
  if(h) {
    Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->getMatrixI().coeffRef(0,pm._oid*12)));
    Mat3T HRT=-hessian*ROT(tLink).transpose();
    for(sizeType r=0; r<3; r++) {
      hBlk.template block<3,3>(r*3,9)+=p*HRT.row(r)-Mat3T::Identity()*normal[r]*D;
      hBlk.template block<3,3>(9,r*3)+=(p*HRT.row(r)-Mat3T::Identity()*normal[r]*D).transpose();
      for(sizeType c=0; c<3; c++)
        hBlk.template block<3,3>(r*3,c*3)+=p*p.transpose()*hessian(r,c);
    }
    hBlk.template block<3,3>(9,9)-=ROT(tLink)*HRT;
  }
}
//instance
PRJ_BEGIN
template class PointCloudRegistrationEnergy<double>;
#ifdef ALL_TYPES
template class PointCloudRegistrationEnergy<__float128>;
template class PointCloudRegistrationEnergy<mpfr::mpreal>;
#endif
PRJ_END
