#include "PrimalDualQInfMetricEnergy.h"
#include "GraspPlanner.h"

USE_PRJ_NAMESPACE

template <typename T>
PrimalDualQInfMetricEnergy<T>::PrimalDualQInfMetricEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const T& alpha,T coef,METRIC_ACTIVATION a,T normalExtrude)
  :MetricEnergy<T>(obj,info,planner,object,0,alpha,coef,Q_INF_CONSTRAINT,a,normalExtrude)
{
  DSSQPObjectiveComponent<T>::_name="PrimalDualQInfMetricEnergy(alpha="+std::to_string(_alpha)+",coef="+std::to_string(_coef)+")";
  for(sizeType i=0; i<values(); i++) {
    DSSQPObjectiveComponent<T>::_gl.push_back(0);
    DSSQPObjectiveComponent<T>::_gu.push_back(DSSQPObjective<T>::infty());
  }
}
template <typename T>
int PrimalDualQInfMetricEnergy<T>::operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>*,ParallelMatrix<Mat12XT>*,Vec* fgrad,STrips*)
{
  e.getValueI()+=x[MetricEnergy<T>::_off]*_coef;
  if(fgrad)
    fgrad->coeffRef(MetricEnergy<T>::_off)+=_coef;
  return 0;
}
template <typename T>
int PrimalDualQInfMetricEnergy<T>::operator()(const Vec& x,Vec& fvec,STrips* fjac)
{
  sizeType nrC=values();
  T area=_planner.area();
  ParallelMatrix<Vec> linkObjCoef(Vec::Zero(_pss.cols()));
  MatT linkObjCoefG;
  if(fjac)
    linkObjCoefG.setZero(_pss.cols()*3,_planner.body().nrJ()*4);
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType oid=0; oid<_pss.cols(); oid++)
      for(sizeType linkId=0; linkId<_planner.body().nrJ(); linkId++) {
        Eigen::Map<Mat3X4T,0,Eigen::OuterStride<>> linkObjCoefGM(fjac?&(linkObjCoefG.coeffRef(oid*3,linkId*4)):NULL,3,4,linkObjCoefG.outerStride());
        linkObjCoef.getMatrixI()[oid]+=addTerm(area,linkId,oid,linkObjCoefGM);
      }
    if(fjac) {
      for(sizeType i=0; i<nrC; i++) {
        Vec cjacRow=Vec::Zero(_planner.body().nrDOF());
        Mat3XT G=Mat3XT::Zero(3,_planner.body().nrJ()*4);
        for(sizeType oid=0; oid<_pss.cols(); oid++)
          G+=linkObjCoefG.block(3*oid,0,3,_planner.body().nrJ()*4)*_object.gij()(oid,i);
        fjac->push_back(STrip(i+DSSQPObjectiveComponent<T>::_offset,MetricEnergy<T>::_off,-1));
        _info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(G),ArticulatedObjective<T>::mapV(cjacRow));
        addBlock(*fjac,i+DSSQPObjectiveComponent<T>::_offset,0,cjacRow.transpose());
      }
    }
  } else {
    //OMP_PARALLEL_FOR_
    for(sizeType oid=0; oid<_pss.cols(); oid++)
      for(sizeType linkId=0; linkId<_planner.body().nrJ(); linkId++) {
        Eigen::Map<Mat3X4T,0,Eigen::OuterStride<>> linkObjCoefGM(fjac?&(linkObjCoefG.coeffRef(oid*3,linkId*4)):NULL,3,4,linkObjCoefG.outerStride());
        linkObjCoef.getMatrixI()[oid]+=addTerm(area,linkId,oid,linkObjCoefGM);
      }
    if(fjac) {
      //OMP_PARALLEL_FOR_
      for(sizeType i=0; i<nrC; i++) {
        Vec cjacRow=Vec::Zero(_planner.body().nrDOF());
        Mat3XT G=Mat3XT::Zero(3,_planner.body().nrJ()*4);
        for(sizeType oid=0; oid<_pss.cols(); oid++)
          G+=linkObjCoefG.block(3*oid,0,3,_planner.body().nrJ()*4)*_object.gij()(oid,i);
        fjac->push_back(STrip(i+DSSQPObjectiveComponent<T>::_offset,MetricEnergy<T>::_off,-1));
        _info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(G),ArticulatedObjective<T>::mapV(cjacRow));
        addBlock(*fjac,i+DSSQPObjectiveComponent<T>::_offset,0,cjacRow.transpose());
      }
    }
  }
  fvec.template segment(DSSQPObjectiveComponent<T>::_offset,nrC)=_object.gij().transpose()*linkObjCoef.getMatrix()-Vec::Constant(nrC,x[MetricEnergy<T>::_off]);
  return 0;
}
template <typename T>
int PrimalDualQInfMetricEnergy<T>::values() const
{
  return _object.gij().cols();
}
template <typename T>
T PrimalDualQInfMetricEnergy<T>::addTerm(T area,sizeType linkId,sizeType oid,Mat3X4TM cjacG) const
{
  T ret=0;
  Vec3T po=_pss.col(oid);
  const Mat3XT& pg=_planner.pnss()[linkId].first;
  for(sizeType linkPId=0; linkPId<pg.cols(); linkPId++) {
    Vec3T pgG=ROTI(_info._TM,linkId)*pg.col(linkPId)+CTRI(_info._TM,linkId);
    T len=std::sqrt((pgG-po).squaredNorm()),D;
    T A=MetricEnergy<T>::activation(len,&D);
    Vec3T dpUnit=(pgG-po)/len;
    Vec3T G=dpUnit*(D*area);
    if(cjacG.data())
      cjacG+=G*Vec4T(pg(0,linkPId),pg(1,linkPId),pg(2,linkPId),1).transpose();
    ret+=A*area;
  }
  return ret;
}
//instance
PRJ_BEGIN
template class PrimalDualQInfMetricEnergy<double>;
#ifdef ALL_TYPES
template class PrimalDualQInfMetricEnergy<__float128>;
template class PrimalDualQInfMetricEnergy<mpfr::mpreal>;
#endif
PRJ_END
