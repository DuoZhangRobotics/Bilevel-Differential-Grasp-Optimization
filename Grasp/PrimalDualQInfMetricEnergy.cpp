#include "PrimalDualQInfMetricEnergy.h"
#include "GraspPlanner.h"

USE_PRJ_NAMESPACE

template <typename T>
PrimalDualQInfMetricEnergy<T>::PrimalDualQInfMetricEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,const T& alpha,T coef,METRIC_ACTIVATION a,T normalExtrude)
  :MetricEnergy<T>(planner,obj,0,alpha,coef,Q_INF_CONSTRAINT,a,normalExtrude) {}
template <typename T>
int PrimalDualQInfMetricEnergy<T>::operator()(const Vec& x,const PBDArticulatedGradientInfo<T>&,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>*,ParallelMatrix<Mat12XT>*,Vec* gFinal,MatT*)
{
  e.getValueI()+=x[off]*_coef;
  if(gFinal)
    gFinal->coeffRef(off)+=_coef;
  return 0;
}
template <typename T>
void PrimalDualQInfMetricEnergy<T>::cons(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType offr,sizeType offc,Vec& c,MatT* cjac)
{
  sizeType nrC=nrCons();
  T area=_planner.area();
  ParallelMatrix<Vec> linkObjCoef(Vec::Zero(_pss.cols()));
  MatT linkObjCoefG;
  if(cjac)
    linkObjCoefG.setZero(_pss.cols()*3,_planner.body().nrJ()*4);
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType oid=0; oid<_pss.cols(); oid++)
      for(sizeType linkId=0; linkId<_planner.body().nrJ(); linkId++) {
        Eigen::Map<Mat3X4T,0,Eigen::OuterStride<>> linkObjCoefGM(cjac?&(linkObjCoefG.coeffRef(oid*3,linkId*4)):NULL,3,4,linkObjCoefG.outerStride());
        linkObjCoef.getMatrixI()[oid]+=addTerm(area,linkId,oid,info,linkObjCoefGM);
      }
    if(cjac) {
      for(sizeType i=0; i<nrC; i++) {
        Vec cjacRow=-Vec::Unit(cjac->cols(),offc);
        Mat3XT G=Mat3XT::Zero(3,_planner.body().nrJ()*4);
        for(sizeType oid=0; oid<_pss.cols(); oid++)
          G+=linkObjCoefG.block(3*oid,0,3,_planner.body().nrJ()*4)*_obj.gij()(oid,i);
        info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(G),ArticulatedObjective<T>::mapV(cjacRow));
        cjac->row(offr+i)=cjacRow;
      }
    }
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType oid=0; oid<_pss.cols(); oid++)
      for(sizeType linkId=0; linkId<_planner.body().nrJ(); linkId++) {
        Eigen::Map<Mat3X4T,0,Eigen::OuterStride<>> linkObjCoefGM(cjac?&(linkObjCoefG.coeffRef(oid*3,linkId*4)):NULL,3,4,linkObjCoefG.outerStride());
        linkObjCoef.getMatrixI()[oid]+=addTerm(area,linkId,oid,info,linkObjCoefGM);
      }
    if(cjac) {
      OMP_PARALLEL_FOR_
      for(sizeType i=0; i<nrC; i++) {
        Vec cjacRow=-Vec::Unit(cjac->cols(),offc);
        Mat3XT G=Mat3XT::Zero(3,_planner.body().nrJ()*4);
        for(sizeType oid=0; oid<_pss.cols(); oid++)
          G+=linkObjCoefG.block(3*oid,0,3,_planner.body().nrJ()*4)*_obj.gij()(oid,i);
        info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(G),ArticulatedObjective<T>::mapV(cjacRow));
        cjac->row(offr+i)=cjacRow;
      }
    }
  }
  c.template segment(offr,nrC)=_obj.gij().transpose()*linkObjCoef.getMatrix()-Vec::Constant(nrC,x[offc]);
}
template <typename T>
sizeType PrimalDualQInfMetricEnergy<T>::nrAdditionalDOF() const
{
  return 1;
}
template <typename T>
sizeType PrimalDualQInfMetricEnergy<T>::nrCons() const
{
  return _obj.gij().cols();
}
template <typename T>
std::string PrimalDualQInfMetricEnergy<T>::name() const
{
  return "PrimalDualQInfMetricEnergy(alpha="+std::to_string(_alpha)+",coef="+std::to_string(_coef)+")";
}
template <typename T>
T PrimalDualQInfMetricEnergy<T>::addTerm(T area,sizeType linkId,sizeType oid,const PBDArticulatedGradientInfo<T>& info,Mat3X4TM cjacG) const
{
  T ret=0;
  Vec3T po=_pss.col(oid);
  const Mat3XT& pg=_planner.pnss()[linkId].first;
  for(sizeType linkPId=0; linkPId<pg.cols(); linkPId++) {
    Vec3T pgG=ROTI(info._TM,linkId)*pg.col(linkPId)+CTRI(info._TM,linkId);
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
