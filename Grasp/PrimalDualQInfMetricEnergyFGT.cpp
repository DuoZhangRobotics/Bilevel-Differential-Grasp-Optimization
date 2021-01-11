#include "PrimalDualQInfMetricEnergyFGT.h"
#include "GraspPlanner.h"
#include "FGTTreeNode.h"

USE_PRJ_NAMESPACE

template <typename T>
PrimalDualQInfMetricEnergyFGT<T>::PrimalDualQInfMetricEnergyFGT(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,const T& alpha,T coef,T normalExtrude,T FGTThres)
  :PrimalDualQInfMetricEnergy<T>(planner,obj,alpha,coef,SQR_EXP_ACTIVATION,normalExtrude),_FGTThres(FGTThres)
{
#define LEAF_THRES_FGT 32
  _objectFGT.reset(new FGTTreeNode<T>(NULL,_pss,Vec2i(0,_pss.cols()),LEAF_THRES_FGT));
  _gripperFGT.resize(_planner.body().nrJ());
  for(sizeType i=0; i<_planner.body().nrJ(); i++) {
    Mat3XT& yl=const_cast<Mat3XT&>(_planner.pnss()[i].first);
    Mat3XT& yln=const_cast<Mat3XT&>(_planner.pnss()[i].second);
    if(yl.cols()==0)
      continue;
    _gripperFGT[i].reset(new FGTTreeNode<T>(NULL,yl,Vec2i(0,yl.cols()),LEAF_THRES_FGT,&yln));
  }
#undef LEAF_THRES_FGT
}
template <typename T>
void PrimalDualQInfMetricEnergyFGT<T>::cons(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType offr,sizeType offc,Vec& c,MatT* cjac)
{
  T invHSqr=1/_alpha;
  Vec G=Vec::Zero(_pss.cols());
  std::vector<MatX4T,Eigen::aligned_allocator<MatX4T>> DGDT;
  if(cjac)
    DGDT.assign(_planner.body().nrJ(),MatX4T::Zero(3*_pss.cols(),4));
  for(sizeType i=0; i<_planner.body().nrJ(); i++) {
    const Mat3XT& yl=const_cast<Mat3XT&>(_planner.pnss()[i].first);
    if(yl.cols()==0)
      continue;
    Mat3XT y=ROTI(info._TM,i)*yl+CTRI(info._TM,i)*Vec::Ones(yl.cols()).transpose();
    _gripperFGT[i]->transform(ROTI(info._TM,i),CTRI(info._TM,i));
    FGTTreeNode<T>::FGT(G,cjac?&DGDT[i]:NULL,NULL,y,&yl,_pss,*_gripperFGT[i],*_objectFGT,invHSqr,_FGTThres);
  }
  sizeType nrC=nrCons();
  T area=_planner.area();
  if(cjac) {
    Vec cjacRow;
    Mat3XT DGDTc;
    DGDTc.resize(3,_planner.body().nrJ()*4);
    for(sizeType r=0; r<nrC; r++) {
      DGDTc.setZero();
      OMP_PARALLEL_FOR_
      for(sizeType j=0; j<_planner.body().nrJ(); j++) {
        if(!_gripperFGT[j])
          continue;
        for(sizeType i=0; i<_pss.cols(); i++)
          TRANSI(DGDTc,j)+=_obj.gij()(i,r)*DGDT[j].template block<3,4>(i*3,0);
      }
      DGDTc*=area;
      cjacRow=-Vec::Unit(cjac->cols(),offc);
      info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(DGDTc),ArticulatedObjective<T>::mapV(cjacRow));
      cjac->row(r+offr)=cjacRow;
    }
  }
  c.template segment(offr,nrC)=_obj.gij().transpose()*G*area-Vec::Constant(nrC,x[offc]);
}
template <typename T>
std::string PrimalDualQInfMetricEnergyFGT<T>::name() const
{
  return "PrimalDualQInfMetricEnergyFGT(alpha="+std::to_string(_alpha)+",coef="+std::to_string(_coef)+",FGTThres="+std::to_string(_FGTThres)+")";
}
//instance
PRJ_BEGIN
template class PrimalDualQInfMetricEnergyFGT<double>;
#ifdef ALL_TYPES
template class PrimalDualQInfMetricEnergyFGT<__float128>;
template class PrimalDualQInfMetricEnergyFGT<mpfr::mpreal>;
#endif
PRJ_END
