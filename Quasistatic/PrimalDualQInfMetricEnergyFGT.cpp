#include "PrimalDualQInfMetricEnergyFGT.h"
#include "GraspPlanner.h"
#include "FGTTreeNode.h"
#include <Utils/SparseUtils.h>
#include <chrono>
USE_PRJ_NAMESPACE

template <typename T>
PrimalDualQInfMetricEnergyFGT<T>::PrimalDualQInfMetricEnergyFGT(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const T& alpha,T coef,T normalExtrude,T FGTThres)
  :PrimalDualQInfMetricEnergy<T>(obj,info,planner,object,alpha,coef,SQR_EXP_ACTIVATION,normalExtrude),_FGTThres(FGTThres)
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
  DSSQPObjectiveComponent<T>::_name="PrimalDualQInfMetricEnergyFGT(alpha="+std::to_string(_alpha)+",coef="+std::to_string(_coef)+",FGTThres="+std::to_string(_FGTThres)+")";
}
template <typename T>
int PrimalDualQInfMetricEnergyFGT<T>::operator()(const Vec& x,Vec& fvec,STrips* fjac)
{
  T invHSqr=1/_alpha;
  Vec G=Vec::Zero(_pss.cols());
  std::vector<MatX4T,Eigen::aligned_allocator<MatX4T>> DGDT;
  if(fjac)
    DGDT.assign(_planner.body().nrJ(),MatX4T::Zero(3*_pss.cols(),4));
  for(sizeType i=0; i<_planner.body().nrJ(); i++) {
    const Mat3XT& yl=const_cast<Mat3XT&>(_planner.pnss()[i].first);
    if(yl.cols()==0)
      continue;
    Mat3XT y=ROTI(_info._TM,i)*yl+CTRI(_info._TM,i)*Vec::Ones(yl.cols()).transpose();
    _gripperFGT[i]->transform(ROTI(_info._TM,i),CTRI(_info._TM,i));
    FGTTreeNode<T>::FGT(G,fjac?&DGDT[i]:NULL,NULL,y,&yl,_pss,*_gripperFGT[i],*_objectFGT,invHSqr,_FGTThres);
  }

  sizeType nrC=values();
  T area=_planner.area();
  if(fjac) {
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
          TRANSI(DGDTc,j)+=_object.gij()(i,r)*DGDT[j].template block<3,4>(i*3,0);
      }
      DGDTc*=area;

      cjacRow=Vec::Zero(_planner.body().nrDOF());
      fjac->push_back(STrip(r+DSSQPObjectiveComponent<T>::_offset,MetricEnergy<T>::_off,-1));
      _info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(DGDTc),ArticulatedObjective<T>::mapV(cjacRow));
      addBlock(*fjac,r+DSSQPObjectiveComponent<T>::_offset,0,cjacRow.transpose());
    }
  }
  fvec.template segment(DSSQPObjectiveComponent<T>::_offset,nrC)=_object.gij().transpose()*G*area-Vec::Constant(nrC,x[MetricEnergy<T>::_off]);
  return 0;
}
//instance
PRJ_BEGIN
template class PrimalDualQInfMetricEnergyFGT<double>;
#ifdef ALL_TYPES
template class PrimalDualQInfMetricEnergyFGT<__float128>;
template class PrimalDualQInfMetricEnergyFGT<mpfr::mpreal>;
#endif
PRJ_END
