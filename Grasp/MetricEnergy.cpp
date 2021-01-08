#include "MetricEnergy.h"
#include "GraspPlanner.h"
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
MetricEnergy<T>::MetricEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T d0,const T& alpha,T coef,METRIC_TYPE m,METRIC_ACTIVATION a,T normalExtrude)
  :ArticulatedObjective<T>(planner,obj),_d0(d0),_coef(coef),_alpha(alpha),_type(m),_activation(a),_pss(obj.pss(normalExtrude)) {}
template <typename T>
int MetricEnergy<T>::operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* gFinal,MatT*)
{
  T area=_planner.area();
  Vec w=Vec::Zero(_obj.gij().rows());
  //point on object
  for(sizeType k=0; k<_pss.cols(); k++) {
    const Vec3T po=_pss.col(k);
    //point on hand
    for(sizeType i=0; i<_planner.body().nrJ(); i++) {
      //point on link i
      const std::pair<Mat3XT,Mat3XT>& pn=_planner.pnss()[i];
      Mat3XT p=ROTI(info._TM,i)*pn.first+CTRI(info._TM,i)*Vec::Ones(pn.first.cols()).transpose();
      for(sizeType j=0; j<p.cols(); j++) {
        T len=std::sqrt((p.col(j)-po).squaredNorm());
        w[k]+=activation(len)*area;
      }
    }
  }
  //Q-Metric
  Vec gw;
  if(_type==Q_1)
    e+=_obj.computeQ1(w,g?&gw:NULL)*_coef;
  else if(_type==Q_INF)
    e+=_obj.computeQInf(w,g?&gw:NULL)*_coef;
  else if(_type==Q_INF_BARRIER)
    e+=_obj.computeQInfBarrier(w,x[off],_d0,g?&gw:NULL)*_coef;
  else {
    ASSERT_MSGV(false,"Unknown metric type: %d",_type)
  }
  if(!std::isfinite(e.getValueI()))
    return -1;
  //assemble gradient/hessian
  if(g || h) {
    gw*=_coef;
    if(std::is_same<T,mpfr::mpreal>::value) {
      for(sizeType oid=0; oid<_pss.cols(); oid++)
        for(sizeType linkId=0; linkId<_planner.body().nrJ(); linkId++)
          for(sizeType linkPId=0; linkPId<_planner.pnss()[linkId].first.cols(); linkPId++)
            addTerm(area,gw,linkId,linkPId,oid,info,g,h);
    } else {
      OMP_PARALLEL_FOR_
      for(sizeType oid=0; oid<_pss.cols(); oid++)
        for(sizeType linkId=0; linkId<_planner.body().nrJ(); linkId++)
          for(sizeType linkPId=0; linkPId<_planner.pnss()[linkId].first.cols(); linkPId++)
            addTerm(area,gw,linkId,linkPId,oid,info,g,h);
    }
    if(_type==Q_INF_BARRIER)
      gFinal->coeffRef(off)+=gw[_obj.gij().rows()];
  }
  return 0;
}
template <typename T>
sizeType MetricEnergy<T>::nrAdditionalDOF() const
{
  return _type==Q_INF_BARRIER?1:0;
}
template <typename T>
T MetricEnergy<T>::activation(T param,T* D,T* DD) const
{
  T ret,d,coef;
  switch(_activation)
  {
  case EXP_ACTIVATION:
    ret=std::exp(-param/_alpha);
    if(D)
      *D=-ret/_alpha;
    if(DD)
      *DD=ret/(_alpha*_alpha);
    break;
  case SQR_EXP_ACTIVATION:
    ret=std::exp(-param*param/_alpha);
    if(D)
      *D=-ret*2*param/_alpha;
    if(DD)
      *DD=ret*(4*param*param)/(_alpha*_alpha)-ret*2/_alpha;
    break;
  case INVERSE_ACTIVATION:
    d=1/(param+_alpha);
    ret=_alpha*d;
    if(D)
      *D=-_alpha*d*d;
    if(DD)
      *DD=2*_alpha*d*d*d;
    break;
  case SQR_INVERSE_ACTIVATION:
    coef=-1/(2*M_PI*std::log(_alpha));
    d=1/(param*param+_alpha*_alpha);
    ret=coef*d;
    if(D)
      *D=-2*coef*d*d*param;
    if(DD)
      *DD=2*coef*d*d*(d*(4*param*param)-1);
    break;
  default:
    ASSERT_MSG(false,"Unknown activation")
  }
  return ret;
}
template <typename T>
std::string MetricEnergy<T>::name() const
{
  return "MetricEnergy(d0="+std::to_string(_d0)+",alpha="+std::to_string(_alpha)+",coef="+std::to_string(_coef)+",type="+std::to_string(_type)+")";
}
template <typename T>
void MetricEnergy<T>::addTerm(T area,const Vec& gw,sizeType linkId,sizeType linkPId,sizeType oid,const PBDArticulatedGradientInfo<T>& info,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const
{
  if(gw[oid]==0)
    return;
  Vec3T po=_pss.col(oid);
  Vec3T pg=_planner.pnss()[linkId].first.col(linkPId);
  Vec3T pgG=ROTI(info._TM,linkId)*pg+CTRI(info._TM,linkId);
  T len=std::sqrt((pgG-po).squaredNorm()),D,DD;
  activation(len,&D,&DD);
  Vec3T dpUnit=(pgG-po)/len;
  T coef=D*area*gw[oid];
  Vec3T G=dpUnit*coef;
  if(g) {
    CTRI(g->getMatrixI(),linkId)+=G;
    ROTI(g->getMatrixI(),linkId)+=G*pg.transpose();
  }
  if(h) {
    Mat3T H=dpUnit*dpUnit.transpose()*(DD*area*gw[oid]-coef/len)+Mat3T::Identity()*(coef/len);
    Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->getMatrixI().coeffRef(0,linkId*12)));
    Vec4T cH(pg[0],pg[1],pg[2],1);
    for(sizeType r=0; r<4; r++)
      for(sizeType c=0; c<4; c++)
        hBlk.template block<3,3>(r*3,c*3)+=cH[r]*cH[c]*H;
  }
}
//instance
PRJ_BEGIN
template class MetricEnergy<double>;
#ifdef ALL_TYPES
template class MetricEnergy<__float128>;
template class MetricEnergy<mpfr::mpreal>;
#endif
PRJ_END
