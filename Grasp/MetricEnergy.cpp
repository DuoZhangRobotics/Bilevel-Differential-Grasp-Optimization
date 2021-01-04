#include "MetricEnergy.h"
#include <stack>

USE_PRJ_NAMESPACE
/**
 * *Input:
 * @param planner: Grasp planner //Type: GraspPlanner
 * @param obj: Target object //Type: GraspQualityMetric 
 * @param alpha: Coefficient of distance, exp(-alpha * d) in Q metric //Type: T
 * @param coef: Coefficient of objective function //Type: T
 * *Intialize log barrier energy function
*/
template <typename T>
MetricEnergy<T>::MetricEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T alpha,METRIC_TYPE m,T coef)
  :ArticulatedObjective<T>(planner,obj),_alpha(alpha),_coef(coef),_type(m) {}
template <typename T>
int MetricEnergy<T>::operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h)
{
  T area=_planner.area();
  Vec w=Vec::Zero(_obj.gij().rows());
  //point on object
  const Mat3XT& poss=_obj.pss();
  for(sizeType k=0; k<poss.cols(); k++) {
    const Vec3T po=poss.col(k);
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
  if(_type==Q_INF_MEAN)
    e+=_obj.computeQInfMean(w,g?&gw:NULL)*_coef;
  else if(_type==Q_INF)
    e+=_obj.computeQInf(w,g?&gw:NULL)*_coef;
  else if(_type==Q_1)
    e+=_obj.computeQ1(w,g?&gw:NULL)*_coef;
  else {
    ASSERT_MSGV(false,"Unknown metric type: %d",_type)
  }
  //assemble gradient/hessian
  if(g || h) {
    gw*=_coef;
    for(sizeType k=0; k<poss.cols(); k++) {
      const Vec3T po=poss.col(k);
      if(gw[k]==0)
        continue;
      //point on hand
      for(sizeType i=0; i<_planner.body().nrJ(); i++) {
        //point on link i
        const std::pair<Mat3XT,Mat3XT>& pn=_planner.pnss()[i];
        Mat3XT p=ROTI(info._TM,i)*pn.first+CTRI(info._TM,i)*Vec::Ones(pn.first.cols()).transpose();
        for(sizeType j=0; j<p.cols(); j++) {
          T len=std::sqrt((p.col(j)-po).squaredNorm()),D,DD;
          activation(len,&D,&DD);
          T coef=D*area*gw[k];
          Vec3T dpUnit=(p.col(j)-po)/len;
          Vec3T G=dpUnit*coef;
          if(g) {
            CTRI(*g,i)+=G;
            ROTI(*g,i)+=G*pn.first.col(j).transpose();
          }
          if(h) {
            Mat3T H=dpUnit*dpUnit.transpose()*(DD*area*gw[k]-coef/len)+Mat3T::Identity()*(coef/len);
            Vec4T cH(pn.first.col(j)[0],pn.first.col(j)[1],pn.first.col(j)[2],1);
            Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->coeffRef(0,i*12)));
            for(sizeType r=0; r<4; r++)
              for(sizeType c=0; c<4; c++)
                hBlk.template block<3,3>(r*3,c*3)+=cH[r]*cH[c]*H;
          }
        }
      }
    }
  }
  return 0;
}
template <typename T>
T MetricEnergy<T>::activation(T param,T* D,T* DD) const
{
//#define EXP_ACTIVATION
#ifdef EXP_ACTIVATION
  T ret=std::exp(-_alpha*param);
  if(D)
    *D=-_alpha*ret;
  if(DD)
    *DD=_alpha*_alpha*ret;
  return ret;
#else
  param+=_planner.rad();
  T ret=1/_alpha/param;
  if(D)
    *D=-1/_alpha/(param*param);
  if(DD)
    *DD=2/_alpha/(param*param*param);
  return ret;
#endif
}

//instance
PRJ_BEGIN
template class MetricEnergy<double>;
#ifdef ALL_TYPES
template class MetricEnergy<__float128>;
template class MetricEnergy<mpfr::mpreal>;
#endif
PRJ_END
