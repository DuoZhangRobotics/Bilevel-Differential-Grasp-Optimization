#ifndef JOINT_FUNC_H
#define JOINT_FUNC_H

#include "Joint.h"

PRJ_BEGIN

template <typename T>
struct JointFunc
{
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  //PBTO
  static Mat3X4T TDTDDT(const Joint& J,VecCM x,Mat3XTM DT,Mat3XTM DDT);
  static void DDT(const Joint& J,std::function<void(sizeType,sizeType,T)> hess,Mat3XTCM DDT,const Vec3T& coefRss);
  static void DDT(const Joint& J,MatTM hess,Mat3XTCM DDT,const Vec3T& coefRss);
  static Vec3T DDT(const Joint& J,Mat3XTCM DDT,sizeType R,sizeType C);
  static void DDDTLambda(const Joint& J,VecCM x,VecCM lambda,Mat3XTM DDDTLambda);
  static void DDTLambda(const Joint& J,Mat3XTM hess,Mat3XTCM DDT,const VecCM& lambda);
  //NE
  static void JCALC(const Joint& J,const VecCM q,const VecCM dq,const VecCM ddq,
                    Mat6TM XJ,Mat6XTM S,Vec6TM vJ,Vec6TM dvJ,Mat6XTM DvJDq,Mat6XTM DdvJDq,Mat6XTM DdvJDdq);
  static void DSTDqf(const Joint& J,const VecCM q,MatTM DtauDq,const Vec6T& f);
  static void DSDqf(const Joint& J,const VecCM q,MatTM DtauDq,VecCM f);
};

PRJ_END

#endif
