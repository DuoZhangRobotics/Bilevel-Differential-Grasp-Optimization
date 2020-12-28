#include "LogBarrierObjEnergy.h"
#include <Environment/ObjMeshGeomCellExact.h>
#include <Articulated/MultiPrecisionSeparatingPlane.h>
#include <stack>

USE_PRJ_NAMESPACE

//LogBarrierEnergy
template <typename T>
LogBarrierObjEnergy<T>::LogBarrierObjEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T d0,T mu):ArticulatedObjective<T>(planner,obj),_d0(d0),_mu(mu) {}
template <typename T>
int LogBarrierObjEnergy<T>::operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h)
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  const std::vector<Node<sizeType,BBox<scalar>>>& bvhObj=_obj.getBVH();
  std::vector<KDOP18<scalar>> bbs=updateBVH(info);
  for(sizeType i=0; i<(sizeType)bbs.size(); i++)
    bbs[i].enlarged(std::to_double(_d0));
  //loop
  std::stack<std::pair<sizeType,sizeType>> ss;
  ss.push(std::make_pair((sizeType)bvhHand.size()-1,(sizeType)bvhObj.size()-1));
  while(!ss.empty()) {
    sizeType idHand=ss.top().first;
    sizeType idObj=ss.top().second;
    ss.pop();
    if(!bbs[idHand].intersect(bvhObj[idObj]._bb))
      continue;
    else if(bvhHand[idHand]._cell && bvhObj[idObj]._cell>=0) {
      //compute derivative
      Mat3X4T tLink=TRANSI(info._TM,idHand);
      Vec3T p=_obj.pss().col(bvhObj[idObj]._cell);
      Vec3T pL=ROT(tLink).transpose()*(p-CTR(tLink));
      Vec3T n,normal;
      Mat3T hessian;
      Vec2i feat;
      const ObjMeshGeomCellExact& distCalc=_planner.dist(idHand);
      if(distCalc.empty())
        continue;
      T dist=distCalc.template closest<T>(pL,n,normal,hessian,feat);
      if(dist<=0)
        return -1;
      else {
        T D,DD,E=MultiPrecisionSeparatingPlane<T>::clog(dist,g?&D:NULL,h?&DD:NULL,_d0,_mu);
        if(!std::isfinite(E))
          return false;
        else if(E==0)
          continue;
        e+=E;
        if(g) {
          p-=CTR(tLink);
          CTRI(*g,idHand)+=-D*ROT(tLink)*normal;
          ROTI(*g,idHand)+=D*p*normal.transpose();
        }
        if(h) {
          Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->coeffRef(0,idHand*12)));
          hessian=(hessian*D+normal*normal.transpose()*DD).eval();
          Mat3T HRT=-hessian*ROT(tLink).transpose();
          for(sizeType r=0; r<3; r++) {
            hBlk.template block<3,3>(r*3,9)+=p*HRT.row(r)-Mat3T::Identity()*normal[r]*D;
            hBlk.template block<3,3>(9,r*3)+=(p*HRT.row(r)-Mat3T::Identity()*normal[r]*D).transpose();
            for(sizeType c=0; c<3; c++) {
              hBlk.template block<3,3>(r*3,c*3)+=p*p.transpose()*hessian(r,c);
            }
          }
          hBlk.template block<3,3>(9,9)-=ROT(tLink)*HRT;
        }
      }
    } else if(bvhHand[idHand]._cell) {
      ss.push(std::make_pair(idHand,bvhObj[idObj]._l));
      ss.push(std::make_pair(idHand,bvhObj[idObj]._r));
    } else if(bvhObj[idObj]._cell>=0) {
      ss.push(std::make_pair(bvhHand[idHand]._l,idObj));
      ss.push(std::make_pair(bvhHand[idHand]._r,idObj));
    } else {
      ss.push(std::make_pair(bvhHand[idHand]._l,bvhObj[idObj]._l));
      ss.push(std::make_pair(bvhHand[idHand]._l,bvhObj[idObj]._r));
      ss.push(std::make_pair(bvhHand[idHand]._r,bvhObj[idObj]._l));
      ss.push(std::make_pair(bvhHand[idHand]._r,bvhObj[idObj]._r));
    }
  }
  return 0;
}
//instance
PRJ_BEGIN
template class LogBarrierObjEnergy<double>;
#ifdef ALL_TYPES
template class LogBarrierObjEnergy<__float128>;
template class LogBarrierObjEnergy<mpfr::mpreal>;
#endif
PRJ_END
