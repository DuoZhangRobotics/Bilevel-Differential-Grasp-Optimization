#include "LogBarrierObjEnergy.h"
#include "GraspPlanner.h"
#include <Utils/CLog.h>
#include <Environment/ConvexHullExact.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <Articulated/MultiPrecisionSeparatingPlane.h>
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
LogBarrierObjEnergy<T>::LogBarrierObjEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T d0,T mu,const bool& useGJK)
  :ArticulatedObjective<T>(obj,"LogBarrierObjEnergy(d0="+std::to_string(d0)+",mu="+std::to_string(mu)+")",info,planner,object),
   _useGJK(useGJK),_updateCache(false),_d0(d0),_d1(d0/2),_mu(mu) {}
template <typename T>
int LogBarrierObjEnergy<T>::operator()(const Vec&,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec*,STrips*)
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  const std::vector<Node<sizeType,BBox<scalarD>>>& bvhObj=_object.getBVH();
  std::vector<KDOP18<scalar>> bbs=updateBVH();
  for(sizeType i=0; i<(sizeType)bbs.size(); i++)
    bbs[i].enlarged(std::to_double(_d0));
  //loop
  std::stack<std::pair<sizeType,sizeType>> ss;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> pairs,feats;
  ss.push(std::make_pair((sizeType)bvhHand.size()-1,(sizeType)bvhObj.size()-1));
  while(!ss.empty()) {
    sizeType idHand=ss.top().first;
    sizeType idObj=ss.top().second;
    ss.pop();
    if(!bbs[idHand].intersect(bvhObj[idObj]._bb))
      continue;
    else if(bvhHand[idHand]._cell && bvhObj[idObj]._cell>=0)
      pairs.push_back(Vec2i(idHand,bvhObj[idObj]._cell));
    else if(bvhHand[idHand]._cell) {
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
  //compute derivative
  bool valid=true;
  feats.resize(pairs.size());
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType i=0; i<(sizeType)pairs.size(); i++)
      addTerm(valid,pairs[i],feats[i],e,g,h);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)pairs.size(); i++)
      addTerm(valid,pairs[i],feats[i],e,g,h);
  }
  if(valid && _updateCache)
    for(sizeType i=0; i<(sizeType)pairs.size(); i++)
      _cache[pairs[i]]=feats[i];
  return valid?0:-1;
}
template <typename T>
void LogBarrierObjEnergy<T>::setUpdateCache(const Vec& x,bool update)
{
  ArticulatedObjective<T>::setUpdateCache(x,update);
  _updateCache=update;
}
template <typename T>
void LogBarrierObjEnergy<T>::addTerm(bool& valid,const Vec2i& termId,Vec2i& feat,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const
{
  const std::vector<Node<sizeType,BBox<scalarD>>>& bvhObj=_object.getBVH();
  if(!valid)
    return;
  sizeType idHand=termId[0];
  sizeType idObj=termId[1];
  Mat3X4T tLink=TRANSI(_info._TM,idHand);
  Vec3T p=_object.pss().col(bvhObj[idObj]._cell);
  Vec3T pL=ROT(tLink).transpose()*(p-CTR(tLink));
  Vec3T n,normal;
  Mat3T hessian;
  const ObjMeshGeomCellExact& distCalc=_planner.dist(idHand);
  const ConvexHullExact* distCalcHull=dynamic_cast<const ConvexHullExact*>(&distCalc);
  if(distCalc.empty())
    return;
  T dist;
  {
    //we will use cache
    Vec2i cacheId(idHand,bvhObj[idObj]._cell);
    if(_useGJK && distCalcHull) {
      dist=distCalcHull->closestGJK<T>(pL,n,normal);
      hessian.setZero();
    } else if(_cache.find(cacheId)==_cache.end())
      dist=distCalc.template closest<T>(pL,n,normal,hessian,feat);
    else dist=distCalc.template closest<T>(pL,n,normal,hessian,feat=_cache.find(cacheId)->second,true);
  }
  if(dist<=0)
    valid=false;
  else {
    T D,DD,E=clog(dist-_d1,g?&D:NULL,h?&DD:NULL,_d0-_d1,_mu);
    e+=E;
    if(!std::isfinite(E)) {
      valid=false;
      return;
    } else if(E==0)
      return;
    if(g) {
      p-=CTR(tLink);
      CTRI(g->getMatrixI(),idHand)+=-D*ROT(tLink)*normal;
      ROTI(g->getMatrixI(),idHand)+=D*p*normal.transpose();
    }
    if(h) {
      Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->getMatrixI().coeffRef(0,idHand*12)));
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
}
//instance
PRJ_BEGIN
template class LogBarrierObjEnergy<double>;
#ifdef ALL_TYPES
template class LogBarrierObjEnergy<__float128>;
template class LogBarrierObjEnergy<mpfr::mpreal>;
#endif
PRJ_END
