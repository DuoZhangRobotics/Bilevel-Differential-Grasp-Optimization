#include "LogBarrierSelfEnergy.h"
#include <Articulated/MultiPrecisionSeparatingPlane.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <qpOASES.hpp>
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
LogBarrierSelfEnergy<T>::LogBarrierSelfEnergy(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj,T d0,T mu,bool allPairs)
  :ArticulatedObjective<T>(planner,obj),_allPairs(allPairs),_d0(d0),_mu(mu)
{
  for(sizeType i=0; i<planner.body().nrJ(); i++) {
    const Joint& J=planner.body().joint(i);
    if(J._parent>=0) {
      _exclude.insert(Vec2i(i,J._parent));
      _exclude.insert(Vec2i(J._parent,i));
    }
    _exclude.insert(Vec2i(i,i));
  }
}
template <typename T>
int LogBarrierSelfEnergy<T>::operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g,MatT* h)
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  std::vector<KDOP18<scalar>> bbs=updateBVH(info);
  //update plane
  std::stack<std::pair<sizeType,sizeType>> ss;
  ss.push(std::make_pair((sizeType)bvhHand.size()-1,(sizeType)bvhHand.size()-1));
  while(!ss.empty()) {
    sizeType idHand=ss.top().first;
    sizeType idHand2=ss.top().second;
    ss.pop();
    if(!bbs[idHand].intersect(bbs[idHand2]) && !_allPairs)
      continue;
    else if(bvhHand[idHand]._cell && bvhHand[idHand2]._cell) {
      if(idHand>idHand2)
        continue;
      else if(_exclude.find(Vec2i(idHand,idHand2))!=_exclude.end())
        continue;
      //insert energy
      else if(_plane.find(Vec2i(idHand,idHand2))!=_plane.end())
        continue;
      else if(!initializePlane(info,idHand,idHand2))
        _exclude.insert(Vec2i(idHand,idHand2));
    } else if(bvhHand[idHand]._cell) {
      ss.push(std::make_pair(idHand,bvhHand[idHand2]._l));
      ss.push(std::make_pair(idHand,bvhHand[idHand2]._r));
    } else if(bvhHand[idHand2]._cell) {
      ss.push(std::make_pair(bvhHand[idHand]._l,idHand));
      ss.push(std::make_pair(bvhHand[idHand]._r,idHand));
    } else {
      ss.push(std::make_pair(bvhHand[idHand]._l,bvhHand[idHand2]._l));
      ss.push(std::make_pair(bvhHand[idHand]._l,bvhHand[idHand2]._r));
      ss.push(std::make_pair(bvhHand[idHand]._r,bvhHand[idHand2]._l));
      ss.push(std::make_pair(bvhHand[idHand]._r,bvhHand[idHand2]._r));
    }
  }
  //compute plane gradient
  for(const std::pair<Vec2i,SeparatingPlane>& sp:_plane)
    for(sizeType pass=0; pass<2; pass++) {
      T sgn=pass==0?1:-1;
      const Mat3XT& pss=sp.second._pss[pass];
      Mat3X4T tLink=TRANSI(info._TM,sp.first[pass]);
      Mat3XT Pss=ROT(tLink)*pss+CTR(tLink)*Vec::Ones(pss.cols()).transpose();
      for(sizeType pid=0; pid<pss.cols(); pid++) {
        T tmp=(Pss.col(pid).dot(sp.second._plane.template segment<3>(0))+sp.second._plane[3])*sgn;
        T D,DD,E=MultiPrecisionSeparatingPlane<T>::clog(tmp,g?&D:NULL,h?&DD:NULL,_d0,_mu);
        if(!std::isfinite(E))
          return -1;
        else if(E==0)
          continue;
        e+=E;
        if(g) {
          CTRI(*g,sp.first[pass])+=D*sp.second._plane.template segment<3>(0)*sgn;
          ROTI(*g,sp.first[pass])+=D*sp.second._plane.template segment<3>(0)*pss.col(pid).transpose()*sgn;
        }
        if(h) {
          Mat3T hessian=DD*sp.second._plane.template segment<3>(0)*sp.second._plane.template segment<3>(0).transpose();
          Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->coeffRef(0,sp.first[pass]*12)));
          Vec4T cH(pss(0,pid),pss(1,pid),pss(2,pid),1);
          for(sizeType r=0; r<4; r++)
            for(sizeType c=0; c<4; c++)
              hBlk.template block<3,3>(r*3,c*3)+=cH[r]*cH[c]*hessian;
        }
      }
    }
  return 0;
}
template <typename T>
void LogBarrierSelfEnergy<T>::updatePlanes(const PBDArticulatedGradientInfo<T>& info)
{
  Options ops;
  bool succ;
  Vec4T p;
  MultiPrecisionSeparatingPlane<T> sol(ops,p,_d0);
  ops.setOptions<MultiPrecisionLQP<T>,T>("muInit",1);
  ops.setOptions<MultiPrecisionLQP<T>,T>("muFinal",1);
  ops.setOptions<MultiPrecisionLQP<T>,bool>("callback",false);
  sol.reset(ops);
  for(const std::pair<Vec2i,SeparatingPlane>& pss:_plane) {
    p=pss.second._plane;
    Mat3XT pss0=ROTI(info._TM,pss.first[0])*pss.second._pss[0]+CTRI(info._TM,pss.first[0])*Vec::Ones(pss.second._pss[0].cols()).transpose();
    Mat3XT pss1=ROTI(info._TM,pss.first[1])*pss.second._pss[1]+CTRI(info._TM,pss.first[1])*Vec::Ones(pss.second._pss[1].cols()).transpose();
    sol.clearPoints();
    sol.resetPoints(pss0,pss1);
    sol.solve(succ);
    _plane.find(pss.first)->second._plane=p;
  }
}
template <typename T>
bool LogBarrierSelfEnergy<T>::initializePlane(const PBDArticulatedGradientInfo<T>& info,sizeType idL,sizeType idR)
{
  SeparatingPlane sp;
  {
    const ObjMeshGeomCell& cL=*std::static_pointer_cast<ObjMeshGeomCell>(_planner.body().getGeom().getBVH()[idL]._cell);
    ObjMesh mL;
    cL.getMesh(mL);
    sp._pss[0].resize(3,(sizeType)mL.getV().size());
    for(sizeType i=0; i<(sizeType)mL.getV().size(); i++)
      sp._pss[0].col(i)=mL.getV(i).template cast<T>();
    if(sp._pss[0].size()==0)
      return false;
  }
  {
    const ObjMeshGeomCell& cR=*std::static_pointer_cast<ObjMeshGeomCell>(_planner.body().getGeom().getBVH()[idR]._cell);
    ObjMesh mR;
    cR.getMesh(mR);
    sp._pss[1].resize(3,mR.getV().size());
    for(sizeType i=0; i<(sizeType)mR.getV().size(); i++)
      sp._pss[1].col(i)=mR.getV(i).template cast<T>();
    if(sp._pss[1].size()==0)
      return false;
  }
  //problem
  Eigen::Matrix<double,-1,-1,Eigen::RowMajor> A;
  A.setZero(sp._pss[0].cols()+sp._pss[1].cols(),4);
  Cold lbA=Cold::Zero(A.rows());
  Cold ubA=Cold::Zero(A.rows());
  sizeType k=0;
  A.col(3).setOnes();
  for(sizeType i=0; i<sp._pss[0].cols(); i++) {
    A.block<1,3>(k,0)=(ROTI(info._TM,idL)*sp._pss[0].col(i)+CTRI(info._TM,idL)).transpose().unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    });
    lbA[k]=1;
    ubA[k]= qpOASES::INFTY;
    k++;
  }
  for(sizeType i=0; i<sp._pss[1].cols(); i++) {
    A.block<1,3>(k,0)=(ROTI(info._TM,idR)*sp._pss[1].col(i)+CTRI(info._TM,idR)).transpose().unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    });
    lbA[k]=-qpOASES::INFTY;
    ubA[k]=-1;
    k++;
  }
  Vec4d g=Vec4d::Zero();
  Mat4d h=Mat4d::Identity();
  qpOASES::int_t nWSR=10000;
  qpOASES::SQProblem prob(4,A.rows());
  prob.init(h.data(),g.data(),A.data(),NULL,NULL,lbA.data(),ubA.data(),nWSR);
  if(!prob.isSolved()) {
    return false;
  } else {
    Vec4d dwd;
    prob.getPrimalSolution(dwd.data());
    dwd/=dwd.segment<3>(0).norm();
    sp._plane=dwd.template cast<T>();
    _plane[Vec2i(idL,idR)]=sp;
    return true;
  }
}
//instance
PRJ_BEGIN
template class LogBarrierSelfEnergy<double>;
#ifdef ALL_TYPES
template class LogBarrierSelfEnergy<__float128>;
template class LogBarrierSelfEnergy<mpfr::mpreal>;
#endif
PRJ_END
