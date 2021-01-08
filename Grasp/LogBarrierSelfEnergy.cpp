#include "LogBarrierSelfEnergy.h"
#include "GraspPlanner.h"
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
int LogBarrierSelfEnergy<T>::operator()(const Vec&,const PBDArticulatedGradientInfo<T>& info,sizeType,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec*,MatT*)
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  std::vector<KDOP18<scalar>> bbs=updateBVH(info);
  //update plane
  std::stack<std::pair<sizeType,sizeType>> ss;
  ss.push(std::make_pair((sizeType)bvhHand.size()-1,(sizeType)bvhHand.size()-1));
  while(_updateCache && !ss.empty()) {
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
  std::vector<std::tuple<Vec2i,sizeType,sizeType>> terms;
  for(const std::pair<Vec2i,SeparatingPlane>& sp:_plane)
    for(sizeType pass=0; pass<2; pass++)
      for(sizeType pid=0; pid<sp.second._pss[pass].cols(); pid++)
        terms.push_back(std::make_tuple(sp.first,pass,pid));

  bool valid=true;
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType i=0; i<(sizeType)terms.size(); i++)
      addTerm(valid,terms[i],info,e,g,h);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)terms.size(); i++)
      addTerm(valid,terms[i],info,e,g,h);
  }
  return valid?0:-1;
}
template <typename T>
void LogBarrierSelfEnergy<T>::updatePlanes(const PBDArticulatedGradientInfo<T>& info)
{
  std::vector<std::pair<Vec2i,SeparatingPlane>> pss(_plane.begin(),_plane.end());
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType i=0; i<(sizeType)pss.size(); i++)
      pss[i].second._plane=updatePlane(info,pss[i].first,pss[i].second,true);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)pss.size(); i++)
      pss[i].second._plane=updatePlane(info,pss[i].first,pss[i].second,false);
  }
  _plane.clear();
  _plane.insert(pss.begin(),pss.end());
}
template <typename T>
void LogBarrierSelfEnergy<T>::setUpdateCache(bool update)
{
  _updateCache=update;
}
template <typename T>
std::string LogBarrierSelfEnergy<T>::name() const
{
  return "LogBarrierSelfEnergy(d0="+std::to_string(_d0)+",mu="+std::to_string(_mu)+")";
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
template <typename T>
typename LogBarrierSelfEnergy<T>::Vec4T LogBarrierSelfEnergy<T>::updatePlane(const PBDArticulatedGradientInfo<T>& info,const Vec2i& linkId,const SeparatingPlane& sp,bool callback) const
{
  bool succ;
  Options ops;
  Vec4T p=sp._plane;
  MultiPrecisionSeparatingPlane<T> sol(ops,p,_d0);
  ops.setOptions<MultiPrecisionLQP<T>,T>("muInit",1);
  ops.setOptions<MultiPrecisionLQP<T>,T>("muFinal",1);
  //ops.setOptions<MultiPrecisionLQP<T>,bool>("highPrec",true);
  ops.setOptions<MultiPrecisionLQP<T>,bool>("callback",callback);
  sol.reset(ops);

  Mat3XT pss0=ROTI(info._TM,linkId[0])*sp._pss[0]+CTRI(info._TM,linkId[0])*Vec::Ones(sp._pss[0].cols()).transpose();
  Mat3XT pss1=ROTI(info._TM,linkId[1])*sp._pss[1]+CTRI(info._TM,linkId[1])*Vec::Ones(sp._pss[1].cols()).transpose();
  sol.clearPoints();
  sol.resetPoints(pss0,pss1);
  sol.solve(succ);
  return p;
}
template <typename T>
void LogBarrierSelfEnergy<T>::addTerm(bool& valid,const std::tuple<Vec2i,sizeType,sizeType>& termId,const PBDArticulatedGradientInfo<T>& info,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const
{
  if(!valid)
    return;
  Vec2i links=std::get<0>(termId);
  sizeType pass=std::get<1>(termId);
  sizeType pid=std::get<2>(termId);

  T sgn=pass==0?1:-1;
  const SeparatingPlane& p=_plane.find(links)->second;
  Mat3X4T tLink=TRANSI(info._TM,links[pass]);
  Vec3T PLocal=p._pss[pass].col(pid),P=ROT(tLink)*PLocal+CTR(tLink);
  T tmp=(P.dot(p._plane.template segment<3>(0))+p._plane[3])*sgn;
  T D,DD,E=MultiPrecisionSeparatingPlane<T>::clog(tmp,g?&D:NULL,h?&DD:NULL,_d0,_mu);
  if(!std::isfinite(E)) {
    valid=false;
    return;
  } else if(E==0)
    return;
  e+=E;
  if(g) {
    CTRI(g->getMatrixI(),links[pass])+=D*p._plane.template segment<3>(0)*sgn;
    ROTI(g->getMatrixI(),links[pass])+=D*p._plane.template segment<3>(0)*p._pss[pass].col(pid).transpose()*sgn;
  }
  if(h) {
    Mat3T hessian=DD*p._plane.template segment<3>(0)*p._plane.template segment<3>(0).transpose();
    Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->getMatrixI().coeffRef(0,links[pass]*12)));
    Vec4T cH(PLocal[0],PLocal[1],PLocal[2],1);
    for(sizeType r=0; r<4; r++)
      for(sizeType c=0; c<4; c++)
        hBlk.template block<3,3>(r*3,c*3)+=cH[r]*cH[c]*hessian;
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
