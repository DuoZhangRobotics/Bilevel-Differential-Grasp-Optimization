#include "ConvexLogBarrierSelfEnergy.h"
#include <Articulated/MultiPrecisionSeparatingPlane.h>
#include <Environment/ObjMeshGeomCellExact.h>
#include <Optimizer/QCQPSolverQPOASES.h>
#include "GraspPlanner.h"
#include <Utils/CLog.h>
#include <stack>

USE_PRJ_NAMESPACE

template <typename T>
ConvexLogBarrierSelfEnergy<T>::ConvexLogBarrierSelfEnergy(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,T d0,T mu,bool allPairs)
  :ArticulatedObjective<T>(obj,"ConvexLogBarrierSelfEnergy(d0="+std::to_string(d0)+",mu="+std::to_string(mu)+")",info,planner,object),_allPairs(allPairs),_d0(d0),_mu(mu)
{
  for(sizeType i=0; i<planner.body().nrJ(); i++) {
    sizeType j=planner.body().joint(i)._parent;
    while(j>=0) {
      if(planner.body().joint(j)._M>0)
        break;
      j=planner.body().joint(j)._parent;
    }
    if(j>=0) {
      _exclude.insert(Vec2i(i,j));
      _exclude.insert(Vec2i(j,i));
    }
    _exclude.insert(Vec2i(i,i));
  }
}
template <typename T>
int ConvexLogBarrierSelfEnergy<T>::operator()(const Vec&,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec*,STrips*)
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  std::vector<KDOP18<scalar>> bbs=updateBVH();
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
      else if(!initializePlane(idHand,idHand2))
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
      addTerm(valid,terms[i],e,g,h);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)terms.size(); i++)
      addTerm(valid,terms[i],e,g,h);
  }
  return valid?0:-1;
}
template <typename T>
void ConvexLogBarrierSelfEnergy<T>::updatePlanes()
{
  std::vector<std::pair<Vec2i,SeparatingPlane>> pss(_plane.begin(),_plane.end());
  if(std::is_same<T,mpfr::mpreal>::value) {
    for(sizeType i=0; i<(sizeType)pss.size(); i++)
      pss[i].second._plane=updatePlane(pss[i].first,pss[i].second,true);
  } else {
    OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)pss.size(); i++)
      pss[i].second._plane=updatePlane(pss[i].first,pss[i].second,false);
  }
  _plane.clear();
  _plane.insert(pss.begin(),pss.end());
}
template <typename T>
void ConvexLogBarrierSelfEnergy<T>::setUpdateCache(const Vec& x,bool update)
{
  ArticulatedObjective<T>::setUpdateCache(x,update);
  _updateCache=update;
}
template <typename T>
bool ConvexLogBarrierSelfEnergy<T>::initializePlane(sizeType idL,sizeType idR)
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
    //mL
    //mL.getT()=ROTI(_info._TM,idL).template cast<scalar>();
    //mL.getPos()=CTRI(_info._TM,idL).template cast<scalar>();
    //mL.applyTrans(Vec3::Zero());
    //mL.writeVTK("mL.vtk",true);
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
    //mR
    //mR.getT()=ROTI(_info._TM,idR).template cast<scalar>();
    //mR.getPos()=CTRI(_info._TM,idR).template cast<scalar>();
    //mR.applyTrans(Vec3::Zero());
    //mR.writeVTK("mR.vtk",true);
  }
  //problem
  MatT A=MatT::Zero(sp._pss[0].cols()+sp._pss[1].cols(),4);
  Vec lbA=Vec::Zero(A.rows());
  Vec ubA=Vec::Zero(A.rows());
  sizeType k=0;
  A.col(3).setOnes();
  for(sizeType i=0; i<sp._pss[0].cols(); i++) {
    A.template block<1,3>(k,0)=(ROTI(_info._TM,idL)*sp._pss[0].col(i)+CTRI(_info._TM,idL)).transpose();
    lbA[k]=1;
    ubA[k]= DSSQPObjective<T>::infty();
    k++;
  }
  for(sizeType i=0; i<sp._pss[1].cols(); i++) {
    A.template block<1,3>(k,0)=(ROTI(_info._TM,idR)*sp._pss[1].col(i)+CTRI(_info._TM,idR)).transpose();
    lbA[k]=-DSSQPObjective<T>::infty();
    ubA[k]=-1;
    k++;
  }
  Vec g=Vec4T::Zero(),dwd=Vec4T::Zero();
  MatT h=Mat4T::Identity();
  QCQPSolverQPOASES<T> sol;
  if(sol.solveQP(dwd,h,g,&A,NULL,NULL,&lbA,&ubA,std::vector<Coli,Eigen::aligned_allocator<Coli>>())!=QCQPSolver<T>::SOLVED)
    return false;
  else {
    dwd/=std::sqrt(dwd.template segment<3>(0).squaredNorm());
    sp._plane=dwd;
    _plane[Vec2i(idL,idR)]=sp;
    return true;
  }
}
template <typename T>
typename ConvexLogBarrierSelfEnergy<T>::Vec4T ConvexLogBarrierSelfEnergy<T>::updatePlane(const Vec2i& linkId,const SeparatingPlane& sp,bool callback) const
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

  Mat3XT pss0=ROTI(_info._TM,linkId[0])*sp._pss[0]+CTRI(_info._TM,linkId[0])*Vec::Ones(sp._pss[0].cols()).transpose();
  Mat3XT pss1=ROTI(_info._TM,linkId[1])*sp._pss[1]+CTRI(_info._TM,linkId[1])*Vec::Ones(sp._pss[1].cols()).transpose();
  sol.clearPoints();
  sol.resetPoints(pss0,pss1);
  sol.solve(succ);
  return p;
}
template <typename T>
void ConvexLogBarrierSelfEnergy<T>::addTerm(bool& valid,const std::tuple<Vec2i,sizeType,sizeType>& termId,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h) const
{
  if(!valid)
    return;
  Vec2i links=std::get<0>(termId);
  sizeType pass=std::get<1>(termId);
  sizeType pid=std::get<2>(termId);

  T sgn=pass==0?1:-1;
  const SeparatingPlane& p=_plane.find(links)->second;
  Mat3X4T tLink=TRANSI(_info._TM,links[pass]);
  Vec3T PLocal=p._pss[pass].col(pid),P=ROT(tLink)*PLocal+CTR(tLink);
  T tmp=(P.dot(p._plane.template segment<3>(0))+p._plane[3])*sgn;
  T D,DD,E=clog<T>(tmp,g?&D:NULL,h?&DD:NULL,_d0,_mu);
  e+=E;
  if(!std::isfinite(E)) {
    valid=false;
    return;
  } else if(E==0)
    return;
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
template class ConvexLogBarrierSelfEnergy<double>;
#ifdef ALL_TYPES
template class ConvexLogBarrierSelfEnergy<__float128>;
template class ConvexLogBarrierSelfEnergy<mpfr::mpreal>;
#endif
PRJ_END
