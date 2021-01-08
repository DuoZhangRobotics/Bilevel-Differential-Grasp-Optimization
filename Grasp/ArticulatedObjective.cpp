#include "ArticulatedObjective.h"
#include "GraspPlanner.h"
#include <Utils/DebugGradient.h>
#include <CommonFile/geom/ObjMeshGeomCell.h>

USE_PRJ_NAMESPACE

//ArticulatedObjective
template <typename T>
ArticulatedObjective<T>::ArticulatedObjective(const GraspPlanner<T>& planner,const GraspQualityMetric<T>& obj):_planner(planner),_obj(obj) {}
template <typename T>
int ArticulatedObjective<T>::operator()(const Vec&,const PBDArticulatedGradientInfo<T>&,sizeType,ParallelMatrix<T>&,ParallelMatrix<Mat3XT>*,ParallelMatrix<Mat12XT>*,Vec*,MatT*)
{
  return 0;
}
template <typename T>
int ArticulatedObjective<T>::operator()(const Vec& x,const PBDArticulatedGradientInfo<T>& info,sizeType off,T& e,Vec* g,MatT* h)
{
  sizeType nDOF=info._xM.size();
  ParallelMatrix<Mat3XT> G;
  ParallelMatrix<Mat12XT> H;
  ParallelMatrix<T> E(0);
  Mat3XT tmpG;
  Mat12XT tmpH;
  if(g)
    G.assign(Mat3XT::Zero(3,_planner.body().nrJ()*4));
  if(h)
    H.assign(Mat12XT::Zero(12,_planner.body().nrJ()*12));
  int ret=operator()(x,info,off,E,g?&G:NULL,h?&H:NULL,g,h);
  if(ret<0)
    return ret;
  e+=E.getValue();
  if(g) {
    Eigen::Map<Vec> gMap(g->data(),nDOF);
    info.DTG(_planner.body(),mapM(tmpG=G.getMatrix()),gMap);
  }
  if(h) {
    tmpH=H.getMatrix();
    Eigen::Map<MatT,0,Eigen::OuterStride<>> hMap(h->data(),nDOF,nDOF,h->outerStride());
    Eigen::Map<const MatT,0,Eigen::OuterStride<>> HMap(tmpH.data(),tmpH.rows(),tmpH.cols(),tmpH.outerStride());
    info.toolAB(_planner.body(),HMap,mapM(tmpG=G.getMatrix()),hMap);
  }
  return ret;
}
template <typename T>
void ArticulatedObjective<T>::cons(const Vec&,const PBDArticulatedGradientInfo<T>&,sizeType,sizeType,Vec&,MatT*) {}
template <typename T>
void ArticulatedObjective<T>::setUpdateCache(bool) {}
template <typename T>
sizeType ArticulatedObjective<T>::nrAdditionalDOF() const
{
  return 0;
}
template <typename T>
sizeType ArticulatedObjective<T>::nrCons() const
{
  return 0;
}
template <typename T>
void ArticulatedObjective<T>::debug(Vec x)
{
  sizeType nDOFBody=_planner.body().nrDOF();
  sizeType nDOF=nDOFBody+nrAdditionalDOF();
  DEFINE_NUMERIC_DELTA_T(T)
  Vec dx=Vec::Random(nDOF);
  if(x.size()<nDOF)
    x=concat<Vec,Vec>(x,Vec::Zero(nDOF-x.size()));

  //first evaluate
  T e=0;
  Vec g=Vec::Zero(nDOF),c=Vec::Zero(nrCons());
  MatT h=MatT::Zero(nDOF,nDOF),cjac=MatT::Zero(nrCons(),nDOF);
  PBDArticulatedGradientInfo<T> info(_planner.body(),x);
  setUpdateCache(true);
  if(operator()(x,info,nDOFBody,e,&g,&h)<0) {
    std::ostringstream oss;
    for(sizeType i=0; i<x.size(); i++)
      oss << x[i] << " ";
    WARNINGV("Invalid debug position (x=%s)",oss.str().c_str())
    return;
  }
  if(nrCons()>0)
    cons(x,info,0,nDOFBody,c,&cjac);

  //second evaluate
  T e2=0;
  Vec g2=Vec::Zero(nDOF),c2=Vec::Zero(nrCons());
  PBDArticulatedGradientInfo<T> info2(_planner.body(),x+dx*DELTA);
  setUpdateCache(false);
  operator()(info2._xM,info2,nDOFBody,e2,&g2);
  if(nrCons()>0)
    cons(info2._xM,info2,0,nDOFBody,c2);

  //compare
  DEBUG_GRADIENT(name()+"-G",g.dot(dx),g.dot(dx)-(e2-e)/DELTA)
  DEBUG_GRADIENT(name()+"-H",std::sqrt((h*dx).squaredNorm()),std::sqrt((h*dx-(g2-g)/DELTA).squaredNorm()))
  if(nrCons()>0) {
    DEBUG_GRADIENT("ArticulatedObjective-CJac",std::sqrt((cjac*dx).squaredNorm()),std::sqrt((cjac*dx-(c2-c)/DELTA).squaredNorm()))
  }
}
template <typename T>
std::vector<KDOP18<scalar>> ArticulatedObjective<T>::updateBVH(const PBDArticulatedGradientInfo<T>& info) const
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  std::vector<KDOP18<scalar>> ret(bvhHand.size());
  for(sizeType i=0; i<(sizeType)bvhHand.size(); i++)
    if(bvhHand[i]._cell) {
      Vec3 pos;
      BBox<scalar> bb;
      Mat3 R=ROTI(info._TM,i).unaryExpr([&](const T& in) {
        return (scalar)std::to_double(in);
      });
      Vec3 t=CTRI(info._TM,i).unaryExpr([&](const T& in) {
        return (scalar)std::to_double(in);
      });
      for(sizeType x=0; x<2; x++) {
        pos[0]=x==0?bvhHand[i]._bb._minC[0]:bvhHand[i]._bb._maxC[0];
        for(sizeType y=0; y<2; y++) {
          pos[1]=y==0?bvhHand[i]._bb._minC[1]:bvhHand[i]._bb._maxC[1];
          for(sizeType z=0; z<2; z++) {
            pos[2]=z==0?bvhHand[i]._bb._minC[2]:bvhHand[i]._bb._maxC[2];
            ret[i].setUnion(R*pos+t);
          }
        }
      }
    } else {
      ret[i]=ret[bvhHand[i]._l];
      ret[i].setUnion(ret[bvhHand[i]._r]);
    }
  return ret;
}
//instance
PRJ_BEGIN
template class ArticulatedObjective<double>;
#ifdef ALL_TYPES
template class ArticulatedObjective<__float128>;
template class ArticulatedObjective<mpfr::mpreal>;
#endif
PRJ_END
