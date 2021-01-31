#include "ArticulatedObjective.h"
#include "GraspPlanner.h"
#include <Utils/DebugGradient.h>
#include <CommonFile/geom/ObjMeshGeomCell.h>

USE_PRJ_NAMESPACE

//ArticulatedObjective
template <typename T>
ArticulatedObjective<T>::ArticulatedObjective(DSSQPObjectiveCompound<T>& obj,const std::string& name,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object)
  :DSSQPObjectiveComponent<T>(obj,name),_planner(planner),_object(object),_info(info)
{
  for(sizeType i=0; i<_planner.body().nrDOF(); i++) {
    ASSERT(i==obj.addVar("GripperDOF"+std::to_string(i),-DSSQPObjectiveCompound<T>::infty(),DSSQPObjectiveCompound<T>::infty(),DSSQPObjectiveCompound<T>::NEW_OR_EXIST)._id)
  }
}
template <typename T>
int ArticulatedObjective<T>::operator()(const Vec&,ParallelMatrix<T>&,ParallelMatrix<Mat3XT>*,ParallelMatrix<Mat12XT>*,Vec*,STrips*)
{
  return 0;
}
template <typename T>
int ArticulatedObjective<T>::operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,SMat* fhess)
{
  int ret=operator()(x,e,g,h,fgrad,fhess?&_tmp:NULL);
  if(fhess) {
    SMat tmp;
    tmp.resize(fhess->rows(),fhess->cols());
    tmp.setFromTriplets(_tmp.begin(),_tmp.end());
    *fhess+=tmp;
  }
  return ret;
}
template <typename T>
int ArticulatedObjective<T>::operator()(const Vec& x,ParallelMatrix<T>& e,ParallelMatrix<Mat3XT>* g,ParallelMatrix<Mat12XT>* h,Vec* fgrad,DMat* fhess)
{
  int ret=operator()(x,e,g,h,fgrad,fhess?&_tmp:NULL);
  if(fhess) {
    SMat tmp;
    tmp.resize(fhess->rows(),fhess->cols());
    tmp.setFromTriplets(_tmp.begin(),_tmp.end());
    *fhess+=tmp;
  }
  return ret;
}
template <typename T>
T ArticulatedObjective<T>::operator()(const Vec& x,Vec* fgrad,STrips* fhess)
{
  sizeType nDOF=_info._xM.size();
  ParallelMatrix<Mat12XT> H;
  ParallelMatrix<Mat3XT> G;
  ParallelMatrix<T> E(0);
  Mat3XT tmpG;
  Mat12XT tmpH;
  if(fgrad)
    G.assign(Mat3XT::Zero(3,_planner.body().nrJ()*4));
  if(fhess)
    H.assign(Mat12XT::Zero(12,_planner.body().nrJ()*12));
  int ret=operator()(x,E,fgrad?&G:NULL,fhess?&H:NULL,fgrad,fhess);
  if(ret<0)
    return E.getValue();
  if(fgrad) {
    Eigen::Map<Vec> gMap(fgrad->data(),nDOF);
    _info.DTG(_planner.body(),mapM(tmpG=G.getMatrix()),gMap);
  }
  if(fhess) {
    tmpH=H.getMatrix();
    MatT fhessD=MatT::Zero(nDOF,nDOF);
    Eigen::Map<const MatT,0,Eigen::OuterStride<>> HMap(tmpH.data(),tmpH.rows(),tmpH.cols(),tmpH.outerStride());
    _info.toolAB(_planner.body(),HMap,mapM(tmpG=G.getMatrix()),mapM(fhessD));
    addBlock(*fhess,0,0,fhessD);
  }
  return E.getValue();
}
//whether modifying objective function expression is allowed
template <typename T>
void ArticulatedObjective<T>::setUpdateCache(const Vec& x,bool)
{
  PBDArticulatedGradientInfo<T>& info=
    const_cast<PBDArticulatedGradientInfo<T>&>(_info);

  sizeType nDOF=_planner.body().nrDOF();
  if(info._xM.size()!=nDOF)
    info.reset(_planner.body(),x.segment(0,nDOF));
  else if(info._xM.segment(0,nDOF)!=x.segment(0,nDOF))
    info.reset(_planner.body(),x.segment(0,nDOF));
}
template <typename T>
const PBDArticulatedGradientInfo<T>& ArticulatedObjective<T>::info() const
{
  return _info;
}
template <typename T>
PBDArticulatedGradientInfo<T>& ArticulatedObjective<T>::info()
{
  return const_cast<PBDArticulatedGradientInfo<T>&>(_info);
}
template <typename T>
const PointCloudObject<T>& ArticulatedObjective<T>::object() const
{
  return _object;
}
template <typename T>
std::vector<KDOP18<scalar>> ArticulatedObjective<T>::updateBVH() const
{
  const std::vector<Node<std::shared_ptr<StaticGeomCell>,BBox<scalar>>>& bvhHand=_planner.body().getGeom().getBVH();
  std::vector<KDOP18<scalar>> ret(bvhHand.size());
  for(sizeType i=0; i<(sizeType)bvhHand.size(); i++)
    if(bvhHand[i]._cell) {
      Vec3 pos;
      BBox<scalar> bb;
      Mat3 R=ROTI(_info._TM,i).unaryExpr([&](const T& in) {
        return (scalar)std::to_double(in);
      });
      Vec3 t=CTRI(_info._TM,i).unaryExpr([&](const T& in) {
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
