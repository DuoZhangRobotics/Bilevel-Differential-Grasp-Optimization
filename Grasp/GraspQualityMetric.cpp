#include "GraspQualityMetric.h"
#include <Utils/DebugGradient.h>
#include <Utils/CrossSpatialUtil.h>
#include <CommonFile/MakeMesh.h>
#include <CommonFile/geom/BVHBuilder.h>
#include <CommonFile/geom/ObjMeshGeomCell.h>
#include <CommonFile/ParallelPoissonDiskSampling.h>
#include <Environment/ObjMeshGeomCellExact.h>

USE_PRJ_NAMESPACE

template <typename T>
GraspQualityMetric<T>::GraspQualityMetric() {}
template <typename T>
void GraspQualityMetric<T>::reset(ObjMesh& obj,T rad,sizeType dRes,const Mat6T& M,T mu)
{
  //move to centroid
  obj.makeUniform();
  if(obj.getVolume()<0)
    obj.insideOut();
  obj.getPos()=-obj.getVolumeCentroid();
  obj.applyTrans();
  _m=obj;
  //sample
  _rad=rad;
  ParallelPoissonDiskSampling sampler(3);
  sampler.setRadius(std::to_double(_rad));
  sampler.sample(obj);
  //pss
  _pss.resize(3,sampler.getPSet().size());
  _nss.resize(3,sampler.getPSet().size());
  for(sizeType i=0; i<sampler.getPSet().size(); i++) {
    _pss.col(i)=sampler.getPSet()[i]._pos.template cast<T>();
    _nss.col(i)=sampler.getPSet()[i]._normal.template cast<T>();
  }
  //gij
  ObjMesh m;
  MakeMesh::makeSphere3D(m,1,dRes);
  std::vector<Vec6T,Eigen::aligned_allocator<Vec6T>> dss;
  for(sizeType i=0; i<(sizeType)m.getV().size(); i++) {
    dss.push_back(concatRow<Vec>(m.getV(i).template cast<T>(),Vec3T::Zero()));
    dss.push_back(concatRow<Vec>(Vec3T::Zero(),m.getV(i).template cast<T>()));
  }
  _gij.resize(_pss.cols(),(sizeType)dss.size());
  for(sizeType r=0; r<_pss.cols(); r++)
    for(sizeType c=0; c<(sizeType)dss.size(); c++)
      _gij(r,c)=computeGij(_pss.col(r),_nss.col(r),dss[c],M,mu);
  //build BVH
  _bvh.resize(_pss.cols());
  for(sizeType i=0; i<_pss.cols(); i++) {
    _bvh[i]._cell=i;
    _bvh[i]._nrCell=1;
    _bvh[i]._bb.setUnion(_pss.col(i).unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
  }
  buildBVH<sizeType>(_bvh,3,-1);
  _distExact.reset(new ObjMeshGeomCellExact(ObjMeshGeomCell(Mat4::Identity(),obj,0,true)));
}
template <typename T>
bool GraspQualityMetric<T>::read(std::istream& is,IOData* dat)
{
  registerType<ObjMeshGeomCellExact>(dat);
  readBinaryData(_bvh,is);
  readBinaryData(_distExact,is,dat);
  Mat3Xd pss;
  readBinaryData(pss,is);
  _pss=pss.template cast<T>();
  Mat3Xd nss;
  readBinaryData(nss,is);
  _nss=nss.template cast<T>();
  _m.readBinary(is);
  Matd gij;
  readBinaryData(gij,is);
  _gij=gij.template cast<T>();
  scalarD rad;
  readBinaryData(rad,is);
  _rad=rad;
  return is.good();
}
template <typename T>
bool GraspQualityMetric<T>::write(std::ostream& os,IOData* dat) const
{
  registerType<ObjMeshGeomCellExact>(dat);
  writeBinaryData(_bvh,os);
  writeBinaryData(_distExact,os,dat);
  Mat3Xd pss=_pss.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(pss,os);
  Mat3Xd nss=_nss.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(nss,os);
  _m.writeBinary(os);
  Matd gij=_gij.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  writeBinaryData(gij,os);
  scalarD rad=std::to_double(_rad);
  writeBinaryData(rad,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> GraspQualityMetric<T>::copy() const
{
  return std::shared_ptr<SerializableBase>(new GraspQualityMetric<T>);
}
template <typename T>
std::string GraspQualityMetric<T>::type() const
{
  return typeid(GraspQualityMetric<T>).name();
}
template <typename T>
void GraspQualityMetric<T>::writeVTK(const std::string& path,T len) const
{
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  for(sizeType i=0; i<_pss.cols(); i++) {
    vss.push_back(_pss.col(i).unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
    vss.push_back((_pss.col(i)+_nss.col(i)*len*_rad).unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
  }
  _m.writeVTK(path+"/mesh.vtk",true);
  VTKWriter<scalar> os("particles",path+"/sample.vtk",true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                 VTKWriter<scalar>::POINT);
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                 VTKWriter<scalar>::LINE);
}
template <typename T>
const std::vector<Node<sizeType,BBox<scalar>>>& GraspQualityMetric<T>::getBVH() const
{
  return _bvh;
}
template <typename T>
T GraspQualityMetric<T>::computeQInfMean(const Vec& w,Vec* g) const
{
  T area=_rad*_rad*M_PI;
  Vec wrench=(_gij.transpose()*w.asDiagonal()).rowwise().sum()*area;
  if(g)
    *g=_gij*Vec::Ones(_gij.cols())*area/wrench.size();
  return wrench.mean();
}
template <typename T>
T GraspQualityMetric<T>::computeQInf(const Vec& w,Vec* g) const
{
  sizeType id;
  T area=_rad*_rad*M_PI;
  Vec wrench=(_gij.transpose()*w.asDiagonal()).rowwise().sum()*area;
  T ret=wrench.minCoeff(&id);
  if(g)
    *g=_gij.col(id)*area;
  return ret;
}
template <typename T>
T GraspQualityMetric<T>::computeQ1(const Vec& w,Vec* g) const
{
  sizeType idMax,idMin;
  Vec wrench=(_gij.transpose()*w.asDiagonal()).rowwise().maxCoeff();
  T ret=wrench.minCoeff(&idMin);
  if(g) {
    Vec tmp=(_gij.col(idMin).array()*w.array()).matrix();
    tmp.maxCoeff(&idMax);
    *g=Vec::Unit(_gij.rows(),idMax)*_gij(idMax,idMin);
  }
  return ret;
}
template <typename T>
const ObjMeshGeomCellExact& GraspQualityMetric<T>::dist() const
{
  return *_distExact;
}
template <typename T>
const typename GraspQualityMetric<T>::Mat3XT& GraspQualityMetric<T>::pss() const
{
  return _pss;
}
template <typename T>
const typename GraspQualityMetric<T>::Mat3XT& GraspQualityMetric<T>::nss() const
{
  return _nss;
}
template <typename T>
const typename GraspQualityMetric<T>::MatT& GraspQualityMetric<T>::gij() const
{
  return _gij;
}
template <typename T>
void GraspQualityMetric<T>::debug(sizeType iter)
{
  MatT h;
  Vec g,g2;
  T Q,Q2;
  DEFINE_NUMERIC_DELTA_T(T)
  for(sizeType i=0; i<iter; i++) {
    Vec w=Vec::Random(_gij.rows());
    Vec dw=Vec::Random(_gij.rows());
    Q=computeQInf(w,&g);
    Q2=computeQInf(w+dw*DELTA);
    DEBUG_GRADIENT("QInf-G",g.dot(dw),g.dot(dw)-(Q2-Q)/DELTA)
    Q=computeQInfMean(w,&g);
    Q2=computeQInfMean(w+dw*DELTA);
    DEBUG_GRADIENT("QInfMean-G",g.dot(dw),g.dot(dw)-(Q2-Q)/DELTA)
  }
}
template <typename T>
T GraspQualityMetric<T>::computeGij(const Vec3T& p,const Vec3T& n,const Vec6T& d,const Mat6T& M,T mu)
{
  Mat3T t=Mat3T::Identity()-n*n.transpose();
  Mat6X3T f=concatRow<MatT>(Mat3T::Identity(),cross<T>(p));
  T wPerp=d.dot(M*f*n),wPara=std::sqrt((d.transpose()*M*f*t).squaredNorm());
  if(wPerp*mu>wPara)
    return wPerp+wPara*wPara/wPerp;
  else return std::max<T>(0,wPerp+mu*wPara);
}
//instance
PRJ_BEGIN
template class GraspQualityMetric<double>;
#ifdef ALL_TYPES
template class GraspQualityMetric<__float128>;
template class GraspQualityMetric<mpfr::mpreal>;
#endif
PRJ_END
