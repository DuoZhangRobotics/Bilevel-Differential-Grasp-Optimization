#include <Utils/Scalar.h>
#include "PointCloudObject.h"
#include <Utils/CLog.h>
#include <Utils/Utils.h>
#include <Utils/DebugGradient.h>
#include <Utils/CrossSpatialUtil.h>
#include <CommonFile/Interp.h>
#include <CommonFile/MakeMesh.h>
#include <CommonFile/CameraModel.h>
#include <CommonFile/geom/BVHBuilder.h>
#include <CommonFile/geom/ObjMeshGeomCell.h>
#include <CommonFile/ParallelPoissonDiskSampling.h>
#include <Articulated/MultiPrecisionSeparatingPlane.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Environment/ObjMeshGeomCellExact.h>

USE_PRJ_NAMESPACE

template <typename T>
PointCloudObject<T>::PointCloudObject() {}
template <typename T>
void PointCloudObject<T>::reset(ObjMesh& obj,T rad)
{
  //make uniform
  obj.smooth();
  obj.makeUniform();
  obj.smooth();
  if(obj.getVolume()<0) {
    obj.insideOut();
    obj.smooth();
  }
  //distance
  _m=obj;
  //sample
  samplePoints(rad);
  buildBVH();
}
template <typename T>
void PointCloudObject<T>::resetPointCloud(ArticulatedBody& body,const Vec& x,const Mat4& m,const Mat4& prj,const Vec2i& res)
{
  PBDArticulatedGradientInfo<T> info(body,x);
  _m=body.writeMesh(info._TM.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  }),Joint::MESH);
  //sample
  _rad=1;
  std::vector<Mat4,Eigen::aligned_allocator<Mat4>> tss;
  body.beginUpdateGeom(info._TM.unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  }),tss);
  std::vector<sizeType> idss;
  Vec3 dir,c,X,Y,Z,r,n,normal;
  scalar left=0,right=0,bottom=0,top=0,zNear=0,zFar=0;
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> pss,nss;
  getProjectionMatrixFrustum(prj,left,right,bottom,top,zNear,zFar);
  getViewMatrixFrame<scalar>(m,c,X,Y,Z);
  std::shared_ptr<StaticGeomCell> cell;
  for(sizeType w=0; w<res[0]; w++)
    for(sizeType h=0; h<res[1]; h++) {
      dir=interp1D<Vec3,scalar>(X*left,X*right,(w+0.5f)/res[0])+interp1D<Vec3,scalar>(Y*bottom,Y*top,2*(h+0.5f)/res[1])+Z*zNear;
      dir*=zFar/zNear;
      if(body.getGeom().rayQuery(c,dir,cell,r)) {
        r=transformHomo<scalar>(cell->getT(),r);
        cell->closest(r,n,&normal);
        pss.push_back(r.template cast<T>());
        nss.push_back(normal.template cast<T>());
        for(sizeType j=0; j<body.nrJ(); j++)
          if(body.getGeom().getGPtr(j)==cell)
            idss.push_back((j-2)/2);    //map joint id to object id
      }
    }
  body.endUpdateGeom(tss);
  //assemble and build BVH
  ASSERT(idss.size()==pss.size())
  _pss.resize(3,pss.size());
  _nss.resize(3,nss.size());
  _idss.resize(idss.size());
  for(sizeType c=0; c<(sizeType)pss.size(); c++) {
    _pss.col(c)=pss[c];
    _nss.col(c)=nss[c];
    _idss[c]=idss[c];
  }
  buildBVH();
}
template <typename T>
void PointCloudObject<T>::resetGraspable(ObjMesh& obj,T rad,sizeType dRes,const Mat6T& M,T mu,bool torque)
{
  //make uniform
  obj.smooth();
  obj.makeUniform();
  obj.smooth();
  if(obj.getVolume()<0) {
    obj.insideOut();
    obj.smooth();
  }
  //move to centroid
  obj.getPos()=-obj.getVolumeCentroid();
  obj.applyTrans();
  obj.smooth();
  //distance
  _m=obj;
  _distExact.reset(new ObjMeshGeomCellExact(ObjMeshGeomCell(Mat4::Identity(),obj,0,true)));
  //sample
  samplePoints(rad);
  buildBVH();
  //grasp
  buildGij(dRes,M,mu,torque);
}
template <typename T>
bool PointCloudObject<T>::read(std::istream& is,IOData* dat)
{
  registerType<ObjMeshGeomCellExact>(dat);
  readBinaryData(_bvh,is);
  readBinaryData(_distExact,is,dat);
  readBinaryData(_pss,is);
  readBinaryData(_nss,is);
  readBinaryData(_idss,is);
  _m.readBinary(is);
  readBinaryData(_gij,is);
  readBinaryData(_rad,is);
  return is.good();
}
template <typename T>
bool PointCloudObject<T>::write(std::ostream& os,IOData* dat) const
{
  registerType<ObjMeshGeomCellExact>(dat);
  writeBinaryData(_bvh,os);
  writeBinaryData(_distExact,os,dat);
  writeBinaryData(_pss,os);
  writeBinaryData(_nss,os);
  writeBinaryData(_idss,os);
  _m.writeBinary(os);
  writeBinaryData(_gij,os);
  writeBinaryData(_rad,os);
  return os.good();
}
template <typename T>
std::shared_ptr<SerializableBase> PointCloudObject<T>::copy() const
{
  return std::shared_ptr<SerializableBase>(new PointCloudObject<T>);
}
template <typename T>
std::string PointCloudObject<T>::type() const
{
  return typeid(PointCloudObject<T>).name();
}
template <typename T>
const std::vector<Node<sizeType,BBox<scalarD>>>& PointCloudObject<T>::getBVH() const
{
  return _bvh;
}
template <typename T>
void PointCloudObject<T>::writeVTK(const std::string& path,T len,T normalExtrude) const
{
  Mat3XT PSS=pss(normalExtrude);
  std::vector<scalar> css;
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  for(sizeType i=0; i<PSS.cols(); i++) {
    vss.push_back(PSS.col(i).unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
    vss.push_back((PSS.col(i)+_nss.col(i)*len*_rad).unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
    css.push_back(_idss[i]);
    css.push_back(_idss[i]);
  }
  create(path);
  _m.writeVTK(path+"/mesh.vtk",true);
  VTKWriter<scalar> os("particles",path+"/sample.vtk",true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                 VTKWriter<scalar>::POINT);
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                 VTKWriter<scalar>::LINE);
  os.appendCustomPointData("objectId",css.begin(),css.end());
}
template <typename T>
T PointCloudObject<T>::computeQInfBarrier(const Vec& w,T r,T d0,Vec* g) const
{
  T area=_rad*_rad*M_PI,ret=r,D;
  if(g)
    *g=Vec::Unit(_gij.rows()+1,_gij.rows());
  Vec wrench=(_gij.transpose()*w.asDiagonal()).rowwise().sum()*area;
  for(sizeType i=0; i<wrench.size(); i++) {
    ret-=clog<T>(wrench[i]-r,g?&D:NULL,NULL,d0,1);
    if(!std::isfinite(ret))
      return ret;
    g->segment(0,w.size())-=_gij.col(i)*D*area;
    g->coeffRef(_gij.rows())+=D;
  }
  return ret;
}
template <typename T>
T PointCloudObject<T>::computeQInf(const Vec& w,Vec* g) const
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
T PointCloudObject<T>::computeQ1(const Vec& w,Vec* g) const
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
const ObjMeshGeomCellExact& PointCloudObject<T>::dist() const
{
  return *_distExact;
}
template <typename T>
typename PointCloudObject<T>::Mat3XT PointCloudObject<T>::pss(T normalExtrude) const
{
  return _pss+normalExtrude*_nss;
}
template <typename T>
const typename PointCloudObject<T>::Mat3XT& PointCloudObject<T>::pss() const
{
  return _pss;
}
template <typename T>
const typename PointCloudObject<T>::Mat3XT& PointCloudObject<T>::nss() const
{
  return _nss;
}
template <typename T>
const Coli& PointCloudObject<T>::idss() const
{
  return _idss;
}
template <typename T>
const typename PointCloudObject<T>::MatT& PointCloudObject<T>::gij() const
{
  return _gij;
}
template <typename T>
void PointCloudObject<T>::debug(sizeType iter)
{
  if(_gij.size()==0)
    return;
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
    T r=Q-0.5f;
    Q=computeQInfBarrier(w,r,1,&g);
    Q2=computeQInfBarrier(w+dw*DELTA,r+DELTA,1,&g);
    DEBUG_GRADIENT("QInfBarrier-G",g.dot(concat<Vec,Vec>(dw,Vec::Ones(1))),g.dot(concat<Vec,Vec>(dw,Vec::Ones(1)))-(Q2-Q)/DELTA)
    Q=computeQ1(w,&g);
    Q2=computeQ1(w+dw*DELTA);
    DEBUG_GRADIENT("Q1-G",g.dot(dw),g.dot(dw)-(Q2-Q)/DELTA)
  }
}
//helper
template <typename T>
T PointCloudObject<T>::computeGij(const Vec3T& p,const Vec3T& n,const Vec6T& d,const Mat6T& M,T mu)
{
  Mat3T t=Mat3T::Identity()-n*n.transpose();
  Mat6X3T f=concatRow<MatT>(Mat3T::Identity(),cross<T>(p));
  T wPerp=d.dot(M*f*n),wPara=std::sqrt((d.transpose()*M*f*t).squaredNorm());
  if(wPerp*mu>wPara)
    return wPerp+wPara*wPara/wPerp;
  else return std::max<T>(0,wPerp+mu*wPara);
}
template <typename T>
void PointCloudObject<T>::buildGij(sizeType dRes,const Mat6T& M,T mu,bool torque)
{
  //gij
  ObjMesh m;
  MakeMesh::makeSphere3D(m,1,dRes);
  std::vector<Vec6T,Eigen::aligned_allocator<Vec6T>> dss;
  for(sizeType i=0; i<(sizeType)m.getV().size(); i++) {
    dss.push_back(concatRow<Vec>(m.getV(i).template cast<T>(),Vec3T::Zero()));
    if(torque)
      dss.push_back(concatRow<Vec>(Vec3T::Zero(),m.getV(i).template cast<T>()));
  }
  _gij.resize(_pss.cols(),(sizeType)dss.size());
  for(sizeType r=0; r<_pss.cols(); r++)
    for(sizeType c=0; c<(sizeType)dss.size(); c++)
      _gij(r,c)=computeGij(_pss.col(r),_nss.col(r),dss[c],M,mu);
}
template <typename T>
void PointCloudObject<T>::samplePoints(T rad)
{
  //sample
  _rad=rad;
  ParallelPoissonDiskSampling sampler(3);
  //TODO: I add this line
  sampler.setRadius(std::to_double(_rad));
  sampler.sample(_m, false);
  //pss
  _pss.resize(3,sampler.getPSet().size());
  _nss.resize(3,sampler.getPSet().size());
  for(sizeType i=0; i<sampler.getPSet().size(); i++) {
    _pss.col(i)=sampler.getPSet()[i]._pos.template cast<T>();
    _nss.col(i)=sampler.getPSet()[i]._normal.template cast<T>();
  }
  _idss.setConstant(_pss.cols(),-1);
}
template <typename T>
void PointCloudObject<T>::buildBVH()
{
  //build BVH
  _bvh.assign(_pss.cols(),Node<sizeType,BBox<scalarD>>());
  for(sizeType i=0; i<_pss.cols(); i++) {
    _bvh[i]._cell=i;
    _bvh[i]._nrCell=1;
    _bvh[i]._bb.setUnion(_pss.col(i).unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    }));
  }
  COMMON::buildBVH<sizeType>(_bvh,3,-1);
}
//instance
PRJ_BEGIN
template class PointCloudObject<double>;
#ifdef ALL_TYPES
template class PointCloudObject<__float128>;
template class PointCloudObject<mpfr::mpreal>;
#endif
PRJ_END
