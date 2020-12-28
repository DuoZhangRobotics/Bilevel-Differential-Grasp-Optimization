#include "RigidBodyMass.h"
#include "ArticulatedBody.h"
#include <CommonFile/MakeMesh.h>
#include <CommonFile/geom/StaticGeomCell.h>
#include <Utils/CrossSpatialUtil.h>
#include <Utils/Scalar.h>
#include <Eigen/Eigen>
#include "Joint.h"

USE_PRJ_NAMESPACE

//Joint
Joint::Joint():_parent(-1),_mimic(-1),_color(0.6,0.6,0.6) {}
bool Joint::read(std::istream& is,IOData* dat)
{
  //JointCuda
  readBinaryData(_children,is);
  readBinaryData(_parent,is);
  readBinaryData(_depth,is);
  readBinaryData(_typeJoint,is);
  readBinaryData(_mimic,is);
  readBinaryData(_offDOF,is);
  readBinaryData(_offDDT,is);
  readBinaryData(_limits,is);
  readBinaryData(_control,is);
  readBinaryData(_damping,is);
  readBinaryData(_trans,is);
  //mimic
  readBinaryData(_mult,is);
  readBinaryData(_offset,is);
  //mass
  readBinaryData(_M,is);
  readBinaryData(_MC,is);
  readBinaryData(_MCCT,is);
  readBinaryData(_name,is);
  //sphere approx.
  readBinaryData(_spheres,is);
  readBinaryData(_radGeomColl,is);
  readBinaryData(_radSelfColl,is);
  readBinaryData(_color,is);
  //mesh approx.
  readBinaryData(_mesh,is,dat);
  return is.good();
}
bool Joint::write(std::ostream& os,IOData* dat) const
{
  //JointCuda
  writeBinaryData(_children,os);
  writeBinaryData(_parent,os);
  writeBinaryData(_depth,os);
  writeBinaryData(_typeJoint,os);
  writeBinaryData(_mimic,os);
  writeBinaryData(_offDOF,os);
  writeBinaryData(_offDDT,os);
  writeBinaryData(_limits,os);
  writeBinaryData(_control,os);
  writeBinaryData(_damping,os);
  writeBinaryData(_trans,os);
  //mimic
  writeBinaryData(_mult,os);
  writeBinaryData(_offset,os);
  //mass
  writeBinaryData(_M,os);
  writeBinaryData(_MC,os);
  writeBinaryData(_MCCT,os);
  writeBinaryData(_name,os);
  //sphere approx.
  writeBinaryData(_spheres,os);
  writeBinaryData(_radGeomColl,os);
  writeBinaryData(_radSelfColl,os);
  writeBinaryData(_color,os);
  //mesh approx.
  writeBinaryData(_mesh,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> Joint::copy() const
{
  return std::shared_ptr<SerializableBase>(new Joint);
}
std::string Joint::type() const
{
  return typeid(Joint).name();
}
void Joint::writeVTK(const Mat3X4d& T,VTKWriter<scalar>& os,std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& css,GEOM_TYPE type,const std::vector<unsigned char>* mask) const
{
  { //if(_M > 0) {
    //compute mesh
    ObjMesh mesh=getMesh(type);
    mesh.getT()=ROT(T).cast<scalar>();
    mesh.getPos()=CTR(T).cast<scalar>();
    mesh.applyTrans(Vec3::Zero());
    //write mesh
    os.setRelativeIndex();
    mesh.writeVTK(os);
    if(mask && (type&Joint::MESH)==0) {
      ASSERT_MSG((sizeType)mask->size() == _spheres.cols(),"Mask provided to Joint must have same size as the number of spheres!")
      sizeType nrVPerSphere=mesh.getV().size()/_spheres.cols();
      for(sizeType s=0; s<_spheres.cols(); s++)
        //the 1<<1 bit of mask indicates whether collision will happen, this procedure visualize the mask
        css.insert(css.end(),nrVPerSphere,(mask->at(s)&2)==2?_color:Vec3(Vec3::Ones()-_color));
    } else {
      css.resize(css.size()+mesh.getV().size(),_color);
    }
  }
  if(type&AXIS_INFO) {
    //compute joint variable
    std::vector<Vec3,Eigen::aligned_allocator<Vec3> > vss;
    scalar len=_M > 0 ? _mesh->getBB().getExtent().norm() : 1;
    if(_typeJoint == TRANS_3D) {
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitX()*len).cast<scalar>());
      css.push_back(Vec3::UnitX());
      css.push_back(Vec3::UnitX());
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitY()*len).cast<scalar>());
      css.push_back(Vec3::UnitY());
      css.push_back(Vec3::UnitY());
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitZ()*len).cast<scalar>());
      css.push_back(Vec3::UnitZ());
      css.push_back(Vec3::UnitZ());
    } else if(_typeJoint == TRANS_2D) {
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitX()*len).cast<scalar>());
      css.push_back(Vec3::UnitX());
      css.push_back(Vec3::UnitX());
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitY()*len).cast<scalar>());
      css.push_back(Vec3::UnitY());
      css.push_back(Vec3::UnitY());
    } else if(_typeJoint == TRANS_1D) {
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitX()*len).cast<scalar>());
      css.push_back(Vec3::UnitX());
      css.push_back(Vec3::UnitX());
    } else if(_typeJoint == ROT_3D_EXP || _typeJoint == ROT_3D_XYZ) {
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitX()*len).cast<scalar>());
      css.push_back(Vec3::UnitX());
      css.push_back(Vec3::UnitX());
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitY()*len).cast<scalar>());
      css.push_back(Vec3::UnitY());
      css.push_back(Vec3::UnitY());
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitZ()*len).cast<scalar>());
      css.push_back(Vec3::UnitZ());
      css.push_back(Vec3::UnitZ());
    } else if(_typeJoint == BALL_JOINT) {
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitY()*len).cast<scalar>());
      css.push_back(Vec3::UnitY());
      css.push_back(Vec3::UnitY());
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitZ()*len).cast<scalar>());
      css.push_back(Vec3::UnitZ());
      css.push_back(Vec3::UnitZ());
    } else if(_typeJoint == HINGE_JOINT) {
      vss.push_back(CTR(T).cast<scalar>());
      vss.push_back((CTR(T)+ROT(T)*Vec3d::UnitZ()*len).cast<scalar>());
      css.push_back(Vec3::UnitZ());
      css.push_back(Vec3::UnitZ());
    }
    //write joint variable
    os.setRelativeIndex();
    os.appendPoints(vss.begin(),vss.end());
    os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                   VTKWriter<scalar>::IteratorIndex<Vec3i>((sizeType)vss.size()/2,2,0),
                   VTKWriter<scalar>::LINE,true);
  }
}
ObjMesh Joint::getMesh(GEOM_TYPE type,bool render) const
{
  ObjMesh ball,mesh;
  if(type&MESH) {
    if(_mesh) {
      _mesh->getMesh(mesh,false,render);
      if(_radGeomColl.size() > 0 && _radGeomColl.maxCoeff() < 1E-5f) {
        mesh.smooth();
        mesh=mesh.cutOpen(0.01f);
        mesh.smooth();
      }
    }
    return mesh;
  } else {
    ASSERT((type&SPHERE_GEOM) || (type&SPHERE_SELF))
    MakeMesh::makeSphere3D(ball,1,8);
    mesh.getV().resize(ball.getV().size()*_spheres.cols());
    mesh.getI().resize(ball.getI().size()*_spheres.cols());
    for(sizeType b=0,offV=0,offI=0; b<_spheres.cols(); b++) {
      scalarD rad=(type&SPHERE_GEOM) ? _radGeomColl[b] : _radSelfColl[b];
      Vec3i off=Vec3i::Constant((sizeType)ball.getV().size()*b);
      for(sizeType v=0; v<(sizeType)ball.getV().size(); v++)
        mesh.getV()[offV++]=ball.getV()[v]*rad+_spheres.col(b).cast<scalar>();
      for(sizeType i=0; i<(sizeType)ball.getI().size(); i++)
        mesh.getI()[offI++]=ball.getI()[i]+off;
    }
    return mesh;
  }
}
void Joint::setType(const std::string& type)
{
#define S2T(NAME) if(type == #NAME) {_typeJoint=NAME;return;}
  S2T(TRANS_3D)
  S2T(TRANS_2D)
  S2T(TRANS_1D)
  S2T(ROT_3D_XYZ)
  S2T(ROT_3D_EXP)
  S2T(BALL_JOINT)
  S2T(HINGE_JOINT)
  S2T(FIX_JOINT)
  ASSERT_MSGV(false,"Unknown Type: %s",type.c_str())
#undef S2T
}
std::string Joint::typeToString(JOINT_TYPE type)
{
#define T2S(NAME) if(type == NAME) {return #NAME;}
  T2S(TRANS_3D)
  T2S(TRANS_2D)
  T2S(TRANS_1D)
  T2S(ROT_3D_XYZ)
  T2S(ROT_3D_EXP)
  T2S(BALL_JOINT)
  T2S(HINGE_JOINT)
  T2S(FIX_JOINT)
  ASSERT_MSGV(false,"Unknown Type: %ld",type)
  return "";
#undef T2S
}
void Joint::assemble(scalar rho)
{
  if(!_mesh) {
    _M=0;
    _MC.setZero();
    _MCCT.setZero();
  } else {
    //fill mass
    ObjMesh m;
    _mesh->getMesh(m);
    RigidBodyMass M(m);
    _M=M.getM();
    _MC=M.getMC().cast<scalarD>();
    _MCCT=M.getMCCT().cast<scalarD>();
    if(_M <= 0) {
      m.insideOut();
      RigidBodyMass MInv(m);
      _M=MInv.getM();
      _MC=MInv.getMC().cast<scalarD>();
      _MCCT=MInv.getMCCT().cast<scalarD>();
      Eigen::SelfAdjointEigenSolver<Mat6> eig(MInv.getMassCOM());
      ASSERT(eig.eigenvalues().minCoeff()>=0)
    } else {
      Eigen::SelfAdjointEigenSolver<Mat6> eig(M.getMassCOM());
      ASSERT(eig.eigenvalues().minCoeff()>=0)
    }
    ASSERT(!std::isinf(_M) && !std::isnan(_M))
    ASSERT(isFinite(_MC))
    ASSERT(isFinite(_MCCT))
    //multiply by rho
    _M*=rho;
    _MC*=rho;
    _MCCT*=rho;
  }
}
void Joint::transformMesh(const Mat3X4d& T)
{
  //transform sphere
  if(_spheres.cols() > 0)
    _spheres=getSpheres(T);
  //transform mesh
  if(_mesh)
    _mesh->setT(concatRow(T,Vec4d::Unit(3).transpose()).cast<scalar>()*_mesh->getT());
}
void Joint::transformMesh(const Mat3d& R,const Vec3d& X)
{
  transformMesh((Mat3X4d)concatCol(R,X));
}
void Joint::transformMesh(const Vec3d& X)
{
  transformMesh(Mat3d::Identity(),X);
}
Vec2i Joint::CBegEnd() const
{
  if(_typeJoint == TRANS_3D)
    return Vec2i(_offDOF,_offDOF+3);
  else if(_typeJoint == TRANS_2D)
    return Vec2i(_offDOF,_offDOF+2);
  else if(_typeJoint == TRANS_1D)
    return Vec2i(_offDOF,_offDOF+1);
  else return Vec2i(0,0);
}
Vec2i Joint::RBegEnd() const
{
  if(_typeJoint == ROT_3D_EXP || _typeJoint == ROT_3D_XYZ)
    return Vec2i(_offDOF,_offDOF+3);
  else if(_typeJoint == BALL_JOINT)
    return Vec2i(_offDOF,_offDOF+2);
  else if(_typeJoint == HINGE_JOINT)
    return Vec2i(_offDOF,_offDOF+1);
  else return Vec2i(0,0);
}
sizeType Joint::nrDOF() const
{
  if(_typeJoint == TRANS_3D)
    return 3;
  else if(_typeJoint == TRANS_2D)
    return 2;
  else if(_typeJoint == TRANS_1D)
    return 1;
  else if(_typeJoint == ROT_3D_EXP || _typeJoint == ROT_3D_XYZ)
    return 3;
  else if(_typeJoint == BALL_JOINT)
    return 2;
  else if(_typeJoint == HINGE_JOINT)
    return 1;
  else if(_typeJoint == FIX_JOINT)
    return 0;
  else {
    ASSERT_MSGV(false,"Unknown joint type in %s, name=%s, _typeJoint=%d!",__FUNCTION__,_name.c_str(),_typeJoint)
    return -1;
  }
}
sizeType Joint::nrDDT() const
{
  if(_typeJoint == ROT_3D_EXP || _typeJoint == ROT_3D_XYZ)
    return 9;
  else if(_typeJoint == BALL_JOINT)
    return 1;
  else return 0;
}
Mat6d Joint::getMassC(const Mat3d& R) const
{
  Mat6d ret=getMassC();
  ret.block<3,3>(3,3)=R*(ret.block<3,3>(3,3)*R.transpose()).eval();
  return ret;
}
Mat6d Joint::getMassC() const
{
  Mat6d ret=Mat6d::Zero();
  if(_M <= 0)
    return ret;

  Mat6d M=getMass();
  Vec3d C=getC();
  Mat3d T=C*C.transpose();
  T.diagonal().array()-=C.squaredNorm();
  ret.block<3,3>(0,0)=M.block<3,3>(0,0);
  ret.block<3,3>(3,3)=M.block<3,3>(3,3)+T*_M;
  return ret;
}
Mat6d Joint::getMass(const Mat3d& R) const
{
  Mat6d ret=getMass();
  ret.block<3,3>(3,3)=R*(ret.block<3,3>(3,3)*R.transpose()).eval();
  ret.block<3,3>(0,3)=R*(ret.block<3,3>(0,3)*R.transpose()).eval();
  ret.block<3,3>(3,0)=R*(ret.block<3,3>(3,0)*R.transpose()).eval();
  ret.block<3,3>(0,0)=R*(ret.block<3,3>(0,0)*R.transpose()).eval();
  return ret;
}
Mat6d Joint::getMass() const
{
  Mat6d ret=Mat6d::Zero();
  ret.block<3,3>(0,0).diagonal().setConstant(_M);
  ret.block<3,3>(3,0)=cross<scalarD>(_MC);
  ret.block<3,3>(0,3)=cross<scalarD>(_MC).transpose();
  ret.block<3,3>(3,3)=-_MCCT;
  ret.block<3,3>(3,3).diagonal().array()+=_MCCT.diagonal().sum();
  return ret;
}
Vec3d Joint::getC() const
{
  return _MC/_M;
}
bool Joint::isRotational() const
{
  Vec2i begEnd=RBegEnd();
  return begEnd[1]-begEnd[0];
}
bool Joint::isRoot(const ArticulatedBody& body) const
{
  if(_parent==-1) {
    return true;
  } else {
    //check there are no geometry from this Joint to root Joint
    for(sizeType j=_parent; j>=0; j=body.joint(j)._parent)
      if(body.joint(j)._spheres.cols()>0)
        return false;
    //check Joint type
    Vec2i RP=body.joint(_parent).RBegEnd(),R=RBegEnd();
    Vec2i CP=body.joint(_parent).CBegEnd(),C=CBegEnd();
    if(RP[1]>RP[0] && C[1]>C[0])
      return true;  //this Joint is translational, parent is rotational
    else if(R[1]>R[0] && CP[1]>CP[0])
      return true;  //this Joint is rotational, parent is translational
  }
  return false;
}
Mat3Xd Joint::getSpheres(const Mat3X4d& T) const
{
  Mat3Xd tmp=_spheres;
  getSpheres(T,mapM(tmp));
  return tmp;
}
void Joint::getSpheres(const Mat3X4d& T,Mat3XTM Sss) const
{
  Sss=ROT(T)*_spheres+CTR(T)*Cold::Ones(_spheres.cols()).transpose();
}
std::shared_ptr<StaticGeomCell> Joint::getGeomPtr() const
{
  return _mesh;
}
scalarD Joint::avgRadGeomColl() const
{
  return _radGeomColl.sum()/(scalarD)_radGeomColl.size();
}
scalarD Joint::avgRadSelfColl() const
{
  return _radSelfColl.sum()/(scalarD)_radSelfColl.size();
}
void Joint::loopAllJointTypes(std::function<void(Joint::JOINT_TYPE)> t) {
  sizeType type=NR_JOINT_TYPE-1;
  while(type>0) {
    t((Joint::JOINT_TYPE)type);
    type>>=1;
  }
}
