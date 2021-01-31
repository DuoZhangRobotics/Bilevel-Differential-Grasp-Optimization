#include <CommonFile/Interp.h>
#include <CommonFile/geom/StaticGeomCell.h>
#include "PBDArticulatedGradientInfo.h"
#include "ArticulatedBody.h"
#include "ArticulatedUtils.h"
#include <Utils/RotationUtil.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

//ArticulatedBody
ArticulatedBody::ArticulatedBody() {}
ArticulatedBody::ArticulatedBody(const tinyxml2::XMLElement& pt)
{
  ArticulatedUtils utils(*this);
  utils.assemble(pt);
}
bool ArticulatedBody::read(std::istream& is,IOData* dat)
{
  registerType<StaticGeom>(dat);
  StaticGeom::registerType(dat);
  readBinaryData(_joints,is,dat);
  readBinaryData(_geom,is,dat);
  return is.good();
}
bool ArticulatedBody::write(std::ostream& os,IOData* dat) const
{
  registerType<StaticGeom>(dat);
  StaticGeom::registerType(dat);
  writeBinaryData(_joints,os,dat);
  writeBinaryData(_geom,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> ArticulatedBody::copy() const
{
  return std::shared_ptr<SerializableBase>(new ArticulatedBody);
}
std::string ArticulatedBody::type() const
{
  return typeid(ArticulatedBody).name();
}
void ArticulatedBody::beginUpdateGeom(const Mat3Xd& T,std::vector<Mat4,Eigen::aligned_allocator<Mat4>>& tss)
{
  tss.clear();
  Mat4d TJ=Mat4d::Identity();
  for(sizeType j=0,k=0; j<nrJ(); j++)
    if(_joints[j]._mesh) {
      GETTC(TID,T,j)
      TJ.block<3,4>(0,0)=TID;
      tss.push_back(_joints[j]._mesh->getT());
      _joints[j]._mesh->setT(TJ.cast<scalar>()*_joints[j]._mesh->getT());
      k++;
    }
  _geom->update();
}
void ArticulatedBody::endUpdateGeom(const std::vector<Mat4,Eigen::aligned_allocator<Mat4>>& tss)
{
  for(sizeType j=0,k=0; j<nrJ(); j++)
    if(_joints[j]._mesh) {
      _joints[j]._mesh->setT(tss[k]);
      k++;
    }
  _geom->update();
}
const StaticGeom& ArticulatedBody::getGeom() const
{
  return const_cast<ArticulatedBody&>(*this).getGeom();
}
StaticGeom& ArticulatedBody::getGeom()
{
  if(!_geom)
    _geom.reset(new StaticGeom(3));
  return *_geom;
}
void ArticulatedBody::randomize(sizeType nrLink,bool chain)
{
  sizeType offDOF=0;
  sizeType offDDT=0;
  _joints.resize(nrLink);
  for(sizeType i=0; i<nrLink; i++) {
    if(chain)
      _joints[i]._parent=i-1;
    else _joints[i]._parent=i==0?-1:RandEngine::randI(0,i-1);
    _joints[i]._M=RandEngine::randR(0.5,1);
    _joints[i]._MC=_joints[i]._M*Vec3d::Random();
    _joints[i]._MCCT=_joints[i]._MC*_joints[i]._MC.transpose()/_joints[i]._M;
    _joints[i]._typeJoint=1<<RandEngine::randI(0,7);
    CTR(_joints[i]._trans).setRandom();
    ROT(_joints[i]._trans)=expWGradV<scalarD,Vec3d>(Vec3d::Random()*M_PI);
    _joints[i]._limits.resize(3,_joints[i].nrDOF());
    for(sizeType d=0; d<_joints[i].nrDOF(); d++) {
      _joints[i]._limits(0,d)=-M_PI;
      _joints[i]._limits(1,d)=M_PI;
      _joints[i]._limits(2,d)=1;
    }
    _joints[i]._offDOF=offDOF;
    _joints[i]._offDDT=offDDT;
    offDOF+=_joints[i].nrDOF();
    offDDT+=_joints[i].nrDDT();
  }
  fillChildren();
  //I should assemble _geom in this function,
  //but I decide not to do this to save computation, as randomize is for debug only
}
ObjMesh ArticulatedBody::writeMesh(const Mat3Xd& T,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask) const
{
  ObjMesh m,ret;
  for(sizeType i=0; i<(sizeType)nrJ(); i++)
    if(!jointMask || jointMask->find(i) != jointMask->end()) {
      const Joint& IJ=joint(i);
      if(IJ._M > 0) {
        m=IJ.getMesh(type);
        m.getT()=ROTI(T,i).cast<scalar>();
        m.getPos()=CTRI(T,i).cast<scalar>();
        m.applyTrans(Vec3::Zero());
        ret.addMesh(m,"g"+std::to_string(i));
      }
    }
  return ret;
}
void ArticulatedBody::writeVTK(const Mat3Xd& T,VTKWriter<scalar>& os,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask,const std::vector<std::vector<unsigned char> >* sphereMask) const
{
  std::vector<Vec3,Eigen::aligned_allocator<Vec3> > css;
  ASSERT_MSG(!sphereMask||(sizeType)sphereMask->size()==nrJ(),"Number of sphereMask must be equal to number of joints!")
  for(sizeType i=0; i<(sizeType)nrJ(); i++)
    if(!jointMask || jointMask->find(i) != jointMask->end()) {
      GETTC(JT,T,i)
      _joints[i].writeVTK(JT,os,css,type,sphereMask?&(sphereMask->at(i)):NULL);
    }
  os.appendCustomPointColorData("color",css.begin(),css.end());
}
void ArticulatedBody::writeVTK(const Mat3Xd& T,const std::string& path,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask,const std::vector<std::vector<unsigned char> >* sphereMask) const
{
  VTKWriter<scalar> os("ArticulatedBody",path,true);
  writeVTK(T,os,type,jointMask,sphereMask);
}
void ArticulatedBody::writePov(const Mat3Xd& T,const std::string& path,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask) const
{
  ObjMesh m=writeMesh(T,type,jointMask);
  m.smooth();
  m.writePov(path,true);
}
std::set<sizeType> ArticulatedBody::children(sizeType id,bool direct) const
{
  std::set<sizeType> ret;
  for(sizeType i=0; i<(sizeType)nrJ(); i++)
    if(direct) {
      if(joint(i)._parent==id)
        ret.insert(i);
    } else {
      for(sizeType j=i; j>=0; j=joint(j)._parent)
        if(j == id)
          ret.insert(i);
    }
  return ret;
}
sizeType ArticulatedBody::commonRoot(sizeType id,sizeType id2) const
{
  if(id==id2)
    return id;
  //idChain
  std::vector<sizeType> idChain;
  while(id>=0) {
    idChain.push_back(id);
    id=joint(id)._parent;
  }
  //id2Chain
  std::vector<sizeType> id2Chain;
  while(id2>=0) {
    id2Chain.push_back(id2);
    id2=joint(id2)._parent;
  }
  //compute common root
  sizeType ret=-1;
  while(idChain.back()==id2Chain.back()) {
    ret=idChain.back();
    idChain.pop_back();
    id2Chain.pop_back();
  }
  return ret;
}
sizeType ArticulatedBody::hasJoint(Joint::JOINT_TYPE type) const
{
  for(sizeType i=0; i<nrJ(); i++)
    if(joint(i)._typeJoint==type)
      return i;
  return -1;
}
ArticulatedBody ArticulatedBody::resizeJoints(sizeType nr) const
{
  ArticulatedBody body;
  body._joints=_joints;
  body._joints.resize(nr);
  body.fillChildren();
  return body;
  //I should assemble _geom in this function,
  //but I decide not to do this to save computation, as resizeJoints is used by NESimplifiedDynamics only
}
bool ArticulatedBody::isLeaf(sizeType id) const
{
  return children(id).size() == 1;
}
Cold ArticulatedBody::control() const
{
  Cold ret=Cold::Zero(0);
  for(sizeType j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._control);
  return ret;
}
Cold ArticulatedBody::damping() const
{
  Cold ret=Cold::Zero(0);
  for(sizeType j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._damping);
  return ret;
}
Cold ArticulatedBody::coefLimit() const
{
  Cold ret=Cold::Zero(0);
  for(sizeType j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._limits.row(2).transpose());
  return ret;
}
Cold ArticulatedBody::lowerLimit() const
{
  Cold ret=Cold::Zero(0);
  for(sizeType j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._limits.row(0).transpose());
  return ret;
}
Cold ArticulatedBody::upperLimit() const
{
  Cold ret=Cold::Zero(0);
  for(sizeType j=0; j<nrJ(); j++)
    ret=concat(ret,joint(j)._limits.row(1).transpose());
  return ret;
}
Cold ArticulatedBody::lowerLimit(scalarD infty) const
{
  Cold ret=lowerLimit();
  for(sizeType i=0; i<ret.size(); i++)
    if(!std::isfinite(ret[i]))
      ret[i]=-infty;
  return ret;
}
Cold ArticulatedBody::upperLimit(scalarD infty) const
{
  Cold ret=upperLimit();
  for(sizeType i=0; i<ret.size(); i++)
    if(!std::isfinite(ret[i]))
      ret[i]=infty;
  return ret;
}
const Cold& ArticulatedBody::clampLimit(Cold& x) const
{
  for(sizeType j=0; j<nrJ(); j++) {
    const Joint& J=joint(j);
    for(sizeType i=0; i<J.nrDOF(); i++)
      if(std::isfinite(J._limits(2,i))) {
        scalarD& xEntry=x[J._offDOF+i];
        xEntry=std::min(std::max(xEntry,J._limits(0,i)),J._limits(1,i));
      }
  }
  return x;
}
Cold ArticulatedBody::randomPose(scalarD coef) const
{
  Cold x=Cold::Zero(nrDOF());
  for(sizeType j=0; j<nrJ(); j++) {
    const Joint& J=joint(j);
    for(sizeType i=0; i<J.nrDOF(); i++) {
      scalarD& xEntry=x[J._offDOF+i];
      if(std::isfinite(J._limits(2,i)))
        xEntry=RandEngine::randR(J._limits(0,i),J._limits(1,i));
      else xEntry=RandEngine::randR(-coef,coef);
    }
  }
  return x;
}
Mat3Xd ArticulatedBody::getT(const Cold& x) const
{
  PBDArticulatedGradientInfo<scalarD> info(*this,x);
  return info._TM;
}
void ArticulatedBody::mimic(VecM xMap) const
{
  for(sizeType i=0; i<nrJ(); i++) {
    const Joint& J=joint(i);
    if(J._mimic>=0) {
      const Joint& JM=joint(J._mimic);
      xMap.segment(J._offDOF,J.nrDOF())=xMap.segment(JM._offDOF,JM.nrDOF())*J._mult+Cold::Constant(JM.nrDOF(),J._offset);
    }
  }
}
void ArticulatedBody::mimic(MatT& A,Vec& b,Vec& l,Vec& u) const
{
  sizeType nrDOFReduced=0;
  std::unordered_map<sizeType,sizeType> DOFMap;
  for(sizeType i=0; i<nrJ(); i++)
    if(joint(i)._mimic<0) {
      DOFMap[i]=nrDOFReduced;
      nrDOFReduced+=joint(i).nrDOF();
    }
  A.setZero(nrDOF(),nrDOFReduced);
  b.setZero(nrDOF());
  l.setZero(nrDOFReduced);
  u.setZero(nrDOFReduced);
  for(sizeType i=0; i<nrJ(); i++) {
    const Joint& J=joint(i);
    if(J._mimic>=0) {
      sizeType mimic=J._mimic;
      sizeType offset=J._offset;
      T mult=J._mult;
      while(joint(mimic)._mimic>=0) {
        mimic=joint(mimic)._mimic;
        offset+=mult*joint(mimic)._offset;
        mult*=joint(mimic)._mult;
      }
      A.block(J._offDOF,DOFMap[mimic],J.nrDOF(),J.nrDOF())=MatT::Identity(J.nrDOF(),J.nrDOF())*mult;
      b.segment(J._offDOF,J.nrDOF())=Vec::Constant(J.nrDOF(),offset);
    } else {
      A.block(J._offDOF,DOFMap[i],J.nrDOF(),J.nrDOF())=MatT::Identity(J.nrDOF(),J.nrDOF());
      b.segment(J._offDOF,J.nrDOF())=Vec::Zero(J.nrDOF());
      l.segment(DOFMap[i],J.nrDOF())=lowerLimit().segment(J._offDOF,J.nrDOF());
      u.segment(DOFMap[i],J.nrDOF())=upperLimit().segment(J._offDOF,J.nrDOF());
    }
  }
}
bool ArticulatedBody::movable(sizeType jid,sizeType croot) const
{
  if(jid==croot)
    return false;
  else if(joint(jid)._parent==-1)
    return joint(jid).nrDOF()>0;
  else if(joint(jid).nrDOF()>0)
    return true;
  else return movable(joint(jid)._parent);
}
const Joint& ArticulatedBody::joint(sizeType id) const
{
  return _joints[id];
}
Joint& ArticulatedBody::joint(sizeType id)
{
  return _joints[id];
}
sizeType ArticulatedBody::jointId(const std::string& name) const
{
  for(sizeType k=nrJ()-1; k>=0; k--)
    if(joint(k)._name==name)
      return k;
  ASSERT_MSGV(false,"Cannot file joint with name: %s!",name.c_str())
  return -1;
}
sizeType ArticulatedBody::rootJointId() const
{
  for(sizeType k=nrJ()-1; k>=0; k--)
    if(joint(k).isRoot(*this))
      return k;
  ASSERT_MSG(false,"Cannot find root joint!")
  return -1;
}
sizeType ArticulatedBody::depth() const
{
  sizeType ret=0;
  for(sizeType i=0; i<nrJ(); i++)
    ret=std::max(ret,joint(i)._depth);
  return ret;
}
sizeType ArticulatedBody::nrDOF() const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)nrJ(); i++)
    ret+=joint(i).nrDOF();
  return ret;
}
sizeType ArticulatedBody::nrDDT() const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)nrJ(); i++)
    ret+=joint(i).nrDDT();
  return ret;
}
sizeType ArticulatedBody::nrJ() const
{
  return _joints.size();
}
void ArticulatedBody::fillChildren()
{
  for(sizeType i=0; i<(sizeType)nrJ(); i++)
    joint(i)._children.clear();
  for(sizeType i=0; i<(sizeType)nrJ(); i++)
    if(joint(i)._parent>=0)
      joint(joint(i)._parent)._children.push_back(i);
}
void ArticulatedBody::setRootTrans(const Mat3X4d& t)
{
  sizeType root=rootJointId(),nDOF=0;
  for(sizeType i=0; i<=root; i++)
    nDOF+=joint(i).nrDOF();
  if(nDOF>0) {
    WARNING("Cannot call setT on free-based ArticulatedBody")
  } else {
    joint(root)._trans=t;
  }
}
void ArticulatedBody::debugBase(const std::string& path,T coef) const
{
  recreate(path);
  writeVTK(getT(Cold::Zero(nrDOF())),path+"/ref.vtk",Joint::MESH);
  for(sizeType j=0; j<nrJ(); j++) {
    const Joint& J=joint(j);
    if(J.isRoot(*this))
      for(sizeType d=J._offDOF; d<J._offDOF+J.nrDOF(); d++) {
        Cold DOF=Cold::Unit(nrDOF(),d)*coef;
        writeVTK(getT(DOF),path+"/DOF"+std::to_string(d)+".vtk",Joint::MESH);
      }
  }
}
void ArticulatedBody::addBase(sizeType dim,const Vec3d& planeNormal)
{
  ArticulatedUtils utils(*this);
  utils.addBase(dim,planeNormal);
}
void ArticulatedBody::simplify(sizeType nrDOF)
{
  ArticulatedUtils utils(*this);
  utils.simplify(nrDOF);
}
void ArticulatedBody::eliminateJoint(const std::vector<std::string>& jointName,const Cold& DOF,sizeType nrDebug)
{
  ArticulatedUtils utils(*this);
  utils.eliminateJoint(jointName,DOF,nrDebug);
}
void ArticulatedBody::scaleMass(scalarD coef)
{
  ArticulatedUtils(*this).scaleMass(coef);
}
scalarD ArticulatedBody::totalMass() const
{
  return ArticulatedUtils(const_cast<ArticulatedBody&>(*this)).totalMass();
}
