#include "SimplifiedDynamics.h"
#include <Optimizer/KNInterface.h>
#include <Optimizer/IPOPTInterface.h>
#include <CommonFile/MakeMesh.h>
#include <CommonFile/Interp.h>
#include <CommonFile/Timing.h>
#include <Utils/RotationUtil.h>
#include <Utils/Utils.h>
#include <stack>

USE_PRJ_NAMESPACE

//EndEffectorBounds
EndEffectorBounds::EndEffectorBounds() {}
EndEffectorBounds::EndEffectorBounds(sizeType jid)
{
  _JID.push_back(jid);
}
EndEffectorBounds::EndEffectorBounds(sizeType jid,const Vec3d& localPos,scalarD phi0):_localPos(localPos),_phi0(phi0)
{
  _JID.push_back(jid);
}
bool EndEffectorBounds::exists(const BBox<sizeType>& bb) const
{
  for(sizeType x=bb._minC[0]; x<=bb._maxC[0]; x++)
    for(sizeType y=bb._minC[1]; y<=bb._maxC[1]; y++)
      for(sizeType z=bb._minC[2]; z<=bb._maxC[2]; z++)
        if(_indices.find(Vec3i(x,y,z))==_indices.end())
          return false;
  return true;
}
bool EndEffectorBounds::isBetterBB(const BBox<sizeType>& bb,const std::vector<sizeType>& sym) const
{
  if(bb.getExtent()[0]*bb.getExtent()[1]*bb.getExtent()[2] > _bb.getExtent()[0]*_bb.getExtent()[1]*_bb.getExtent()[2])
    return true;
  BBox<sizeType> _bbCloseFar=closeFar(_bb,sym);
  BBox<sizeType>  bbCloseFar=closeFar( bb,sym);
  for(sizeType d=0; d<3; d++)
    if(_bbCloseFar._minC[d]<bbCloseFar._minC[d])
      return true;
    else return false;
  for(sizeType d=0; d<3; d++)
    if(_bbCloseFar._maxC[d]<bbCloseFar._maxC[d])
      return true;
    else return false;
  return false;
}
bool EndEffectorBounds::sameOrthant(const Vec3i& id,const Vec3d& pos,const std::vector<sizeType>& sym) const
{
  for(sizeType ds=0; ds<(sizeType)sym.size(); ds++) {
    sizeType d=sym[ds];
    if(id[d]*_res*pos[d]<=0)
      return false;
  }
  return true;
}
BBox<sizeType> EndEffectorBounds::closeFar(const BBox<sizeType>& in,const std::vector<sizeType>& sym)
{
  BBox<sizeType> bb=in;
  for(sizeType ds=0; ds<(sizeType)sym.size(); ds++) {
    sizeType d=sym[ds];
    bb._minC[d]=std::min(std::abs(in._minC[d]),std::abs(in._maxC[d]));
    bb._maxC[d]=std::max(std::abs(in._minC[d]),std::abs(in._maxC[d]));
  }
  return bb;
}
sizeType EndEffectorBounds::nrDOF(const ArticulatedBody& body) const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)_JID.size(); i++)
    ret+=body.joint(_JID[i]).nrDOF();
  return ret;
}
Cold EndEffectorBounds::getDOF(const Cold& pt,const Vec3d& ctr) const
{
  Vec3d frac=(pt-ctr)/_res;
  Vec3i id=frac.cast<sizeType>().cwiseMin(_bb._maxC-Vec3i::Ones()).cwiseMax(_bb._minC);
  Cold id000=_indices.find(id+Vec3i(0,0,0))->second;
  Cold id100=_indices.find(id+Vec3i(1,0,0))->second;
  Cold id010=_indices.find(id+Vec3i(0,1,0))->second;
  Cold id110=_indices.find(id+Vec3i(1,1,0))->second;
  Cold id001=_indices.find(id+Vec3i(0,0,1))->second;
  Cold id101=_indices.find(id+Vec3i(1,0,1))->second;
  Cold id011=_indices.find(id+Vec3i(0,1,1))->second;
  Cold id111=_indices.find(id+Vec3i(1,1,1))->second;
  return interp3D<Cold,scalarD>(id000,id100,id010,id110,
                                id001,id101,id011,id111,
                                frac[0]-id[0],frac[1]-id[1],frac[2]-id[2]);
}
Vec3d EndEffectorBounds::globalPosAll(const ArticulatedBody& body,const Cold& xAll) const
{
  Cold x=Cold::Zero(nrDOF(body));
  for(sizeType i=0,id=0; i<(sizeType)_JID.size(); i++) {
    const Joint& J=body.joint(_JID[i]);
    for(sizeType j=0; j<J.nrDOF(); j++,id++)
      x[id]=xAll[J._offDOF+j];
  }
  return globalPos(body,x);
}
Vec3d EndEffectorBounds::globalPos(const ArticulatedBody& body,const Cold& x) const
{
#ifdef OPTIMIZER_SUPPORT
  Cold dof=Cold::Zero(body.nrDOF());
  if(x.size()==dof.size())
    dof=x;
  else dof=InverseKinematicsObj<scalarD>::assignDOF(body,*this,x,NULL);
  PBDArticulatedGradientInfo<scalarD> info(body,dof);
  return ROTI(info._TM,_JID[0])*_localPos+CTRI(info._TM,_JID[0]);
#else
  FUNCTION_NOT_IMPLEMENTED
  return Vec3d::Zero();
#endif
}
Vec3d EndEffectorBounds::randomPos(const Vec3d& ctr) const
{
  BBox<scalarD> bb;
  bb._minC=_bb._minC.cast<scalarD>()*_res+ctr;
  bb._maxC=_bb._maxC.cast<scalarD>()*_res+ctr;
  return Vec3d(interp1D(bb._minC[0],bb._maxC[0],RandEngine::randR01()),
               interp1D(bb._minC[1],bb._maxC[1],RandEngine::randR01()),
               interp1D(bb._minC[2],bb._maxC[2],RandEngine::randR01()));
}
sizeType EndEffectorBounds::jointId() const
{
  return _JID[0];
}
bool EndEffectorBounds::operator==(const EndEffectorBounds& other) const
{
  return _JID[0]==other._JID[0] && _localPos==other._localPos;
}
bool EndEffectorBounds::read(std::istream& is,IOData* dat)
{
  readBinaryData(_indices,is,dat);
  readBinaryData(_JID,is,dat);
  readBinaryData(_bb,is,dat);
  readBinaryData(_localPos,is,dat);
  readBinaryData(_phi0,is,dat);
  readBinaryData(_res,is,dat);
  return is.good();
}
bool EndEffectorBounds::write(std::ostream& os,IOData* dat) const
{
  writeBinaryData(_indices,os,dat);
  writeBinaryData(_JID,os,dat);
  writeBinaryData(_bb,os,dat);
  writeBinaryData(_localPos,os,dat);
  writeBinaryData(_phi0,os,dat);
  writeBinaryData(_res,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> EndEffectorBounds::copy() const
{
  return std::shared_ptr<SerializableBase>(new EndEffectorBounds());
}
std::string EndEffectorBounds::type() const
{
  return typeid(EndEffectorBounds).name();
}
#ifdef OPTIMIZER_SUPPORT
//InverseKinematicsObj
template <typename T>
InverseKinematicsObj<T>::InverseKinematicsObj(DSSQPObjectiveCompound<T>& obj,const ArticulatedBody& body,const EndEffectorBounds& ee,sizeType off,bool allDOF,const Vec* init)
  :DSSQPObjectiveComponent<T>(obj,"InverseKinematicsObj"+std::to_string(ee.jointId()),true),_body(body),_ee(ee)
{
  for(sizeType i=0; i<3; i++) {
    _gl.push_back(0);
    _gu.push_back(0);
  }
  if(allDOF) {
    Cold L=_body.lowerLimit(DSSQPObjective<scalarD>::infty());
    Cold U=_body.upperLimit(DSSQPObjective<scalarD>::infty());
    for(sizeType i=0; i<L.size(); i++) {
      _vid.push_back(obj.addVar(DOF(off,i),L[i],U[i])._id);
      if(init)
        obj.setVarInit(_vid.back(),(*init)[i]);
    }
  } else {
    for(sizeType i=0,id=0; i<(sizeType)ee._JID.size(); i++) {
      const Joint& J=_body.joint(ee._JID[i]);
      for(sizeType j=0; j<J.nrDOF(); j++,id++) {
        _vid.push_back(obj.addVar(DOF(off,allDOF?J._offDOF+j:id),J._limits(0,j),J._limits(1,j))._id);
        if(init)
          obj.setVarInit(_vid.back(),(*init)[J._offDOF+j]);
      }
    }
  }
  _H.setIdentity();
}
template <typename T>
int InverseKinematicsObj<T>::operator()(const Vec& x,Vec& fvec,STrips* fjac)
{
  _info.reset(_body,assignDOF(_body,_ee,x,&_vid));
  Vec3T pos=ROTI(_info._TM,_ee._JID[0])*_ee._localPos.template cast<T>()+CTRI(_info._TM,_ee._JID[0]);
  Vec3T dir=_H*(pos-_target);
  for(sizeType d=0; d<3; d++)
    fvec[_offset+d]=dir[d];
  if(fjac) {
    Mat3XT DTG=Mat3XT::Zero(3,_body.nrDOF());
    for(sizeType d=0; d<3; d++) {
      Vec DTGRow=Vec::Zero(_body.nrDOF());
      Mat3XT G=Mat3XT::Zero(3,_body.nrJ()*4);
      ROTI(G,_ee._JID[0])=Vec3T::Unit(d)*_ee._localPos.transpose().template cast<T>();
      CTRI(G,_ee._JID[0])=Vec3T::Unit(d);
      _info.DTG(_body,mapM(G),mapV(DTGRow));
      DTG.row(d)=DTGRow;
    }
    DTG=(_H*DTG).eval();
    for(sizeType d=0; d<3; d++) {
      Vec DTGRow=DTG.row(d);
      assignDOF(_body,_ee,*fjac,DTGRow,_offset+d,_vid);
    }
  }
  return 0;
}
template <typename T>
typename InverseKinematicsObj<T>::Vec InverseKinematicsObj<T>::assignDOF(const ArticulatedBody& body,const EndEffectorBounds& ee,const Vec& x,const std::vector<sizeType>* vid) {
  Vec dof=Vec::Zero(body.nrDOF());
  if(vid && (sizeType)vid->size()==body.nrDOF()) {
    for(sizeType i=0; i<body.nrDOF(); i++)
      dof[i]=x[vid?(*vid)[i]:i];
  } else {
    for(sizeType i=0,id=0; i<(sizeType)ee._JID.size(); i++) {
      const Joint& J=body.joint(ee._JID[i]);
      for(sizeType j=0; j<J.nrDOF(); j++,id++)
        dof[J._offDOF+j]=x[vid?(*vid)[id]:id];
    }
  }
  return dof;
}
template <typename T>
void InverseKinematicsObj<T>::assignDOF(const ArticulatedBody& body,const EndEffectorBounds& ee,STrips& fjac,const Vec& G,sizeType offset,const std::vector<sizeType>& vid) {
  if((sizeType)vid.size()==body.nrDOF()) {
    for(sizeType i=0; i<body.nrDOF(); i++)
      fjac.push_back(STrip(offset,vid[i],G[i]));
  } else {
    for(sizeType i=0,id=0; i<(sizeType)ee._JID.size(); i++) {
      const Joint& J=body.joint(ee._JID[i]);
      for(sizeType j=0; j<J.nrDOF(); j++,id++)
        fjac.push_back(STrip(offset,vid[id],G[J._offDOF+j]));
    }
  }
}
template <typename T>
std::string InverseKinematicsObj<T>::DOF(sizeType r,sizeType c) {
  return std::string("DOF(")+(r>=0?std::to_string(r):"-")+","+(c>=0?std::to_string(c):"-")+")";
}
//VariationalSimilarity
template <typename T>
VariationalSimilarityObj<T>::VariationalSimilarityObj(DSSQPObjectiveCompound<T>& obj,const ArticulatedBody& body,const Vec& init,sizeType off,bool rotBase,bool transBase)
  :DSSQPObjectiveComponent<T>(obj,"VariationalSimilarityObj"),_body(body)
{
  Cold L=_body.lowerLimit(DSSQPObjective<scalarD>::infty());
  Cold U=_body.upperLimit(DSSQPObjective<scalarD>::infty());
  for(sizeType i=0; i<L.size(); i++) {
    _vid.push_back(obj.addVar(InverseKinematicsObj<T>::DOF(off,i),L[i],U[i])._id);
    obj.setVarInit(_vid.back(),init[i]);
  }
  for(sizeType i=0; i<body.nrJ(); i++) {
    const Joint& J=_body.joint(i);
    for(sizeType j=0; j<J.nrDOF(); j++)  {
      _vidBase.push_back(obj.addVar("Base"+InverseKinematicsObj<T>::DOF(off,J._offDOF+j),L[J._offDOF+j],U[J._offDOF+j])._id);
      obj.setVarInit(_vidBase.back(),init[J._offDOF+j]);
    }
    if(!J.isRotational() && !transBase) {
      for(sizeType j=0; j<J.nrDOF(); j++) {
        obj.addVar(_vid[J._offDOF+j])._l=init[J._offDOF+j];
        obj.addVar(_vid[J._offDOF+j])._u=init[J._offDOF+j];
        obj.addVar(_vidBase[J._offDOF+j])._l=init[J._offDOF+j];
        obj.addVar(_vidBase[J._offDOF+j])._u=init[J._offDOF+j];
      }
    }
    if(J.isRotational() && !rotBase) {
      for(sizeType j=0; j<J.nrDOF(); j++) {
        obj.addVar(_vid[J._offDOF+j])._l=init[J._offDOF+j];
        obj.addVar(_vid[J._offDOF+j])._u=init[J._offDOF+j];
        obj.addVar(_vidBase[J._offDOF+j])._l=init[J._offDOF+j];
        obj.addVar(_vidBase[J._offDOF+j])._u=init[J._offDOF+j];
      }
    }
    if(J._M>0)
      break;
  }
  _info[1].reset(_body,init);
}
template <typename T>
T VariationalSimilarityObj<T>::operator()(const Vec& x,Vec* fgrad)
{
  T ret=0.0;
  Mat3X4T A;
  Mat3XT G,GB,MRR,MRt,MtR,Mtt;
  sizeType nrJ=_body.nrJ();
  //INFO0
  Vec DOF=Vec::Zero(_body.nrDOF());
  for(sizeType i=0; i<(sizeType)_vid.size(); i++)
    DOF[i]=x[_vid[i]];
  _info[0].reset(_body,DOF);
  //INFO1
  DOF=_info[1]._xM;
  for(sizeType i=0; i<(sizeType)_vidBase.size(); i++)
    DOF[i]=x[_vidBase[i]];
  _info[1].reset(_body,DOF);
  if(false) {
    ret+=(_info[0]._xM-_info[1]._xM).squaredNorm();
    DOF=(_info[0]._xM-_info[1]._xM)*2;
    if(fgrad) {
      for(sizeType i=0; i<(sizeType)_vid.size(); i++)
        (*fgrad)[_vid[i]]=DOF[i];
      for(sizeType i=0; i<(sizeType)_vidBase.size(); i++)
        (*fgrad)[_vidBase[i]]=-DOF[i];
    }
  } else {
    if(fgrad) {
      fgrad->setZero(x.size());
      G.setZero(3,4*nrJ);
      GB.setZero(3,4*nrJ);
    }
    for(sizeType k=0; k<nrJ; k++) {
      const Joint& J=_body.joint(k);
      const Mat3T PPT=J._MCCT.template cast<T>();
      const Vec3T P=J._MC.template cast<T>();
      A=TRANSI(_info[0]._TM,k)-TRANSI(_info[1]._TM,k);
      ret+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace();
      if(fgrad) {
        ROTI(G,k)+=ROT(A)*PPT+CTR(A)*P.transpose();
        CTRI(G,k)+=CTR(A)*J._M+ROT(A)*P;
      }
    }
    if(fgrad) {
      G*=2;
      DOF.setZero();
      _info[0].DTG(_body,mapM(GB=G),mapV(DOF));
      for(sizeType i=0; i<(sizeType)_vid.size(); i++)
        (*fgrad)[_vid[i]]=DOF[i];
      DOF.setZero();
      _info[1].DTG(_body,mapM(GB=G),mapV(DOF));
      for(sizeType i=0; i<(sizeType)_vidBase.size(); i++)
        (*fgrad)[_vidBase[i]]=-DOF[i];
    }
  }
  return ret;
}
#endif
//SimplifiedDynamics
SimplifiedDynamics::SimplifiedDynamics() {}
SimplifiedDynamics::SimplifiedDynamics(std::shared_ptr<ArticulatedBody> body,scalarD res,sizeType symDir,const Vec3d& zRange,scalarD zMargin):_body(body)
{
  ASSERT_MSG(body->nrJ()>1 && body->joint(1).isRoot(*body),"Body joint #1 is not root!")
  ASSERT_MSG(body->joint(0)._typeJoint==Joint::TRANS_3D,"Body joint #0 is not translational!")
  ASSERT_MSG(body->joint(1)._typeJoint==Joint::ROT_3D_XYZ||body->joint(1)._typeJoint==Joint::ROT_3D_EXP,"Body joint #1 is not rotational!")
  //detect end-effector
  _rootId=1;
  ObjMesh eeMesh;
  for(sizeType jid=0; jid<body->nrJ(); jid++)
    if(body->children(jid,true).empty()) {
      EndEffectorBounds ee;
      detectEndEffector(*body,jid,ee._localPos,ee._phi0,zRange);
      for(sizeType j=jid; j>1; j=body->joint(j)._parent)
        ee._JID.push_back(j);
      ee._res=ee._phi0*res;
      INFOV("EE%ld res=%f #DOF=%ld!",jid,ee._res,ee.nrDOF(*_body))
      _ee.push_back(ee);
      //draw EE
      ObjMesh s;
      MakeMesh::makeSphere3D(s,ee._phi0,16);
      PBDArticulatedGradientInfo<scalarD> info(*body,Cold::Zero(body->nrDOF()));
      s.getPos()=(ROTI(info._TM,jid)*ee._localPos+CTRI(info._TM,jid)).cast<scalar>();
      s.applyTrans();
      eeMesh.addMesh(s);
    }
  INFOV("Parameter: res=%lf, symDir=%ld, zRange=(%lf,%lf,%lf) zMargin=%lf!",res,symDir,zRange[0],zRange[1],zRange[2],zMargin)
  eeMesh.writeVTK("EndEffector.vtk",true);
  //prune symmetry
  if(_ee.size()==3) {
    Vec3d pos0=_ee[0].globalPos(*_body,Cold::Zero(_ee[0].nrDOF(*_body)));
    Vec3d pos1=_ee[1].globalPos(*_body,Cold::Zero(_ee[1].nrDOF(*_body)));
    Vec3d pos2=_ee[2].globalPos(*_body,Cold::Zero(_ee[2].nrDOF(*_body)));
    scalarD score01=(pos0-pos1).cwiseAbs().minCoeff();
    scalarD score02=(pos0-pos2).cwiseAbs().minCoeff();
    scalarD score12=(pos1-pos2).cwiseAbs().minCoeff();
    if(score01<score02 && score01<score12)
      _ee.erase(_ee.begin()+2);
    else if(score02<score01 && score02<score12)
      _ee.erase(_ee.begin()+1);
    else
      _ee.erase(_ee.begin()+0);
  }
  ASSERT_MSGV(_ee.size()==4 || _ee.size()==2 || _ee.size()==1,"Simplified dynamics only accept #EE=1/2/4, detected #EE=%ld!",_ee.size())
  //discretize workspace
  detectSymmetry();
  for(sizeType i=0; i<(sizeType)_ee.size(); i++) {
    INFOV("Creating EE#%ld!",_ee[i]._JID[0])
    initializeBB(_ee[i]);
    if(symDir>=0 && symDir<3)
      cullSymmetry(i,symDir);
  }
  for(sizeType i=0; i<(sizeType)_ee.size(); i++)
    optimizeBB(_ee[i],zRange,zMargin); //find bb
}
void SimplifiedDynamics::detectEndEffector(const ArticulatedBody& body,sizeType jid,Vec3d& localPos,scalarD& phi0,const Vec3d& zRange)
{
  phi0=0;
  const Joint& joint=body.joint(jid);
  if(joint._spheres.cols()>0) {
    scalarD maxDist=0;
    for(sizeType i=0; i<joint._spheres.cols(); i++)
      if(joint._spheres.col(i).norm()>maxDist) {
        maxDist=joint._spheres.col(i).norm();
        localPos=joint._spheres.col(i);
        phi0=joint._radGeomColl[i];
      }
  }
  if(phi0==0 && zRange.isZero()) {
    ObjMesh m=joint.getMesh(Joint::MESH);
    localPos=m.getCentroid().cast<scalarD>();
    for(sizeType i=0; i<(sizeType)m.getV().size(); i++)
      phi0+=(m.getV(i).cast<scalarD>()-localPos).norm();
    phi0/=(sizeType)m.getV().size();
  }
  if(phi0==0 && !zRange.isZero()) {
    PBDArticulatedGradientInfo<scalarD> info(body,Cold::Zero(body.nrDOF()));
    scalarD minZ=ScalarUtil<scalarD>::scalar_max();
    ObjMesh m=joint.getMesh(Joint::MESH);
    for(sizeType i=0; i<(sizeType)m.getV().size(); i++) {
      Vec3d v=ROTI(info._TM,jid)*m.getV(i).cast<scalarD>()+CTRI(info._TM,jid);
      minZ=std::min<scalarD>(minZ,v.dot(zRange.normalized()));
    }
    std::set<sizeType> footV;
    for(sizeType i=0; i<(sizeType)m.getV().size(); i++) {
      Vec3d v=ROTI(info._TM,jid)*m.getV(i).cast<scalarD>()+CTRI(info._TM,jid);
      scalar z=v.dot(zRange.normalized());
      if(z<minZ+zRange.norm())
        footV.insert(i);
    }
    localPos.setZero();
    for(sizeType i:footV)
      localPos+=m.getV(i).cast<scalarD>();
    localPos/=(sizeType)footV.size();
    for(sizeType i:footV)
      phi0+=(m.getV(i).cast<scalarD>()-localPos).norm();
    phi0/=(sizeType)footV.size();
  }
  ASSERT_MSG(phi0>0,"Failed detecting end-effector!")
}
Cold SimplifiedDynamics::inverseKinematics
(const ArticulatedBody& body,const std::vector<EndEffectorBounds>& ee,
 const std::vector<Vec3d,Eigen::aligned_allocator<Vec3d>>& pos,
 const Mat3d& H,const Cold* init,bool rotBase,bool transBase)
{
#ifdef OPTIMIZER_SUPPORT
  Cold x;
  bool succ=true;
  DSSQPObjectiveCompound<scalarD> obj;
  //constraint
  for(sizeType i=0; i<(sizeType)ee.size(); i++) {
    std::shared_ptr<InverseKinematicsObj<scalarD>> comp(new InverseKinematicsObj<scalarD>(obj,body,ee[i],0,true,init));
    comp->_target=pos[i];
    comp->_H=H;
    obj.addComponent(comp);
  }
  //objective
  if(init) {
    std::shared_ptr<VariationalSimilarityObj<scalarD>> comp(new VariationalSimilarityObj<scalarD>(obj,body,*init,0,rotBase,transBase));
    obj.addComponent(comp);
  }
  //solve
  if(init)
    x=concat(*init,init->segment(0,obj.inputs()-init->size()));
  else x.setZero(obj.inputs());
#if defined(KNITRO_SUPPORT)
  KNInterface<scalarD>(obj).solve(false,0,1e-8f,1e-8f,true,1000,-1,x);
#elif defined(IPOPT_SUPPORT)
  IPOPTInterface<scalarD>::optimize(x,obj,1e-8f,1000,0,-1,true,0);
#else
  FUNCTION_NOT_IMPLEMENTED
#endif
  PBDArticulatedGradientInfo<scalarD> info(body,x);
  body.writeVTK(info._TM,"test.vtk",Joint::MESH);
  //succ
  x=x.segment(0,body.nrDOF()).eval();
  for(sizeType i=0; i<(sizeType)ee.size(); i++) {
    INFOV("EE%ld-Error: %f",i,(H*(ee[i].globalPos(body,x)-pos[i])).cwiseAbs().maxCoeff())
    if((H*(ee[i].globalPos(body,x)-pos[i])).cwiseAbs().maxCoeff()>1e-4f)
      succ=false;
  }
  return succ?x:Cold::Zero(0);
#else
  FUNCTION_NOT_IMPLEMENTED
  return Cold::Zero(0);
#endif
}
Cold SimplifiedDynamics::inverseKinematics(const EndEffectorBounds& ee,const Vec3d& pos,const Cold* init) const
{
#ifdef OPTIMIZER_SUPPORT
  DSSQPObjectiveCompound<scalarD> obj;
  std::shared_ptr<InverseKinematicsObj<scalarD>> comp(new InverseKinematicsObj<scalarD>(obj,*_body,ee,0));
  comp->_target=pos;
  obj.addComponent(comp);
  //solve
  Cold x=(init&&init->size()==obj.inputs())?*init:Cold::Zero(obj.inputs());
#if defined(KNITRO_SUPPORT)
  KNInterface<scalarD>(obj).solve(false,0,1e-8f,1e-8f,true,1000,-1,x);
#elif defined(IPOPT_SUPPORT)
  IPOPTInterface<scalarD>::optimize(x,obj,1e-8f,1000,0,-1,true,0);
#else
  FUNCTION_NOT_IMPLEMENTED
#endif
  if((ee.globalPos(*_body,x)-pos).cwiseAbs().maxCoeff()<1e-4f)
    return x;
  else return Cold::Zero(0);
#else
  FUNCTION_NOT_IMPLEMENTED
  return Cold::Zero(0);
#endif
}
PBDArticulatedGradientInfo<scalarD> SimplifiedDynamics::getDOF(const Vec6d& base,const Mat3Xd& xss) const
{
#ifdef OPTIMIZER_SUPPORT
  ASSERT_MSG(xss.cols()==(sizeType)_ee.size(),"End Effector size mismatch in getDOF!")
  Cold x=Cold::Zero(_body->nrDOF());
  Mat3d R=expWGradV<scalarD,Vec3d>(base.segment<3>(3));
  for(sizeType i=0; i<(sizeType)_ee.size(); i++) {
    Vec3d xssi=xss.col(i);
    xssi=R.transpose()*(xssi-base.segment<3>(0));
    Cold init=_ee[i].getDOF(xssi,_ctr);
    x+=InverseKinematicsObj<scalarD>::assignDOF(*_body,_ee[i],init,NULL);
  }
  x.segment(0,6)=base;
  return PBDArticulatedGradientInfo<scalarD>(*_body,x);
#else
  FUNCTION_NOT_IMPLEMENTED
  return PBDArticulatedGradientInfo<scalarD>();
#endif
}
PBDArticulatedGradientInfo<scalarD> SimplifiedDynamics::getDOF(const Mat3X4d& base,const Mat3Xd& xss) const
{
  return getDOF((Vec6d)concat<Cold>(CTR(base),invExpW<scalarD>(ROT(base))),xss);
}
void SimplifiedDynamics::debugDOF(const std::string& path,sizeType nrPose) const
{
  Mat3Xd xss=Mat3Xd::Zero(3,_ee.size());
  Mat3Xd xss1=Mat3Xd::Zero(3,_ee.size());
  Mat3Xd xss2=Mat3Xd::Zero(3,_ee.size());
  for(sizeType i=0; i<nrPose; i++) {
    Vec6d base=Vec6d::Random()*M_PI;
    Mat3d R=expWGradV<scalarD,Vec3d>(base.segment<3>(3));
    for(sizeType d=0; d<(sizeType)_ee.size(); d++) {
      xss1.col(d)=_ee[d].randomPos(_ctr);
      xss.col(d)=R*xss1.col(d)+base.segment<3>(0);
    }
    PBDArticulatedGradientInfo<scalarD> info=getDOF(base,xss);
    _body->writeVTK(info._TM,path+"/pose"+std::to_string(i)+".vtk",Joint::MESH);
    for(sizeType d=0; d<(sizeType)_ee.size(); d++)
      xss2.col(d)=_ee[d].globalPosAll(*_body,info._xM);
    INFOV("poseErr: %lf!",(xss1-xss2).norm())
  }
}
void SimplifiedDynamics::optimizeBB(EndEffectorBounds& ee,const Vec3d& zRange,scalarD zMargin) const
{
  Vec3d pos=ee.globalPos(*_body,Cold::Zero(ee.nrDOF(*_body)))-_ctr;
  std::unordered_set<std::pair<Vec3i,Vec3i>,Hash> exists,existsExpanded;
  for(const std::pair<Vec3i,Cold>& id:ee._indices) {
    Vec3d pos=ee.globalPos(*_body,Cold::Zero(ee.nrDOF(*_body)));
    if(zRange.norm()>0 && zMargin>0 && (id.first.cast<scalarD>()*ee._res+_ctr).dot(zRange.normalized())>pos.dot(zRange.normalized())+ee._phi0*zMargin)
      continue;
    BBox<sizeType> bb(id.first,id.first+Vec3i::Ones());
    if(!ee.sameOrthant(bb._minC,pos,_sym) || !ee.sameOrthant(bb._maxC,pos,_sym))
      continue;
    if(ee.exists(bb))
      exists.insert(std::make_pair(bb._minC,bb._maxC));
  }
  while(true) {
    existsExpanded.clear();
    for(const std::pair<Vec3i,Vec3i>& bb:exists)
      for(sizeType d=0; d<3; d++) {
        BBox<sizeType> bbd(bb.first,bb.second);
        bbd._minC[d]=bbd._maxC[d]=bb.second[d]+1;
        if(ee.exists(bbd)) {
          bbd._minC[d]=bb.first[d];
          if(!ee.sameOrthant(bbd._minC,pos,_sym) || !ee.sameOrthant(bbd._maxC,pos,_sym))
            continue;
          existsExpanded.insert(std::make_pair(bbd._minC,bbd._maxC));
        }
      }
    INFOV("DynamicProgrammingPass: %ld",existsExpanded.size())
    if(existsExpanded.empty())
      break;
    exists.swap(existsExpanded);
  }
  bool found=false;
  for(const std::pair<Vec3i,Vec3i>& bb:exists) {
    if(!found || ee.isBetterBB(BBox<sizeType>(bb.first,bb.second),_sym))
      ee._bb=BBox<sizeType>(bb.first,bb.second);
    found=true;
  }
  INFOV("Found BBSize=%ldx%ldx%ld=%ld",ee._bb.getExtent()[0],ee._bb.getExtent()[1],ee._bb.getExtent()[2],ee._bb.getExtent().prod())
}
void SimplifiedDynamics::initializeBB(EndEffectorBounds& ee) const
{
  disableTiming();
  ee._indices.clear();
  sizeType initRange=3;
  std::vector<std::pair<Vec3i,Cold>> ss,ss2;
  std::unordered_set<Vec3i,Hash> invalids;
  Vec3d pos=ee.globalPos(*_body,Cold::Zero(ee.nrDOF(*_body)));
  for(sizeType x=-initRange; x<=initRange; x++)
    for(sizeType y=-initRange; y<=initRange; y++)
      for(sizeType z=-initRange; z<=initRange; z++)
        ss.push_back(std::make_pair(Vec3i(x,y,z)+((pos-_ctr)/ee._res).cast<sizeType>(),Cold::Zero(0)));
  while(!ss.empty()) {
    //pass-1 parallel inverse kinematics
    INFOV("InitializingPass: %ld",ss.size())
    std::vector<Cold,Eigen::aligned_allocator<Cold>> xss(ss.size());
    //OMP_PARALLEL_FOR_
    for(sizeType i=0; i<(sizeType)ss.size(); i++) {
      Vec3i id=ss[i].first;
      Cold init=ss[i].second;
      xss[i]=inverseKinematics(ee,id.template cast<scalarD>()*ee._res+_ctr,&init);
    }
    //pass-2 expand
    ss2.clear();
    for(sizeType i=0; i<(sizeType)ss.size(); i++)
      if(xss[i].size()==0)
        invalids.insert(ss[i].first);
      else {
        ee._indices[ss[i].first]=xss[i];
        for(sizeType x=-1; x<=1; x++)
          for(sizeType y=-1; y<=1; y++)
            for(sizeType z=-1; z<=1; z++) {
              Vec3i id2=Vec3i(x,y,z)+ss[i].first;
              if(ee._indices.find(id2)==ee._indices.end() && invalids.find(id2)==invalids.end())
                ss2.push_back(std::make_pair(id2,Cold::Zero(0)));
            }
      }
    std::unordered_map<Vec3i,Cold,Hash> ss2Set(ss2.begin(),ss2.end());
    ss.assign(ss2Set.begin(),ss2Set.end());
  }
}
void SimplifiedDynamics::detectSymmetry()
{
  _ctr.setZero();
  _sym.clear();
  if((sizeType)_ee.size()==4) {
    Vec3d p[4];
    sizeType index;
    for(sizeType d=0; d<4; d++)
      p[d]=_ee[d].globalPos(*_body,Cold::Zero(_ee[d].nrDOF(*_body)));
    for(sizeType d=0; d<4; d++)
      for(sizeType d2=0; d2<d; d2++) {
        (p[d]-p[d2]).cwiseAbs().maxCoeff(&index);
        if(std::find(_sym.begin(),_sym.end(),index)==_sym.end())
          _sym.push_back(index);
      }
    _ctr=(p[0]+p[1]+p[2]+p[3])/4;
    ASSERT_MSG(_sym.size()==2,"Four-foot symmetry detection failed!")
  } else if((sizeType)_ee.size()==2) {
    sizeType index;
    Vec3d p0=_ee[0].globalPos(*_body,Cold::Zero(_ee[0].nrDOF(*_body)));
    Vec3d p1=_ee[1].globalPos(*_body,Cold::Zero(_ee[1].nrDOF(*_body)));
    (p1-p0).cwiseAbs().maxCoeff(&index);
    _sym.push_back(index);
    _ctr=(p0+p1)/2;
  } else {
    ASSERT((sizeType)_ee.size()==1)
  }
}
void SimplifiedDynamics::cullSymmetry(sizeType symDir)
{
  for(sizeType i=0; i<(sizeType)_ee.size(); i++)
    if(symDir>=0 && symDir<3)
      cullSymmetry(i,symDir);
}
void SimplifiedDynamics::cullSymmetry(sizeType i,sizeType symDir)
{
  //find symmetric pose
  EndEffectorBounds& eei=_ee[i];
  Vec3d pos=eei.globalPos(*_body,Cold::Zero(eei.nrDOF(*_body)));
  pos[symDir]=_ctr[symDir]*2-pos[symDir];
  //find symmetric ee
  sizeType jSym=-1;
  scalarD minDist=0;
  for(sizeType j=0; j<(sizeType)_ee.size(); j++) {
    Vec3d pos2=_ee[j].globalPos(*_body,Cold::Zero(_ee[j].nrDOF(*_body)));
    if(jSym==-1 || (pos2-pos).norm()<minDist) {
      jSym=j;
      minDist=(pos2-pos).norm();
    }
  }
  ASSERT_MSG(jSym!=i,"Symmetry direction detection failed!")
  INFOV("Symmetric pruning: EE#%ld-EE#%ld!",i,jSym)
  if(jSym<i) {
    std::unordered_set<Vec3i,Hash> delCache;
    EndEffectorBounds& eej=_ee[jSym];
    for(const std::pair<Vec3i,Cold>& id:eei._indices) {
      Vec3i other=id.first;
      other[symDir]*=-1;
      if(eej._indices.find(other)==eej._indices.end())
        delCache.insert(id.first);
    }
    INFOV("Culling %ld indices!",delCache.size())
    for(const Vec3i& del:delCache)
      eei._indices.erase(del);
  }
  if(jSym<i) {
    std::unordered_set<Vec3i,Hash> delCache;
    EndEffectorBounds& eej=_ee[jSym];
    for(const std::pair<Vec3i,Cold>& id:eej._indices) {
      Vec3i other=id.first;
      other[symDir]*=-1;
      if(eei._indices.find(other)==eei._indices.end())
        delCache.insert(id.first);
    }
    INFOV("Culling %ld indices!",delCache.size())
    for(const Vec3i& del:delCache)
      eej._indices.erase(del);
  }
}
void SimplifiedDynamics::writeVTK(const std::string& path) const
{
  for(sizeType i=0; i<(sizeType)_ee.size(); i++) {
    std::vector<Vec3d,Eigen::aligned_allocator<Vec3d>> bbs;
    std::vector<Vec3d,Eigen::aligned_allocator<Vec3d>> pts;
    bbs.push_back(_ee[i]._bb._minC.template cast<scalarD>()*_ee[i]._res+_ctr);
    bbs.push_back(_ee[i]._bb._maxC.template cast<scalarD>()*_ee[i]._res+_ctr);
    for(const std::pair<Vec3i,Cold>& id:_ee[i]._indices)
      pts.push_back(id.first.cast<scalarD>()*_ee[i]._res+_ctr);
    //VTKWriter<scalarD>
    VTKWriter<scalarD> os("SimplifiedBB",path+"/EE"+std::to_string(_ee[i].jointId())+".vtk",true);
    os.appendVoxels(bbs.begin(),bbs.end(),true);
    os.setRelativeIndex();
    os.appendPoints(pts.begin(),pts.end());
    os.appendCells(VTKWriter<scalarD>::IteratorIndex<Vec3i>(0,0,0),
                   VTKWriter<scalarD>::IteratorIndex<Vec3i>((sizeType)pts.size(),0,0),
                   VTKWriter<scalarD>::POINT,true);
  }
}
std::shared_ptr<ArticulatedBody> SimplifiedDynamics::getBody() const
{
  return _body;
}
EndEffectorBounds SimplifiedDynamics::getForceInfo(sizeType jid) const
{
  for(sizeType i=0; i<(sizeType)_ee.size(); i++)
    if(_ee[i]._JID[0]==jid)
      return _ee[i];
  ASSERT_MSGV(false,"Cannot find end-effector for joint%ld",jid)
  return EndEffectorBounds();
}
const std::vector<EndEffectorBounds>& SimplifiedDynamics::getForceInfos() const
{
  return _ee;
}
Vec3d SimplifiedDynamics::lb(sizeType eeid) const
{
  const EndEffectorBounds& ee=_ee[eeid];
  return ee._bb._minC.template cast<scalarD>()*ee._res+_ctr;
}
Vec3d SimplifiedDynamics::ub(sizeType eeid) const
{
  const EndEffectorBounds& ee=_ee[eeid];
  return ee._bb._maxC.template cast<scalarD>()*ee._res+_ctr;
}
bool SimplifiedDynamics::read(std::istream& is,IOData* dat)
{
  registerType<ArticulatedBody>(dat);
  registerType<EndEffectorBounds>(dat);
  readBinaryData(_body,is,dat);
  readBinaryData(_ee,is,dat);
  readBinaryData(_sym,is,dat);
  readBinaryData(_rootId,is,dat);
  readBinaryData(_ctr,is,dat);
  return is.good();
}
bool SimplifiedDynamics::write(std::ostream& os,IOData* dat) const
{
  registerType<ArticulatedBody>(dat);
  registerType<EndEffectorBounds>(dat);
  writeBinaryData(_body,os,dat);
  writeBinaryData(_ee,os,dat);
  writeBinaryData(_sym,os,dat);
  writeBinaryData(_rootId,os,dat);
  writeBinaryData(_ctr,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> SimplifiedDynamics::copy() const
{
  return std::shared_ptr<SerializableBase>(new SimplifiedDynamics);
}
std::string SimplifiedDynamics::type() const
{
  return typeid(SimplifiedDynamics).name();
}
//instance
#ifdef OPTIMIZER_SUPPORT
PRJ_BEGIN
template class InverseKinematicsObj<double>;
#ifdef ALL_TYPES
template class InverseKinematicsObj<__float128>;
template class InverseKinematicsObj<mpfr::mpreal>;
#endif
PRJ_END
#endif
