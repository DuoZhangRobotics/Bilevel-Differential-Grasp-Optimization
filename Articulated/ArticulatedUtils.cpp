#include "ArticulatedUtils.h"
#include "PBDArticulatedGradientInfo.h"
#include "Joint.h"
#include "JointFunc.h"
#include <CommonFile/Interp.h>
#include <CommonFile/geom/StaticGeomCell.h>
#include <Utils/DebugGradient.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

ArticulatedUtils::ArticulatedUtils(ArticulatedBody& body):_body(body) {}
void ArticulatedUtils::assembleGlobalVars(const tinyxml2::XMLElement& pt)
{
  _globalRho=get<scalarF>(pt,"globalRho",1);
  _globalLimitCoef=get<scalarF>(pt,"globalLimitCoef",1000);
  _globalControlCoef=get<scalarF>(pt,"globalControlCoef",1);
  _globalDampingCoef=get<scalarF>(pt,"globalDampingCoef",0);
  _globalSampleDist=get<scalarF>(pt,"globalSampleDist",ScalarUtil<scalarF>::scalar_max());
  _globalSampleDistLocal=get<scalarF>(pt,"globalSampleDistLocal",1);
  _globalRes=get<sizeType>(pt,"globalRes",8);
}
void ArticulatedUtils::assemble(const tinyxml2::XMLElement& pt)
{
  assembleGlobalVars(pt);
  //build joints
  std::vector<TransInfo> infos;
  assembleJoints(*getChild(pt,"joint"),infos);
  //apply X0
  for(sizeType i=0; i<(sizeType)_body.nrJ(); i++) {
    Joint& J=_body.joint(i);
    const ObjMesh mesh=J.getMesh(Joint::MESH);
    if(!mesh.getV().empty())
      J.transformMesh(infos[i]._X0);
    if(J._parent >= 0) {
      const TransInfo& infoP=infos[J._parent];
      APPLY_TRANS(J._trans,infoP._X0,J._trans)
    }
  }
  //assemble joints
  sizeType offDOF=0,offDDT=0;
  _body._geom.reset(new StaticGeom(3));
  for(sizeType i=0; i<(sizeType)_body.nrJ(); i++) {
    const TransInfo& I=infos[i];
    Joint& J=_body.joint(i);
    J._offDOF=offDOF;
    J._offDDT=offDDT;
    offDOF+=J.nrDOF();
    offDDT+=J.nrDDT();
    J.assemble(I._rho);
    if(J._mesh)
      _body._geom->addGeomCell(std::dynamic_pointer_cast<StaticGeomCell>(J._mesh->copy()));
  }
  _body._geom->assemble();
  _body.fillChildren();
}
void ArticulatedUtils::sampleGeom3DBox(Joint& J) const
{
  J._radGeomColl.setZero(1);
  J._radSelfColl.setZero(1);
  const ObjMesh mesh=J.getMesh(Joint::MESH);
  J._spheres.resize(3,mesh.getV().size());
  for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++)
    J._spheres.col(i)=mesh.getV()[i].cast<scalarD>();
  J._radGeomColl.setZero(J._spheres.cols());
  J._radSelfColl.setZero(J._spheres.cols());
}
void ArticulatedUtils::sampleGeom2Sphere(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const
{
  std::shared_ptr<TwoSphereMeshCell> twoSphere=std::dynamic_pointer_cast<TwoSphereMeshCell>(J._mesh);
  Vec3 ctr1=twoSphere->ctr1();
  Vec3 ctr2=twoSphere->ctr2();
  scalar minRad=std::min(I._rad,I._rad2);
  scalar maxRad=std::max(I._rad,I._rad2);
  scalarD dist=get<scalarF>(pt,"sampleDist",_globalSampleDist);
  dist=std::min<scalarD>(dist,minRad*get<scalarF>(pt,"sampleDistLocal",_globalSampleDistLocal));
  sizeType nr=std::max<sizeType>(ceil(I._lenY.norm()/dist),1)+1;
  J._spheres.setZero(3,nr);
  J._radGeomColl.setZero(nr);
  J._radSelfColl.setZero(nr);
  for(sizeType y=0; y<nr; y++) {
    scalar alpha=y/scalarD(nr-1);
    J._spheres.col(y)=interp1D<Vec3>(ctr1,ctr2,alpha).cast<scalarD>();
    J._radGeomColl[y]=interp1D(I._rad,I._rad2,alpha);
    J._radSelfColl[y]=Vec2d(maxRad,dist/2).norm();
  }
}
void ArticulatedUtils::sampleGeom3Sphere(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const
{
  std::shared_ptr<ThreeSphereMeshCell> threeSphere=std::dynamic_pointer_cast<ThreeSphereMeshCell>(J._mesh);
  Vec3 ctr1=threeSphere->ctr1();
  Vec3 ctr2=threeSphere->ctr2();
  Vec3 ctr3=threeSphere->ctr3();
  scalar minRad=std::min(std::min(I._rad,I._rad2),I._rad3);
  scalar maxRad=std::max(std::max(I._rad,I._rad2),I._rad3);
  scalarD dist=get<scalarF>(pt,"sampleDist",_globalSampleDist);
  dist=std::min<scalarD>(dist,minRad*get<scalarF>(pt,"sampleDistLocal",_globalSampleDistLocal));
  //count subd
  sizeType nrSubd=0;
  scalarD maxDist=std::max(std::max(I._lenX.norm(),I._lenY.norm()),I._lenZ.norm());
  while(maxDist > dist) {
    maxDist/=2;
    nrSubd++;
  }
  //create sample point
  ObjMesh m;
  m.getV().push_back(ctr1);
  m.getV().push_back(ctr2);
  m.getV().push_back(ctr3);
  m.getI().push_back(Vec3i(0,1,2));
  m.smooth();
  m.subdivide(nrSubd);
  sizeType nr=(sizeType)m.getV().size();
  J._spheres.setZero(3,nr);
  J._radGeomColl.setZero(nr);
  J._radSelfColl.setZero(nr);
  TriangleTpl<scalar> tri(ctr1,ctr2,ctr3);
  for(sizeType i=0; i<nr; i++) {
    Vec3 bary=tri.bary(m.getV(i));
    J._spheres.col(i)=m.getV(i).cast<scalarD>();
    J._radGeomColl[i]=bary.dot(Vec3(I._rad,I._rad2,I._rad3));
    J._radSelfColl[i]=Vec2d(maxRad,maxDist/2).norm();
  }
}
void ArticulatedUtils::sampleGeom2D(Joint& J,const Vec3d& base,const Vec3d& lenX,const Vec3d& lenY,scalarD rad,const tinyxml2::XMLElement& pt) const
{
  scalarD dist=get<scalarF>(pt,"sampleDist",_globalSampleDist);
  dist=std::min(dist,rad*get<scalarF>(pt,"sampleDistLocal",_globalSampleDistLocal));
  sizeType nrX=std::max<sizeType>(ceil(lenX.norm()/dist),1)+1;
  sizeType nrY=std::max<sizeType>(ceil(lenY.norm()/dist),1)+1;
  scalarD distX=lenX.norm()/scalarD(nrX-1);
  scalarD distY=lenY.norm()/scalarD(nrY-1);
  sizeType i=J._spheres.cols();
  J._spheres=concatCol(J._spheres,Mat3Xd::Zero(3,nrX*nrY));
  for(sizeType x=0; x<nrX; x++)
    for(sizeType y=0; y<nrY; y++)
      J._spheres.col(i++)=lenX.normalized()*scalarD(x)*distX+lenY.normalized()*scalarD(y)*distY+base;
  J._radGeomColl.setConstant(i,rad);
  J._radSelfColl.setConstant(i,Vec3d(rad,distX/2,distY/2).norm());
}
void ArticulatedUtils::sampleGeom3D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const
{
  J._spheres.setZero(3,0);
  sampleGeom2D(J,-I._lenX/2-I._lenZ/2,        I._lenX,I._lenZ,I._rad,pt);
  sampleGeom2D(J,-I._lenX/2-I._lenZ/2+I._lenY,I._lenX,I._lenZ,I._rad,pt);

  sampleGeom2D(J,-I._lenX/2-I._lenZ/2,I._lenX,I._lenY,I._rad,pt);
  sampleGeom2D(J,-I._lenX/2+I._lenZ/2,I._lenX,I._lenY,I._rad,pt);

  sampleGeom2D(J,-I._lenZ/2-I._lenX/2,I._lenZ,I._lenY,I._rad,pt);
  sampleGeom2D(J,-I._lenZ/2+I._lenX/2,I._lenZ,I._lenY,I._rad,pt);
}
void ArticulatedUtils::sampleGeom2D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const
{
  J._spheres.setZero(3,0);
  sampleGeom2D(J,-I._lenX/2,I._lenX,I._lenY,I._rad,pt);
}
void ArticulatedUtils::sampleGeom1D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const
{
  scalarD dist=get<scalarF>(pt,"sampleDist",_globalSampleDist);
  dist=std::min<scalarD>(dist,I._rad*get<scalarF>(pt,"sampleDistLocal",_globalSampleDistLocal));
  sizeType nr=std::max<sizeType>(ceil(I._lenY.norm()/dist),1)+1;
  dist=I._lenY.norm()/scalarD(nr-1);
  J._spheres.setZero(3,nr);
  for(sizeType y=0; y<nr; y++)
    J._spheres.col(y)=I._lenY.normalized()*scalarD(y)*dist;
  J._radGeomColl.setConstant(nr,I._rad);
  J._radSelfColl.setConstant(nr,Vec2d(I._rad,dist/2).norm());
}
void ArticulatedUtils::sampleGeom0D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement&) const
{
  J._spheres.setZero(3,1);
  J._radGeomColl.setConstant(1,I._rad);
  J._radSelfColl.setConstant(1,I._rad);
}
void ArticulatedUtils::assembleJoints(const tinyxml2::XMLElement& pt,std::vector<TransInfo>& infos)
{
  std::stack<std::pair<const tinyxml2::XMLElement*,sizeType> > ss;
  ss.push(std::make_pair(&pt,-1));
  while(!ss.empty()) {
    std::pair<const tinyxml2::XMLElement*,sizeType> p=ss.top();
    ss.pop();
    //assemble this joint
    assembleJoint(*(p.first),p.second,infos);
    //find children
    for(const tinyxml2::XMLElement* v=p.first->FirstChildElement(); v; v=v->NextSiblingElement()) {
      if(std::string(v->Name()) == "joint")
        ss.push(std::make_pair(v,_body.nrJ()-1));
    }
  }
}
void ArticulatedUtils::assembleJoint(const tinyxml2::XMLElement& pt,sizeType parent,std::vector<TransInfo>& infos)
{
  Joint J;
  Mat3d R;
  Mat3X4d JT;
  Mat3Xd DT,DDT;
  TransInfo info;
  std::string defaultName;
  TransInfo infoParent;
  if(parent >= 0)
    infoParent=infos[parent];
  else {
    infoParent._TLast.setIdentity();
    infoParent._X0.setIdentity();
    infoParent._lenX.setZero();
    infoParent._lenY.setZero();
    infoParent._lenZ.setZero();
    infoParent._origin.setZero();
  }
  //basics
  if(parent < 0)
    defaultName="joint";
  else defaultName=_body._joints[parent]._name+"_"+std::to_string(_body._joints.size());
  J._color=parsePtreeDef<Vec3>(pt,"color",_globalColor);
  J._name=get<std::string>(pt,"name",defaultName);
  J.setType(get<std::string>(pt,"type"));
  J._parent=parent;
  J._depth=parent<0?0:_body._joints[parent]._depth+1;
  J._offDDT=J._offDOF=0;
  J._mult=J._offset=0;
  DT.resize(3,J.nrDOF());
  DDT.resize(3,J.nrDDT());
  Eigen::Map<Mat3Xd,0,Eigen::OuterStride<> > DTMap(DT.data(),DT.rows(),DT.cols(),Eigen::OuterStride<>(DT.rows()));
  Eigen::Map<Mat3Xd,0,Eigen::OuterStride<> > DDTMap(DDT.data(),DDT.rows(),DDT.cols(),Eigen::OuterStride<>(DDT.rows()));
  //geometry
  info._lenX.setZero();
  info._lenY.setZero();
  info._lenZ.setZero();
  std::vector<Joint> geomSpheres;
  std::vector<std::shared_ptr<StaticGeomCell> > geoms;
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(!beginsWith(std::string(v->Name()),"geom"))
      continue;
    const tinyxml2::XMLElement& c=*v;
    if(std::string(v->Name()) == "geom3D") {
      info._rad=get<scalarF>(c,"rad");
      info._lenX[0]=get<scalarF>(c,"lenX");
      info._lenY[1]=get<scalarF>(c,"lenY");
      info._lenZ[2]=get<scalarF>(c,"lenZ");
      if(info._rad == 0)
        J._mesh.reset(new BoxGeomCell(Mat4::Identity(),3,Vec3(info._lenX.norm()/2,info._lenY.norm()/2,info._lenZ.norm()/2)));
      else J._mesh.reset(new SphericalBoxGeomCell(Mat4::Identity(),3,Vec4(info._lenX.norm()/2,info._lenY.norm()/2,info._lenZ.norm()/2,info._rad)));
      J.transformMesh(Vec3d(0,info._lenY.norm()/2,0));
      if(info._rad == 0)
        sampleGeom3DBox(J);
      else sampleGeom3D(J,info,pt);
    } else if(std::string(v->Name()) == "geom2D") {
      info._rad=get<scalarF>(c,"rad");
      info._lenX[0]=get<scalarF>(c,"lenX");
      info._lenY[1]=get<scalarF>(c,"lenY");
      info._lenZ[2]=info._rad;
      ASSERT(info._rad>0)
      J._mesh.reset(new SphericalBoxGeomCell(Mat4::Identity(),3,Vec4(info._lenX.norm()/2,info._lenY.norm()/2,0,info._rad)));
      J.transformMesh(Vec3d(0,info._lenY.norm()/2,0));
      sampleGeom2D(J,info,pt);
    } else if(std::string(v->Name()) == "geom1D") {
      info._rad=get<scalarF>(c,"rad");
      info._lenX[0]=info._rad;
      info._lenY[1]=get<scalarF>(c,"lenY");
      info._lenZ[2]=info._rad;
      ASSERT(info._rad>0)
      J._mesh.reset(new CapsuleGeomCell(Mat4::Identity(),3,info._rad,info._lenY.norm()/2));
      J.transformMesh(Vec3d(0,info._lenY.norm()/2,0));
      sampleGeom1D(J,info,pt);
    } else if(std::string(v->Name()) == "geom0D") {
      info._rad=get<scalarF>(c,"rad");
      info._lenX[0]=info._rad;
      info._lenY[1]=info._rad;
      info._lenZ[2]=info._rad;
      ASSERT(info._rad>0)
      J._mesh.reset(new SphereGeomCell(Mat4::Identity(),3,info._rad));
      sampleGeom0D(J,info,pt);
    } else if(std::string(v->Name()) == "geom2Sphere") {
      info._rad=get<scalarF>(c,"rad1");
      info._rad2=get<scalarF>(c,"rad2");
      Vec3 ctr1=parsePtree<Vec3>(c,"ctr1");
      Vec3 ctr2=parsePtree<Vec3>(c,"ctr2");
      info._lenX.setZero();
      info._lenY=(ctr2-ctr1).cast<scalarD>();
      info._lenZ.setZero();
      J._mesh.reset(new TwoSphereMeshCell(3,ctr1,ctr2,info._rad,info._rad2));
      sampleGeom2Sphere(J,info,pt);
    } else if(std::string(v->Name()) == "geom3Sphere") {
      info._rad=get<scalarF>(c,"rad1");
      info._rad2=get<scalarF>(c,"rad2");
      info._rad3=get<scalarF>(c,"rad3");
      Vec3 ctr1=parsePtree<Vec3>(c,"ctr1");
      Vec3 ctr2=parsePtree<Vec3>(c,"ctr2");
      Vec3 ctr3=parsePtree<Vec3>(c,"ctr3");
      info._lenX=(ctr2-ctr1).cast<scalarD>();
      info._lenY=(ctr3-ctr1).cast<scalarD>();
      info._lenZ=(ctr3-ctr2).cast<scalarD>();
      J._mesh.reset(new ThreeSphereMeshCell(3,ctr1,ctr2,ctr3,info._rad,info._rad2,info._rad3));
      sampleGeom3Sphere(J,info,pt);
    }
    if(J._mesh)
      J._mesh->setRes(get<sizeType>(c,"res",_globalRes));
    geomSpheres.push_back(J);
    geoms.push_back(J._mesh);
  }
  if(geoms.size()>1) {
    J._spheres.resize(3,0);
    J._radGeomColl.resize(0);
    J._radSelfColl.resize(0);
    for(sizeType i=0; i<(sizeType)geomSpheres.size(); i++) {
      J._spheres=concatCol(J._spheres,geomSpheres[i]._spheres);
      J._radGeomColl=concat(J._radGeomColl,geomSpheres[i]._radGeomColl);
      J._radSelfColl=concat(J._radSelfColl,geomSpheres[i]._radSelfColl);
    }
    J._mesh.reset(new CompositeGeomCell(Mat4::Identity(),geoms));
  }
  //transformation: translation
  CTR(J._trans)=infoParent._origin;
  CTR(J._trans)+=infoParent._lenX*get<scalarF>(pt,"transX",0);
  CTR(J._trans)+=infoParent._lenY*get<scalarF>(pt,"transY",0);
  CTR(J._trans)+=infoParent._lenZ*get<scalarF>(pt,"transZ",0);
  CTR(J._trans)+=parsePtreeDef<Vec3d>(pt,"trans",0);
  //transformation: rotation-joint
  if(hasAttribute(pt,"rotQuat")) {
    Vec4d r=parsePtree<Vec4d>(pt,"rotQuat");
    R=Quatd(r[0],r[1],r[2],r[3]).toRotationMatrix();
  } else if(J._typeJoint == Joint::ROT_3D_EXP || J._typeJoint == Joint::ROT_3D_XYZ) {
    R.col(1)=parsePtreeDef<Vec3d>(pt,"jointDirY",Vec3d::UnitY()).normalized();
    R.col(2)=parsePtreeDef<Vec3d>(pt,"jointDirZ",Vec3d::UnitZ()).normalized();
    R.col(0)=R.col(1).cross(R.col(2)).normalized();
    R.col(2)=R.col(0).cross(R.col(1)).normalized();
  } else if(J._typeJoint == Joint::BALL_JOINT) {
    if(!hasAttribute(pt,"jointDirY")) {
      Vec3d Zto=parsePtreeDef<Vec3d>(pt,"jointDirZ",Vec3d::UnitZ()).normalized();
      R=Quatd::FromTwoVectors(Vec3d::UnitZ(),Zto).toRotationMatrix();
    } else {
      R.col(1)=parsePtree<Vec3d>(pt,"jointDirY").normalized();
      R.col(2)=parsePtree<Vec3d>(pt,"jointDirZ").normalized();
      R.col(0)=R.col(1).cross(R.col(2)).normalized();
      R.col(2)=R.col(0).cross(R.col(1)).normalized();
    }
  } else if(J._typeJoint == Joint::HINGE_JOINT) {
    R=Quatd::FromTwoVectors(Vec3d::UnitZ(),parsePtreeDef<Vec3d>(pt,"jointDirZ",Vec3d::UnitZ()).normalized()).toRotationMatrix();
  } else if(J._typeJoint == Joint::TRANS_1D) {
    R=Quatd::FromTwoVectors(Vec3d::UnitX(),parsePtreeDef<Vec3d>(pt,"jointDirX",Vec3d::UnitX()).normalized()).toRotationMatrix();
  } else {
    R=ROT(infoParent._TLast);
  }
  ROT(J._trans)=ROT(infoParent._TLast).transpose()*R;
  if(hasAttribute(pt,"rotQuatAbs")) {
    Vec4d r=parsePtree<Vec4d>(pt,"rotQuatAbs");
    ROT(J._trans)=Quatd(r[0],r[1],r[2],r[3]).toRotationMatrix();
  }
  //transformation: rotation-mesh
  R.setIdentity();
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(!beginsWith(std::string(v->Name()),"geom"))
      continue;
    if(std::string(v->Name()) == "geom2D" || std::string(v->Name()) == "geom3D") {
      R.col(0)=parsePtree<Vec3d>(pt,"meshDirX").normalized();
      R.col(1)=parsePtree<Vec3d>(pt,"meshDirY").normalized();
      R.col(2)=R.col(0).cross(R.col(1)).normalized();
      R.col(1)=R.col(2).cross(R.col(0)).normalized();
    } else if(std::string(v->Name()) == "geom1D") {
      R=Quatf::FromTwoVectors(Vec3f::UnitY(),parsePtree<Vec3f>(pt,"meshDirY").normalized()).toRotationMatrix().cast<scalarD>();
    } else if(std::string(v->Name()) == "geom0D" || std::string(v->Name()) == "geom2Sphere" || std::string(v->Name()) == "geom3Sphere") {
      R.setIdentity();
    }
  }
  R=((ROT(infoParent._TLast)*ROT(J._trans)).transpose()*R).eval();
  if(hasAttribute(pt,"rotMeshQuatAbs")) {
    Vec4d r=parsePtree<Vec4d>(pt,"rotMeshQuatAbs");
    R=Quatd(r[0],r[1],r[2],r[3]).toRotationMatrix();
  }
  J.transformMesh(R,Vec3d::Zero());
  info._lenX=R*info._lenX;
  info._lenY=R*info._lenY;
  info._lenZ=R*info._lenZ;
  //transformation: translation-mesh
  Vec3d MT=parsePtreeDef<Vec3d>(pt,"meshTrans",Vec3d::Zero());
  info._origin=info._lenX/std::max<scalar>(ScalarUtil<scalar>::scalar_eps(),info._lenX.norm())*MT[0]+
               info._lenY/std::max<scalar>(ScalarUtil<scalar>::scalar_eps(),info._lenY.norm())*MT[1]+
               info._lenZ/std::max<scalar>(ScalarUtil<scalar>::scalar_eps(),info._lenZ.norm())*MT[2];
  J.transformMesh(info._origin);
  //transformation: origin
  if(hasAttribute(pt,"x0")) {
    Cold x0=parsePtree<Cold>(pt,"x0",J.nrDOF());
    Eigen::Map<const Cold> x0Map(x0.data(),x0.size());
    info._X0=JointFunc<scalarD>::TDTDDT(J,x0Map,DTMap,DDTMap);
  } else info._X0.setIdentity();
  //limit/control/damping
  J._limits.setConstant(3,J.nrDOF(),ScalarUtil<scalarD>::scalar_inf());
  J._control.setConstant(J.nrDOF(),0);
  J._damping.setConstant(J.nrDOF(),0);
  if(!J.isRoot(_body) && J.nrDOF()>0) {
    const tinyxml2::XMLElement* c=getChild(pt,"limit");
    if(c) {
      J._limits.row(0)=parsePtree<Cold>(*c,"lower",J.nrDOF());
      J._limits.row(1)=parsePtree<Cold>(*c,"upper",J.nrDOF());
      J._limits.row(2)=parsePtreeDef<Cold>(*c,"coef",_globalLimitCoef,J.nrDOF());
    }
    J._control=parsePtreeDef<Cold>(pt,"controlCoef",_globalControlCoef,J.nrDOF());
    J._damping=parsePtreeDef<Cold>(pt,"dampingCoef",_globalDampingCoef,J.nrDOF());
    for(sizeType i=0; i<J.nrDOF(); i++) {
      if(std::isfinite(J._limits(2,i))) {
        ASSERT_MSGV(J._limits(0,i) <= J._limits(1,i),"Joint: %s, limitRange error: (%f>%f)!",J._name.c_str(),J._limits(0,i),J._limits(1,i))
      }
      ASSERT_MSGV(J._control[i] >= 0,"Joint: %s, controlCoef negative (%f<0)!",J._name.c_str(),J._control[i])
      ASSERT_MSGV(J._damping[i] >= 0,"Joint: %s, dampingCoef negative (%f<0)!",J._name.c_str(),J._damping[i])
    }
  }
  //info
  Cold x0=Cold::Zero(J.nrDOF());
  Eigen::Map<const Cold> x0Map(x0.data(),x0.size());
  JT=JointFunc<scalarD>::TDTDDT(J,x0Map,DTMap,DDTMap);
  APPLY_TRANS(info._TLast,infoParent._TLast,J._trans);
  APPLY_TRANS(info._TLast,info._TLast,JT)
  //descend
  const ObjMesh mesh=J.getMesh(Joint::MESH);
  if(!mesh.getV().empty())
    info._rho=get<scalarF>(pt,"rho",_globalRho);
  _body._joints.push_back(J);
  infos.push_back(info);
}
void ArticulatedUtils::addBase(sizeType dim,const Vec3d& planeNormal)
{
  //check if we already have root
  for(sizeType i=1; i<_body.nrJ(); i++)
    if(_body.joint(i).isRoot(_body)) {
      INFOV("Joint #%ld is already a root, skip addBase!",i)
      return;
    }
  if(_body.joint(0)._typeJoint!=Joint::FIX_JOINT) {
    INFO("Joint #0 is not fixed, skip addBase!")
    return;
  }
  //update index,depth
  _body._joints.insert(_body._joints.begin(),Joint());
  _body._joints.insert(_body._joints.begin(),Joint());
  _body._joints[0]._name="baseTrans";
  _body._joints[1]._name="baseRot";
  for(sizeType i=2; i<_body.nrJ(); i++) {
    Joint& joint=_body.joint(i);
    joint._parent+=2;
    joint._depth+=2;
    if(joint._mimic>=0)
      joint._mimic+=2;
  }
  //type joint
  if(dim==2) {
    _body.joint(1)._typeJoint=Joint::HINGE_JOINT;
    _body.joint(0)._typeJoint=Joint::TRANS_2D;
  } else {
    _body.joint(1)._typeJoint=Joint::ROT_3D_XYZ;
    //_body.joint(1)._typeJoint=Joint::ROT_3D_EXP;
    _body.joint(0)._typeJoint=Joint::TRANS_3D;
  }
  for(sizeType i=0; i<2; i++) {
    Joint& joint=_body.joint(i);
    joint._parent=i-1;
    joint._depth=i;
    //init limit,control,damping
    joint._limits.setConstant(3,joint.nrDOF(),ScalarUtil<scalarD>::scalar_inf());
    joint._control.setConstant(joint.nrDOF(),0);
    joint._damping.setConstant(joint.nrDOF(),0);
    joint._trans.setIdentity();
    //assemble mass (=0)
    joint.assemble(1);
  }
  sizeType offDOF=0,offDDT=0;
  for(sizeType i=0; i<_body.nrJ(); i++) {
    Joint& joint=_body.joint(i);
    joint._offDOF=offDOF;
    joint._offDDT=offDDT;
    offDOF+=joint.nrDOF();
    offDDT+=joint.nrDDT();
  }
  _body.fillChildren();
  //2D plane
  if(dim==2) {
    ScalarUtil<scalarD>::ScalarQuat quat=ScalarUtil<scalarD>::ScalarQuat::FromTwoVectors(planeNormal.normalized(),Vec3d::UnitZ());
    std::set<sizeType> children=_body.children(1,true);
    for(sizeType i:children) {
      Joint& J=_body.joint(i);
      J._trans=(quat.toRotationMatrix()*J._trans).eval();
    }
    Joint& J0=_body.joint(0);
    J0._trans=(quat.toRotationMatrix().transpose()*J0._trans).eval();
  }
}
void ArticulatedUtils::simplify(std::function<bool(const Joint&)> canSimplify,const Cold& DOF,sizeType nrDebug)
{
  ArticulatedBody bodyOriginal=_body;
  PBDArticulatedGradientInfo<scalarD> info(_body,DOF);
  //remove fixed
  for(sizeType i=1; i<_body.nrJ(); i++) {
    const Joint& joint=_body.joint(i);
    if(canSimplify(joint)) {
      Mat3X4d trans=TRANSI(info._TK_1KM,i);
      sizeType parentId=joint._parent;
      while(parentId>0 && canSimplify(_body.joint(parentId))) {
        Mat3X4d currTrans=trans,pTrans=TRANSI(info._TK_1KM,parentId);
        APPLY_TRANS(trans,pTrans,currTrans)
        parentId=_body.joint(parentId)._parent;
      }
      Joint& parent=_body.joint(parentId);
      ASSERT(parentId==0 || !canSimplify(parent))
      //merge M,MC,MCCT
      parent._M+=joint._M;
      parent._MC+=ROT(trans)*joint._MC+CTR(trans)*joint._M;
      parent._MCCT+=ROT(trans)*joint._MCCT*ROT(trans).transpose();
      parent._MCCT+=ROT(trans)*joint._MC*CTR(trans).transpose();
      parent._MCCT+=CTR(trans)*joint._MC.transpose()*ROT(trans).transpose();
      parent._MCCT+=CTR(trans)*CTR(trans).transpose()*joint._M;
      parent._MCCT=((parent._MCCT+parent._MCCT.transpose())/2).eval();
      //merge name
      parent._name+="+"+joint._name;
      //merge sphere
      Mat3Xd spheres=ROT(trans)*joint._spheres+CTR(trans)*Cold::Ones(joint._spheres.cols()).transpose();
      parent._spheres=concatCol(parent._spheres,spheres);
      parent._radGeomColl=concat(parent._radGeomColl,joint._radGeomColl);
      parent._radSelfColl=concat(parent._radSelfColl,joint._radSelfColl);
      //merge mesh
      ObjMesh m=joint.getMesh(Joint::MESH);
      m.getT()=ROT(trans).cast<scalar>();
      m.getPos()=CTR(trans).cast<scalar>();
      m.applyTrans(Vec3::Zero());
      ObjMesh mParent=parent.getMesh(Joint::MESH);
      mParent.addMesh(m);
      if(mParent.getV().empty())
        parent._mesh=NULL;
      else parent._mesh.reset(new ObjMeshGeomCell(Mat4::Identity(),mParent,0,true,false));
    }
  }
  //re-index
  ArticulatedBody body;
  std::map<sizeType,sizeType> idMap;
  body._geom.reset(new StaticGeom(3));
  body._geomEnv=_body._geomEnv;
  for(sizeType i=0; i<_body.nrJ(); i++) {
    Joint& joint=_body.joint(i);
    if(_body.children(i).empty() && joint.getMesh(Joint::MESH).getV().empty())
      continue;
    else if(i==0 || !canSimplify(joint)) {
      while(joint._parent>0 && canSimplify(_body.joint(joint._parent))) {
        Mat3X4d currTrans=joint._trans,pTrans=TRANSI(info._TK_1KM,joint._parent);
        APPLY_TRANS(joint._trans,pTrans,currTrans)
        joint._parent=_body.joint(joint._parent)._parent;
      }
      ASSERT(joint._parent==0 || !canSimplify(_body.joint(joint._parent)))
      idMap[i]=body._joints.size();
      body._joints.push_back(joint);
    }
  }
  //assemble
  sizeType offDOF=0,offDDT=0;
  for(sizeType i=0; i<(sizeType)body._joints.size(); i++) {
    Joint& joint=body._joints[i];
    joint._parent=i==0?-1:idMap.find(joint._parent)->second;
    if(joint._mimic>=0)
      joint._mimic=idMap.find(joint._mimic)->second;
    joint._offDOF=offDOF;
    joint._offDDT=offDDT;
    offDOF+=joint.nrDOF();
    offDDT+=joint.nrDDT();
    //geom
    if(joint._mesh)
      body._geom->addGeomCell(joint._mesh);
  }
  body._geom->assemble();
  body.fillChildren();
  //replace
  _body=body;
  //debug
  DEFINE_NUMERIC_DELTA_T(scalarD)
  for(sizeType i=0; i<nrDebug; i++) {
    Cold xOriginal=bodyOriginal.randomPose(1),x=Cold::Zero(0);
    for(sizeType j=0; j<bodyOriginal.nrJ(); j++) {
      const Joint& J=bodyOriginal.joint(j);
      if(j==0 || !canSimplify(J))
        x=concat<Cold>(x,xOriginal.segment(J._offDOF,J.nrDOF()));
      else xOriginal.segment(J._offDOF,J.nrDOF())=DOF.segment(J._offDOF,J.nrDOF());
    }
    //EOriginal
    scalarD EOriginal=0;
    PBDArticulatedGradientInfo<scalarD> infoOriginal(bodyOriginal,xOriginal);
    for(sizeType j=0; j<bodyOriginal.nrJ(); j++) {
      const Joint& J=bodyOriginal.joint(j);
      const Mat3d& PPT=J._MCCT;
      const Vec3d& P=J._MC;
      const Mat3X4d A=TRANSI(infoOriginal._TM,j);
      EOriginal+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace();
    }
    //E
    scalarD E=0;
    ASSERT_MSG(x.size()==body.nrDOF(),"Error estimating simplified body DOF!")
    info=PBDArticulatedGradientInfo<scalarD>(body,x);
    for(sizeType j=0; j<body.nrJ(); j++) {
      const Joint& J=body.joint(j);
      const Mat3d& PPT=J._MCCT;
      const Vec3d& P=J._MC;
      const Mat3X4d A=TRANSI(info._TM,j);
      E+=(ROT(A)*PPT*ROT(A).transpose()+2*CTR(A)*P.transpose()*ROT(A).transpose()+CTR(A)*CTR(A).transpose()*J._M).trace();
    }
    DEBUG_GRADIENT("simplifyConsistency",E,E-EOriginal)
  }
}
void ArticulatedUtils::simplify(sizeType nrDebug)
{
  simplify([](const Joint& J) {
    return J._typeJoint==Joint::FIX_JOINT;
  },Cold::Zero(_body.nrDOF()),nrDebug);
}
void ArticulatedUtils::eliminateJoint(const std::vector<std::string>& jointName,const Cold& DOF,sizeType nrDebug)
{
  simplify([&](const Joint& J) {
    return std::find(jointName.begin(),jointName.end(),J._name)!=jointName.end();
  },DOF,nrDebug);
}
void ArticulatedUtils::scaleMass(scalarD coef)
{
  for(sizeType j=0; j<_body.nrJ(); j++) {
    Joint& J=_body.joint(j);
    J._M*=coef;
    J._MC*=coef;
    J._MCCT*=coef;
  }
}
scalarD ArticulatedUtils::totalMass() const
{
  scalarD M=0;
  for(sizeType j=0; j<_body.nrJ(); j++) {
    const Joint& J=_body.joint(j);
    M+=J._M;
  }
  return M;
}
