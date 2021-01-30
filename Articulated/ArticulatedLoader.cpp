#include "PBDArticulatedGradientInfo.h"
#include "ArticulatedLoader.h"
#include "ArticulatedUtils.h"
#include "ConvexHull.h"
#include <Utils/Utils.h>
#include <Utils/RotationUtil.h>
#include <Utils/ArticulatedBodyPragma.h>
#include <assimp/scene.h>
#include <assimp/vector3.h>
#include <assimp/Importer.hpp>
#include <CommonFile/MakeMesh.h>
#include <CommonFile/DisjointSet.h>
#include <CommonFile/geom/BVHBuilder.h>
#include <CommonFile/geom/ObjMeshGeomCell.h>
#include <Eigen/Eigen>

USE_PRJ_NAMESPACE

ObjMesh ArticulatedLoader::createMeshFromScene(const aiScene* scene,const aiNode* node)
{
  ObjMesh mesh;
  for(sizeType i=0; i<node->mNumChildren; i++)
    mesh.addMesh(createMeshFromScene(scene,node->mChildren[i]));
  for(sizeType i=0; i<node->mNumMeshes; i++) {
    ObjMesh comp;
    const aiMesh* m=scene->mMeshes[node->mMeshes[i]];
    for(sizeType i=0; i<m->mNumVertices; i++) {
      const aiVector3D& v=m->mVertices[i];
      comp.getV().push_back(Vec3(v.x,v.y,v.z));
    }
    for(sizeType i=0; i<m->mNumFaces; i++) {
      const aiFace& f=m->mFaces[i];
      for(sizeType j=0; j<f.mNumIndices-2; j++)
        comp.getI().push_back(Vec3i(f.mIndices[0],f.mIndices[j+1],f.mIndices[j+2]));
    }
    mesh.addMesh(comp);
  }
  Mat4 T;
  T(0,0)=node->mTransformation.a1;
  T(0,1)=node->mTransformation.a2;
  T(0,2)=node->mTransformation.a3;
  T(0,3)=node->mTransformation.a4;
  //
  T(1,0)=node->mTransformation.b1;
  T(1,1)=node->mTransformation.b2;
  T(1,2)=node->mTransformation.b3;
  T(1,3)=node->mTransformation.b4;
  //
  T(2,0)=node->mTransformation.c1;
  T(2,1)=node->mTransformation.c2;
  T(2,2)=node->mTransformation.c3;
  T(2,3)=node->mTransformation.c4;
  //
  T(3,0)=node->mTransformation.d1;
  T(3,1)=node->mTransformation.d2;
  T(3,2)=node->mTransformation.d3;
  T(3,3)=node->mTransformation.d4;
  for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++) {
    Vec4 v=T*Vec4(mesh.getV(i)[0],mesh.getV(i)[1],mesh.getV(i)[2],1);
    mesh.getV(i)=v.segment<3>(0)/v[3];
  }
  return mesh;
}
Mat3 ArticulatedLoader::RPY2Mat(const Vec3& rpy)
{
  Mat3 RX=expWGradV<scalar,Vec3>(Vec3::UnitX()*rpy[0]);
  Mat3 RY=expWGradV<scalar,Vec3>(Vec3::UnitY()*rpy[1]);
  Mat3 RZ=expWGradV<scalar,Vec3>(Vec3::UnitZ()*rpy[2]);
  return RZ*RY*RX;
}
void ArticulatedLoader::mergeMesh(ObjMesh& mesh)
{
  mesh.smooth();
  //build BVH
  ObjMesh::EdgeMap eMap;
  mesh.buildEdge(eMap);
  std::vector<Node<std::pair<int,int>>> bvh;
  for(std::map<std::pair<int,int>,ObjMesh::Edge,ObjMesh::EdgeMap::LSS>::const_iterator
      beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++)
  {
    Node<std::pair<int,int>> n;
    n._cell=beg->first;
    n._bb.setUnion(mesh.getV(beg->first.first));
    n._bb.setUnion(mesh.getV(beg->first.second));
    n._bb.enlarged(1E-3f);
    bvh.push_back(n);
  }
  BVHBuilder<Node<std::pair<int,int>>,3> builder;
  builder.buildBVH(bvh);
  for(sizeType i=(sizeType)eMap._ess.size(); i<(sizeType)bvh.size(); i++)
    bvh[i]._cell=std::make_pair(-1,-1);
  //cluster edge
  DisjointSet<sizeType> ds((sizeType)mesh.getV().size());
  for(std::map<std::pair<int,int>,ObjMesh::Edge,ObjMesh::EdgeMap::LSS>::const_iterator
      beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++)
  {
    BBox<scalar> bb;
    bb.setUnion(mesh.getV(beg->first.first));
    bb.setUnion(mesh.getV(beg->first.second));
    ASSERT(beg->second._tris.size()<=2)
    if(beg->second._tris.size()<2) {
      std::stack<sizeType> ss;
      ss.push((sizeType)bvh.size()-1);
      sizeType matchFirst=-1,matchSecond=-1;
      scalar bestScore=ScalarUtil<scalar>::scalar_max(),score1,score2;
      while(!ss.empty()) {
        const Node<std::pair<int,int>>& n=bvh[ss.top()];
        ss.pop();
        if(!n._bb.intersect(bb))
          continue;
        else if(n._cell.first==-1 && n._cell.second==-1) {
          ss.push(n._l);
          ss.push(n._r);
        } else if(beg->first==n._cell)
          continue;
        else {
          score1=(mesh.getV(beg->first.first)-mesh.getV(n._cell.first)).norm()+(mesh.getV(beg->first.second)-mesh.getV(n._cell.second)).norm();
          score2=(mesh.getV(beg->first.first)-mesh.getV(n._cell.second)).norm()+(mesh.getV(beg->first.second)-mesh.getV(n._cell.first)).norm();
          if(score1<bestScore) {
            bestScore=score1;
            matchFirst=n._cell.first;
            matchSecond=n._cell.second;
          }
          if(score2<bestScore) {
            bestScore=score2;
            matchFirst=n._cell.second;
            matchSecond=n._cell.first;
          }
        }
      }
      //use best match
      if(bestScore==0) {
        ds.joinSafe(beg->first.first,matchFirst);
        ds.joinSafe(beg->first.second,matchSecond);
      }
    }
  }
  //realize the merge
  ObjMesh meshM;
  sizeType nrVert=0;
  for(sizeType i=0; i<(sizeType)ds._elts.size(); i++)
    if(ds._elts[i]._p==i) {
      meshM.getV().push_back(mesh.getV(i));
      ds._elts[i]._Int=nrVert++;
    }
  for(sizeType i=0; i<(sizeType)mesh.getI().size(); i++)
    meshM.getI().push_back(Vec3i(ds.weightR(mesh.getI(i)[0]),ds.weightR(mesh.getI(i)[1]),ds.weightR(mesh.getI(i)[2])));
  meshM.smooth();
  mesh=meshM;
  //let's check
  //mesh.writeVTK("mesh.vtk",true);
  mesh.buildEdge(eMap);
}
void ArticulatedLoader::makeConvex(tinyxml2::XMLElement& pt)
{
  for(tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(std::string(v->Name())=="link")
      put<char>(*v,"convex",1);
    else makeConvex(*v);
  }
}
ObjMesh ArticulatedLoader::readMesh(const std::experimental::filesystem::v1::path& path,bool adjust)
{
  Assimp::Importer importer;
  const aiScene *scene=importer.ReadFile(path.string().c_str(),0);//NULL);
  ObjMesh mesh=createMeshFromScene(scene,scene->mRootNode);
  if(adjust) {
    mergeMesh(mesh);
    //make sure positive volume
    mesh.smooth();
    //mesh.makeUniform();
    if(mesh.getVolume()<0)
      mesh.insideOut();
  }
  mesh.smooth();
  return mesh;
}
void ArticulatedLoader::readLinks(Links& links,const tinyxml2::XMLElement& pt,const std::experimental::filesystem::v1::path& path,bool visualMesh)
{
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {

    if(std::string(v->Name())=="gazebo")
    {
      continue;
    }

    
    if(std::string(v->Name())=="link") {
      const tinyxml2::XMLElement& link=*v;
      if(!hasAttribute(link,"<xmlattr>.name"))
      {
        continue;
      }
 
      //name
      Joint joint;
      ObjMesh meshAll;
      joint._name=get<std::string>(link,"<xmlattr>.name");
      //mesh
      for(const tinyxml2::XMLElement* g=link.FirstChildElement(); g; g=g->NextSiblingElement()) {
        if(visualMesh && std::string(g->Name())!="visual")
          continue;
        if(!visualMesh && std::string(g->Name())!="collision")
          continue;
        const tinyxml2::XMLElement& meshPt=*g;
        //T
        Mat4 T=Mat4::Identity();
        ROT(T)=RPY2Mat(parsePtreeDef<Vec3>(meshPt,"origin.<xmlattr>.rpy",Vec3::Zero()));
        CTR(T)=parsePtreeDef<Vec3>(meshPt,"origin.<xmlattr>.xyz",Vec3::Zero());
        //geometry
        ObjMesh mesh;
        if(hasAttribute(meshPt,"geometry.mesh")) {
          if(!hasAttribute(meshPt,"geometry.mesh.<xmlattr>.filename"))
            continue;
          mesh=readMesh(path/get<std::string>(meshPt,"geometry.mesh.<xmlattr>.filename"),!visualMesh);
          if(get<char>(link,"convex",(char)0))
            mesh=::makeConvex(mesh);
          Vec3 scale=parsePtreeDef<Vec3>(meshPt,"geometry.mesh.<xmlattr>.scale",Vec3::Ones());
          mesh.getT()=scale.asDiagonal();
          mesh.applyTrans(Vec3::Zero());
        } else if(hasAttribute(meshPt,"geometry.box")) {
          std::cout << "Box" << std::endl;
          Vec3 size=parsePtree<Vec3>(meshPt,"geometry.box.<xmlattr>.size");
          std::cout << "size = " << size[0] << " " << size[1] << " " << size[2] << std::endl;
          MakeMesh::makeBox3D(mesh,size/2);
          joint._mesh=std::shared_ptr<ObjMeshGeomCell>(new ObjMeshGeomCell(T,mesh,0,true,false));
        } else if(hasAttribute(meshPt,"geometry.sphere")) {
          scalar r=get<scalar>(meshPt,"geometry.sphere.<xmlattr>.radius");
          MakeMesh::makeSphere3D(mesh,r,16);
        } else if(hasAttribute(meshPt,"geometry.cylinder")) {
          scalar r=get<scalar>(meshPt,"geometry.cylinder.<xmlattr>.radius");
          scalar l=get<scalar>(meshPt,"geometry.cylinder.<xmlattr>.length");
          MakeMesh::makeCylinder3D(mesh,r,l/2,16);
          mesh.getT()=expWGradV<scalar,Vec3>(Vec3::UnitX()*M_PI/2);
          mesh.applyTrans(Vec3::Zero());
        }
        if(!mesh.getV().empty()) {
          mesh.getT()=ROT(T);
          mesh.getPos()=CTR(T);
          mesh.applyTrans(Vec3::Zero());
          meshAll.addMesh(mesh);
        }
      }
      //assemble
      if(!meshAll.getV().empty())
        joint._mesh=std::shared_ptr<ObjMeshGeomCell>(new ObjMeshGeomCell(Mat4::Identity(),meshAll,0,true,false));
      if(joint._mesh) {
        if(hasAttribute(link,"inertial.mass.<xmlattr>.value") && !get<char>(link,"convex",(char)0)) { //if convex is required, recompute mass
          INFO("Using mass in URDF file!")
          Mat3d inertia;
          inertia(0,0)=get<scalar>(link,"inertial.inertia.<xmlattr>.ixx");
          inertia(1,1)=get<scalar>(link,"inertial.inertia.<xmlattr>.iyy");
          inertia(2,2)=get<scalar>(link,"inertial.inertia.<xmlattr>.izz");
          inertia(0,1)=inertia(1,0)=get<scalar>(link,"inertial.inertia.<xmlattr>.ixy");
          inertia(0,2)=inertia(2,0)=get<scalar>(link,"inertial.inertia.<xmlattr>.ixz");
          inertia(1,2)=inertia(2,1)=get<scalar>(link,"inertial.inertia.<xmlattr>.iyz");
          Mat3d MCCT=Mat3d::Identity()*inertia.trace()/2-inertia;
          Eigen::SelfAdjointEigenSolver<Mat3d> eig(MCCT);
          if(eig.eigenvalues().minCoeff()<0) {
            std::cout << "Found negative MCCT tensor (link=" << joint._name << "): " << eig.eigenvalues().transpose() << std::endl;
          }
          //convert
          Mat3d R=RPY2Mat(parsePtreeDef<Vec3>(link,"inertial.origin.<xmlattr>.rpy",Vec3::Zero())).cast<scalarD>();
          Vec3d t=parsePtreeDef<Vec3>(link,"inertial.origin.<xmlattr>.xyz",Vec3::Zero()).cast<scalarD>();
          //build 6x6 matrix in body frame
          scalar mass=get<scalar>(link,"inertial.mass.<xmlattr>.value");
          joint._M=mass;
          joint._MC=mass*t;
          joint._MCCT=R*MCCT*R.transpose()+mass*t*t.transpose();
        } else {
          WARNING("Using recomputed mass!")
          std::cout << "Joint name is: " << joint._name << std::endl;
          // if (joint._name == "thbase" || joint._name == "thhub")
          // {
          //   INFO("Using mass in URDF file!")
          // Mat3d inertia;
          // inertia(0,0)=get<scalar>(link,"inertial.inertia.<xmlattr>.ixx");
          // inertia(1,1)=get<scalar>(link,"inertial.inertia.<xmlattr>.iyy");
          // inertia(2,2)=get<scalar>(link,"inertial.inertia.<xmlattr>.izz");
          // inertia(0,1)=inertia(1,0)=get<scalar>(link,"inertial.inertia.<xmlattr>.ixy");
          // inertia(0,2)=inertia(2,0)=get<scalar>(link,"inertial.inertia.<xmlattr>.ixz");
          // inertia(1,2)=inertia(2,1)=get<scalar>(link,"inertial.inertia.<xmlattr>.iyz");
          // Mat3d MCCT=Mat3d::Identity()*inertia.trace()/2-inertia;
          // Eigen::SelfAdjointEigenSolver<Mat3d> eig(MCCT);
          // if(eig.eigenvalues().minCoeff()<0) {
          //   std::cout << "Found negative MCCT tensor (link=" << joint._name << "): " << eig.eigenvalues().transpose() << std::endl;
          // }
          // //convert
          // Mat3d R=RPY2Mat(parsePtreeDef<Vec3>(link,"inertial.origin.<xmlattr>.rpy",Vec3::Zero())).cast<scalarD>();
          // Vec3d t=parsePtreeDef<Vec3>(link,"inertial.origin.<xmlattr>.xyz",Vec3::Zero()).cast<scalarD>();
          // //build 6x6 matrix in body frame
          // scalar mass=get<scalar>(link,"inertial.mass.<xmlattr>.value");
          // joint._M=mass;
          // joint._MC=mass*t;
          // joint._MCCT=R*MCCT*R.transpose()+mass*t*t.transpose();
          // }
          // else
          // {
            joint.assemble(1);
          // }
        }
      } else {
        joint._M=0;
        joint._MC.setZero();
        joint._MCCT.setZero();
      }
      links[joint._name]=joint;
    } else 
    {
      readLinks(links,*v,path,visualMesh);
    }
  }
}
bool ArticulatedLoader::readRelations(Relations& relations,const tinyxml2::XMLElement& pt)
{
  for(const tinyxml2::XMLElement* v=pt.FirstChildElement(); v; v=v->NextSiblingElement()) {
    if(std::string(v->Name())=="gazebo")
      continue;
    if(std::string(v->Name())=="joint") {
      JointInfo info;
      const tinyxml2::XMLElement& joint=*v;
      //name
      if(!hasAttribute(joint,"<xmlattr>.name"))
        continue;
      info._name=get<std::string>(joint,"<xmlattr>.name");
      //type
      if(!hasAttribute(joint,"<xmlattr>.type"))
        continue;
      if(get<std::string>(joint,"<xmlattr>.type")!="revolute" && get<std::string>(joint,"<xmlattr>.type")!="fixed" &&
          get<std::string>(joint,"<xmlattr>.type")!="continuous" && get<std::string>(joint,"<xmlattr>.type")!="prismatic")
        return false;
      info._type=get<std::string>(joint,"<xmlattr>.type");
      //parent
      if(!hasAttribute(joint,"parent.<xmlattr>.link"))
        continue;
      //child
      if(!hasAttribute(joint,"child.<xmlattr>.link"))
        continue;
      info._child=get<std::string>(joint,"child.<xmlattr>.link");
      //param
      info._T.setIdentity();
      ROT(info._T)=RPY2Mat(parsePtreeDef<Vec3>(joint,"origin.<xmlattr>.rpy",Vec3::Zero()));
      CTR(info._T)=parsePtreeDef<Vec3>(joint,"origin.<xmlattr>.xyz",Vec3::Zero());
      info._axis=parsePtreeDef<Vec3>(joint,"axis.<xmlattr>.xyz",Vec3::Zero());
      info._ctrl=get<scalar>(joint,"limit.<xmlattr>.effort",0);
      info._damp=get<scalar>(joint,"dynamics.<xmlattr>.damping",0);
      //lower
      if(hasAttribute(joint,"limit.<xmlattr>.lower"))
        info._lmt[0]=get<scalar>(joint,"limit.<xmlattr>.lower");
      else info._lmt[0]=-ScalarUtil<scalar>::scalar_inf();
      //upper
      if(hasAttribute(joint,"limit.<xmlattr>.upper"))
        info._lmt[1]=get<scalar>(joint,"limit.<xmlattr>.upper");
      else info._lmt[1]=ScalarUtil<scalar>::scalar_inf();
      //mimic
      if(get<std::string>(joint,"<xmlattr>.type")=="revolute") {
        if(hasAttribute(joint,"mimic.<xmlattr>.joint")) {
          info._mimic=get<std::string>(joint,"mimic.<xmlattr>.joint");
          info._multiplier=get<scalar>(joint,"mimic.<xmlattr>.multiplier",0);
          info._offset=get<scalar>(joint,"mimic.<xmlattr>.offset",0);
        }
      }
      relations[get<std::string>(joint,"parent.<xmlattr>.link")].push_back(info);
    } else if(!readRelations(relations,*v))
      return false;
  }
  return true;
}
std::string ArticulatedLoader::findRelationByChild(const std::string& name,const Relations& relations,JointInfo& info)
{
  for(const std::pair<std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      if(c._child==name) {
        info=c;
        return v.first;
      }
    }
  }
  ASSERT_MSGV(false,"Cannot find parent for link: %s!",name.c_str())
  return name;
}
std::string ArticulatedLoader::findRelationByName(const std::string& name,const Relations& relations,JointInfo& info)
{
  for(const std::pair<std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      if(c._name==name) {
        info=c;
        return v.first;
      }
    }
  }
  ASSERT_MSGV(false,"Cannot find parent with name: %s!",name.c_str())
  return name;
}
bool ArticulatedLoader::buildBody(ArticulatedBody& body,const std::string& name,const Links& links,const Relations& relations)
{
  sizeType JID=-1;
  Mat3X4d invA2Z=Mat3X4d::Identity();
  if(body._joints.size()==0) {
    //this is root
    if(links.find(name)==links.end()) {
      WARNINGV("Cannot find link: %s!",name.c_str())
      return false;
    }
    body._joints.push_back(links.find(name)->second);
    Joint& joint=body._joints[JID=0];
    joint._parent=-1;
    joint._depth=0;
    joint._typeJoint=Joint::FIX_JOINT;
    joint._offDOF=0;
    joint._offDDT=0;
    joint._limits.resize(3,0);
    joint._control.resize(0);
    joint._damping.resize(0);
    joint._trans.setIdentity();
  } else {
    //this is non-root, then get info
    JointInfo info;
    sizeType nrDOF=body.nrDOF();
    sizeType nrDDT=body.nrDDT();
    const std::string& parentName=findRelationByChild(name,relations,info);
    //build
    if(links.find(name)==links.end()) {
      WARNINGV("Cannot find link: %s!",name.c_str())
      return false;
    }
    body._joints.push_back(links.find(name)->second);
    Joint& joint=body._joints[JID=body.nrJ()-1];
    joint._parent=body.jointId(parentName);
    joint._depth=body.joint(joint._parent)._depth+1;
    joint._mult=joint._offset=0;
    if(info._type=="revolute" || info._type=="continuous") {
      joint._typeJoint=Joint::HINGE_JOINT;
      if(!info._mimic.empty()) {
        joint._mult=info._multiplier;
        joint._offset=info._offset;
      }
    } else if(info._type=="fixed")
      joint._typeJoint=Joint::FIX_JOINT;
    else if(info._type=="prismatic") {
      joint._typeJoint=Joint::TRANS_1D;
    } else {
      ASSERT_MSGV(false,"Unknown type: %s!",info._type.c_str())
    }
    joint._offDOF=nrDOF;
    joint._offDDT=nrDDT;
    joint._limits=Vec3d(info._lmt[0],info._lmt[1],0);
    if(!std::isfinite(info._lmt[0]) || !std::isfinite(info._lmt[1]))
      joint._limits(2,0)=ScalarUtil<scalarD>::scalar_inf();
    joint._control.resize(1);
    joint._control[0]=info._ctrl;
    joint._damping.resize(1);
    joint._damping[0]=info._damp;
    //Axis2Z
    Mat3X4d A2Z=Mat3X4d::Identity(),trans0=info._T.cast<scalarD>().block<3,4>(0,0);
    if(joint._typeJoint==Joint::HINGE_JOINT)
      ROT(A2Z)=ScalarUtil<scalarD>::ScalarQuat::FromTwoVectors(info._axis.cast<scalarD>(),Vec3d::UnitZ()).toRotationMatrix();
    if(joint._typeJoint==Joint::TRANS_1D)
      ROT(A2Z)=ScalarUtil<scalarD>::ScalarQuat::FromTwoVectors(info._axis.cast<scalarD>(),Vec3d::UnitX()).toRotationMatrix();
    APPLY_TRANS(joint._trans,trans0,A2Z)
    //joint
    if(joint._mesh) {
      Mat4d T=joint._mesh->getT().cast<scalarD>();
      ASSERT(T==Mat4d::Identity())
      ROT(T)=ROT(A2Z).transpose()*ROT(T);
      CTR(T)=ROT(A2Z).transpose()*CTR(T);
      joint._mesh->setT(T.cast<scalar>());
//#define DEBUG_JOINT_TRANSFORM
#ifdef DEBUG_JOINT_TRANSFORM
      Joint tmpJ=joint;
      tmpJ.assemble(1);
      scalarD coef=joint._M/tmpJ._M;
      tmpJ._M*=coef;
      tmpJ._MC*=coef;
      tmpJ._MCCT*=coef;
#endif
      Vec3d MC=joint._MC;
      //update MC
      joint._MC=(ROT(T)*joint._MC+CTR(T)).eval();
      //update MCCT
      joint._MCCT=(ROT(T)*joint._MCCT*ROT(T).transpose()).eval();
      joint._MCCT+=2*CTR(T)*MC.transpose()*ROT(T).transpose();
      joint._MCCT+=CTR(T)*CTR(T).transpose()*joint._M;
#ifdef DEBUG_JOINT_TRANSFORM
      INFOV("MErr: %f",joint._M-tmpJ._M)
      INFOV("MCErr: %f",(joint._MC-tmpJ._MC).norm())
      INFOV("MCCTErr: %f",(joint._MCCT-tmpJ._MCCT).norm())
#endif
    }
    INV(invA2Z,A2Z)
  }
  //recurse
  if(relations.find(name)!=relations.end()) {
    for(const JointInfo& c:relations.find(name)->second) {
      if(!buildBody(body,c._child,links,relations))
        return false;
    }
  }
  //children
  std::set<sizeType> ids=body.children(JID,true);
  for(sizeType c:ids) {
    Joint& joint=body._joints[c];
    Mat3X4d trans0=joint._trans;
    APPLY_TRANS(joint._trans,invA2Z,trans0)
  }
  return true;
}
ArticulatedBody ArticulatedLoader::readURDF(const std::string& file,bool convex,bool visualMesh)
{
  ASSERT_MSGV(exists(file),"Cannot find %s",file.c_str())
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  if(convex)
    makeConvex(*(pt.RootElement()));
  //read links

  Links links;
  readLinks(links,*(pt.RootElement()),std::experimental::filesystem::v1::path(file).parent_path(),visualMesh);
  //read joints
  Relations relations;
  if(!readRelations(relations,*(pt.RootElement())))
    return ArticulatedBody();
  //merge
  std::unordered_set<std::string> rootFilter;
  for(const std::pair<std::string,std::vector<JointInfo>>& v:relations) {
    rootFilter.insert(v.first);
  }
  for(const std::pair<std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      rootFilter.erase(c._child);
    }
  }
  ArticulatedBody body;
  ASSERT(rootFilter.size()==1)
  if(!buildBody(body,*(rootFilter.begin()),links,relations))
    return ArticulatedBody();
  body.fillChildren();
  //geom
  body._geom.reset(new StaticGeom(3));
  for(sizeType i=0; i<body.nrJ(); i++)
    if(body.joint(i)._mesh)
      body._geom->addGeomCell(body.joint(i)._mesh);
  body._geom->assemble();
  //mimic

  for(const std::pair<std::string,std::vector<JointInfo>>& v:relations) {
    for(const JointInfo& c:v.second) {
      if(!c._mimic.empty()) {
        JointInfo info;
        Joint& joint=body.joint(body.jointId(c._child));
        findRelationByName(c._mimic,relations,info);
        joint._mimic=body.jointId(info._child);
      }
    }
  }
  return body;
}
void ArticulatedLoader::compareVisualMeshURDF(const std::string& file)
{
  tinyxml2::XMLDocument pt;
  pt.LoadFile(file.c_str());
  //read links
  Links links,linksRef;
  readLinks(linksRef,*(pt.RootElement()),std::experimental::filesystem::v1::path(file).parent_path(),false);
  readLinks(links,*(pt.RootElement()),std::experimental::filesystem::v1::path(file).parent_path(),true);
  //visualMesh
  recreate("visualMesh");
  for(const std::pair<std::string,Joint>& v:links) {
    linksRef[v.first].getMesh(Joint::MESH).writeVTK("visualMesh/"+v.first+".vtk",true);
    links[v.first].getMesh(Joint::MESH).writeVTK("visualMesh/"+v.first+"Visual.vtk",true);
  }
}
//visualize
void ArticulatedLoader::visualizeJointLimitVTK(const ArticulatedBody& body,const std::string& path,Joint::GEOM_TYPE type)
{
  recreate(path);
  sizeType nrJ=body.nrJ();
  Cold x0=Cold::Zero(body.nrDOF()),x,random;
  body.mimic(mapV(x0));
  //ref
  PBDArticulatedGradientInfo<scalarD> info(body,x0);
  body.writeVTK(info._TM,path+"/ref.vtk",type);
  for(sizeType i=0; i<nrJ; i++) {
    bool hasLimit=false;
    //limit
    random=x0;
    const Joint& J=body.joint(i);
    if(J._mimic>=0)
      continue;
    for(sizeType j=0; j<J.nrDOF(); j++)
      if(std::isfinite(J._limits(2,j))) {
        hasLimit=true;
        //lower
        x=x0;
        x[J._offDOF+j]=J._limits(0,j);
        body.mimic(mapV(x));
        info.reset(body,x);
        body.writeVTK(info._TM,path+"/lower"+std::to_string(J._offDOF+j)+".vtk",type);
        //upper
        x=x0;
        x[J._offDOF+j]=J._limits(1,j);
        body.mimic(mapV(x));
        info.reset(body,x);
        body.writeVTK(info._TM,path+"/upper"+std::to_string(J._offDOF+j)+".vtk",type);
        //random
        random[J._offDOF+j]=RandEngine::randR(J._limits(0,j),J._limits(1,j));
      }
    //random
    if(hasLimit) {
      info.reset(body,random);
      body.writeVTK(info._TM,path+"/random"+std::to_string(i)+".vtk",type);
    }
  }
}
//built-in bodies
sizeType ArticulatedLoader::inferDim(sizeType rootType)
{
  sizeType dim=2;
  Joint::loopAllJointTypes([&](Joint::JOINT_TYPE t) {
    if(t&rootType) {
      if(t==Joint::TRANS_3D || t==Joint::ROT_3D_EXP || t==Joint::ROT_3D_XYZ || t==Joint::BALL_JOINT)
        dim=3;
    }
  });
  return dim;
}
tinyxml2::XMLElement* ArticulatedLoader::addRootJoint(tinyxml2::XMLElement& pt,sizeType rootType)
{
  sizeType nrRotJoint=0;
  sizeType nrTransJoint=0;
  Joint::JOINT_TYPE rotJoint=Joint::NR_JOINT_TYPE;
  Joint::JOINT_TYPE transJoint=Joint::NR_JOINT_TYPE;
  Joint::loopAllJointTypes([&](Joint::JOINT_TYPE t) {
    if(t&rootType) {
      Joint J;
      J._typeJoint=t;
      if(J.RBegEnd()[1]-J.RBegEnd()[0]) {
        rotJoint=t;
        nrRotJoint++;
      } else if(J.CBegEnd()[1]-J.CBegEnd()[0]) {
        transJoint=t;
        nrTransJoint++;
      }
    }
  });
  ASSERT_MSG(nrRotJoint<=1 && nrTransJoint<=1,"You can only have at most 1 rotJoint and 1 transJoint!")
  tinyxml2::XMLElement* ret=NULL;
  if(nrTransJoint>0) {
    ret=addChild(ret?*ret:pt,"joint");
    put<std::string>(*ret,"type",Joint::typeToString(transJoint));
  }
  if(nrRotJoint>0) {
    ret=addChild(ret?*ret:pt,"joint");
    put<std::string>(*ret,"type",Joint::typeToString(rotJoint));
  }
  if(!ret) {
    ret=addChild(ret?*ret:pt,"joint");
    put<std::string>(*ret,"type",Joint::typeToString(Joint::FIX_JOINT));
  }
  return ret;
}
tinyxml2::XMLElement* ArticulatedLoader::createBird(tinyxml2::XMLElement& pt,sizeType rootType,scalar bodySz,scalar bodyLen,scalar neckSz,scalar neckLen,scalar footLen1,scalar footLen2,scalar footLen3,scalar footRad1,scalar footRad2,scalar footRad3,bool fixFoot,bool head)
{
  Mat3 R0=expWGradV<scalar,Vec3>(Vec3(0,0,D2R(-85.0f)));
  Mat3 R1=expWGradV<scalar,Vec3>(Vec3(0,0,D2R(-45.0f)));
  Mat3 R2=expWGradV<scalar,Vec3>(Vec3(0,0,D2R(30.0f)));
  Mat3 R3=expWGradV<scalar,Vec3>(Vec3(0,0,D2R(-10.0f)));
  tinyxml2::XMLElement* body=addRootJoint(pt,rootType);
  put<scalarF>(*body,"geom1D.lenY",bodyLen);
  put<scalarF>(*body,"geom1D.rad",bodySz);
  putPtree<Vec3>(*body,"meshDirY",R0*Vec3::UnitY());
  for(sizeType i=0; i<2; i++) {
    //leg1
    scalarD sgn=i==0?1:-1;
    tinyxml2::XMLElement* leg1=addChild(*body,"joint");
    put<std::string>(*leg1,"type",Joint::typeToString(Joint::BALL_JOINT));
    putPtree(*leg1,"limit.lower",Vec2d(D2R(-60.0f),D2R(-90.0f)));
    putPtree(*leg1,"limit.upper",Vec2d(D2R(60.0f),D2R(30.0f)));
    put<scalarF>(*leg1,"geom1D.lenY",footLen1);
    put<scalarF>(*leg1,"geom1D.rad",footRad1);
    putPtree<Vec3>(*leg1,"trans",Vec3::UnitZ()*bodySz*2/3*sgn-Vec3::UnitY()*bodySz/2+R0*Vec3::UnitY()*bodySz);
    putPtree<Vec3>(*leg1,"jointDirY",Vec3(0,0,1));
    putPtree<Vec3>(*leg1,"jointDirZ",Vec3(1,0,0)*sgn);
    putPtree<Vec3>(*leg1,"meshDirY",R1*Vec3(0,-1,0));

    //leg2
    tinyxml2::XMLElement* leg2=addChild(*leg1,"joint");
    put<std::string>(*leg2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree(*leg2,"limit.lower",Cold::Constant(1,D2R(-90.0f)));
    putPtree(*leg2,"limit.upper",Cold::Constant(1,D2R(90.0f)));
    //putPtree(*leg2,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
    put<scalarF>(*leg2,"geom1D.lenY",footLen2);
    put<scalarF>(*leg2,"geom1D.rad",footRad2);
    put<scalarF>(*leg2,"transY",1);
    putPtree<Vec3>(*leg2,"jointDirZ",Vec3(0,0,1));
    putPtree<Vec3>(*leg2,"meshDirY",R2*Vec3(0,-1,0));

    //leg3
    tinyxml2::XMLElement* leg3=addChild(*leg2,"joint");
    if(fixFoot)
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::FIX_JOINT));
    else {
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::HINGE_JOINT));
      putPtree(*leg3,"limit.lower",Cold::Constant(1,D2R(-45.0f)));
      putPtree(*leg3,"limit.upper",Cold::Constant(1,D2R(45.0f)));
    }
    //putPtree(*leg3,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
    put<scalarF>(*leg3,"geom3D.lenX",footRad3);
    put<scalarF>(*leg3,"geom3D.lenY",footLen3);
    put<scalarF>(*leg3,"geom3D.lenZ",footRad3/2);
    put<scalarF>(*leg3,"geom3D.rad",0);
    put<scalarF>(*leg3,"transY",1);
    putPtree<Vec3>(*leg3,"trans",-Vec3::UnitY()*footRad2);
    putPtree<Vec3>(*leg3,"jointDirZ",Vec3(0,0,1));
    putPtree<Vec3>(*leg3,"meshDirX",Vec3(0,0,1));
    putPtree<Vec3>(*leg3,"meshDirY",Vec3(1,0,0));
    putPtree<Vec3>(*leg3,"meshTrans",-Vec3::UnitY()*footLen3/2);
  }
  if(head) {
    //neck
    tinyxml2::XMLElement* neck=addChild(*body,"joint");
    put<std::string>(*neck,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree(*neck,"limit.lower",Cold::Constant(1,D2R(-90.0f)));
    putPtree(*neck,"limit.upper",Cold::Constant(1,D2R(90.0f)));
    put<scalarF>(*neck,"geom1D.lenY",neckLen);
    put<scalarF>(*neck,"geom1D.rad",neckSz);
    putPtree<Vec3>(*neck,"trans",Vec3::UnitY()*bodySz/2+R0*Vec3::UnitY()*(bodyLen+bodySz/2));
    putPtree<Vec3>(*neck,"jointDirZ",Vec3(0,0,1));
    putPtree<Vec3>(*neck,"meshDirY",Vec3(0,1,0));

    //head
    tinyxml2::XMLElement* head=addChild(*neck,"joint");
    put<std::string>(*head,"type",Joint::typeToString(Joint::FIX_JOINT));
    put<scalarF>(*head,"geom1D.lenY",neckSz*2);
    put<scalarF>(*head,"geom1D.rad",neckSz*1.5f);
    put<scalarF>(*head,"transY",1);
    putPtree<Vec3>(*head,"meshDirY",R3*Vec3(1,0,0));
  }
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createBipedal(tinyxml2::XMLElement& pt,sizeType rootType,scalar bodySz,scalar bodyLen,scalar footLen1,scalar footLen2,scalar footLen3,scalar footRad1,scalar footRad2,scalar footRad3,bool fixFoot)
{
  tinyxml2::XMLElement* body=addRootJoint(pt,rootType);
  put<scalarF>(*body,"geom1D.lenY",bodyLen);
  put<scalarF>(*body,"geom1D.rad",bodySz);
  putPtree<Vec3>(*body,"meshDirY",Vec3::UnitY());
  for(sizeType i=0; i<2; i++) {
    //leg1
    scalarD sgn=i==0?1:-1;
    tinyxml2::XMLElement* leg1=addChild(*body,"joint");
    put<std::string>(*leg1,"type",Joint::typeToString(Joint::BALL_JOINT));
    putPtree(*leg1,"limit.lower",Vec2d(D2R(0.0f),D2R(-60.0f)));
    putPtree(*leg1,"limit.upper",Vec2d(D2R(150.0f),D2R(10.0f)));
    put<scalarF>(*leg1,"geom1D.lenY",footLen1);
    put<scalarF>(*leg1,"geom1D.rad",footRad1);
    putPtree<Vec3>(*leg1,"trans",Vec3::UnitZ()*bodySz*2/3*sgn-Vec3::UnitY()*bodySz/2);
    putPtree<Vec3>(*leg1,"jointDirY",Vec3(0,0,1));
    putPtree<Vec3>(*leg1,"jointDirZ",Vec3(1,0,0)*sgn);
    putPtree<Vec3>(*leg1,"meshDirY",Vec3(0,-1,0));

    //leg2
    tinyxml2::XMLElement* leg2=addChild(*leg1,"joint");
    put<std::string>(*leg2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree(*leg2,"limit.lower",Cold::Constant(1,D2R(-150.0f)));
    putPtree(*leg2,"limit.upper",Cold::Constant(1,D2R(0.0f)));
    //putPtree(*leg2,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
    put<scalarF>(*leg2,"geom1D.lenY",footLen2);
    put<scalarF>(*leg2,"geom1D.rad",footRad2);
    put<scalarF>(*leg2,"transY",1);
    putPtree<Vec3>(*leg2,"jointDirZ",Vec3(0,0,1));
    putPtree<Vec3>(*leg2,"meshDirY",Vec3(0,-1,0));

    //leg3
    tinyxml2::XMLElement* leg3=addChild(*leg2,"joint");
    if(fixFoot)
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::FIX_JOINT));
    else {
      put<std::string>(*leg3,"type",Joint::typeToString(Joint::HINGE_JOINT));
      putPtree(*leg3,"limit.lower",Cold::Constant(1,D2R(-30.0f)));
      putPtree(*leg3,"limit.upper",Cold::Constant(1,D2R(30.0f)));
    }
    //putPtree(*leg3,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
    put<scalarF>(*leg3,"geom3D.lenX",footRad3);
    put<scalarF>(*leg3,"geom3D.lenY",footLen3);
    put<scalarF>(*leg3,"geom3D.lenZ",footRad3/2);
    put<scalarF>(*leg3,"geom3D.rad",0);
    put<scalarF>(*leg3,"transY",1);
    putPtree<Vec3>(*leg3,"trans",-Vec3::UnitY()*footRad2);
    putPtree<Vec3>(*leg3,"jointDirZ",Vec3(0,0,1));
    putPtree<Vec3>(*leg3,"meshDirX",Vec3(0,0,1));
    putPtree<Vec3>(*leg3,"meshDirY",Vec3(1,0,0));
    putPtree<Vec3>(*leg3,"meshTrans",-Vec3::UnitY()*footRad2);
  }
  return body;
}
void ArticulatedLoader::addBipedalHand(tinyxml2::XMLElement& pt,scalar bodySz,scalar bodyLen,scalar handLen1,scalar handLen2,scalar handRad1,scalar handRad2)
{
  tinyxml2::XMLElement* body=&pt;
  for(sizeType i=0; i<2; i++) {
    //leg1
    scalarD sgn=i==0?1:-1;
    tinyxml2::XMLElement* hand1=addChild(*body,"joint");
    put<std::string>(*hand1,"type",Joint::typeToString(Joint::BALL_JOINT));
    putPtree(*hand1,"limit.lower",Vec2d(D2R(-60.0f),D2R(-90.0f)));
    putPtree(*hand1,"limit.upper",Vec2d(D2R(90.0f),D2R(0.0f)));
    put<scalarF>(*hand1,"geom1D.lenY",handLen1);
    put<scalarF>(*hand1,"geom1D.rad",handRad1);
    putPtree<Vec3>(*hand1,"trans",Vec3::UnitZ()*(bodySz+handRad1)*sgn+Vec3::UnitY()*bodyLen);
    putPtree<Vec3>(*hand1,"jointDirY",Vec3(0,0,1));
    putPtree<Vec3>(*hand1,"jointDirZ",Vec3(1,0,0)*sgn);
    putPtree<Vec3>(*hand1,"meshDirY",-Vec3::UnitY());

    //leg2
    tinyxml2::XMLElement* hand2=addChild(*hand1,"joint");
    put<std::string>(*hand2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree(*hand2,"limit.lower",Cold::Constant(1,D2R(0.0f)));
    putPtree(*hand2,"limit.upper",Cold::Constant(1,D2R(150.0f)));
    //putPtree(*hand2,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
    put<scalarF>(*hand2,"geom1D.lenY",handLen2);
    put<scalarF>(*hand2,"geom1D.rad",handRad2);
    put<scalarF>(*hand2,"transY",1);
    putPtree<Vec3>(*hand2,"jointDirZ",Vec3(0,0,1));
    putPtree<Vec3>(*hand2,"meshDirY",Vec3(0,-1,0));
  }
}
tinyxml2::XMLElement* ArticulatedLoader::createChain(tinyxml2::XMLElement& pt,sizeType rootType,sizeType nr,scalar l,scalar rad,scalar rot,scalar rot0,scalar t,scalar t0,sizeType geomDim,scalar yOff,scalar ratioX,scalar ratioZ)
{
  sizeType nrDOF;
  sizeType dim=inferDim(rootType);
  tinyxml2::XMLElement* child=NULL;
  for(sizeType i=0; i<nr; i++) {
    if(i == 0)
      child=addRootJoint(pt,rootType);
    else {
      nrDOF=dim == 2 ? 1 : 2;
      put<std::string>(*child,"type",dim == 2 ? "HINGE_JOINT" : "BALL_JOINT");
      putPtree(*child,"limit.lower",Cold::Constant(nrDOF,-rot));
      putPtree(*child,"limit.upper",Cold::Constant(nrDOF,rot));
      //putPtree(*child,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
      putPtree(*child,"x0",Cold::Constant(nrDOF,rot0));
    }
    //geometry
    put<scalarF>(*child,"geom"+std::to_string(geomDim)+"D.lenX",l*ratioX);
    put<scalarF>(*child,"geom"+std::to_string(geomDim)+"D.lenY",l);
    put<scalarF>(*child,"geom"+std::to_string(geomDim)+"D.lenZ",l*ratioZ);
    put<scalarF>(*child,"geom"+std::to_string(geomDim)+"D.rad",rad);
    //rotation
    put<scalarF>(*child,"transY",1+yOff);
    putPtree(*child,"meshDirX",Vec3d(1,0,0));
    putPtree(*child,"meshDirY",Vec3d(0,-1,0));
    putPtree(*child,"jointDirY",Vec3d(0,0,1));
    putPtree(*child,"jointDirZ",Vec3d(1,0,0));
    //translation
    if(t > 0 && i<nr-1) {
      child=addChild(*child,"joint");
      put<std::string>(*child,"type",Joint::typeToString(Joint::TRANS_1D));
      put<scalarF>(*child,"limit.lower",-t);
      put<scalarF>(*child,"limit.upper",t);
      //put<scalarF>(*child,"limit.coef",1000);  //use global joint limit
      put<scalarF>(*child,"transY",1);
      put<scalarF>(*child,"x0",t0);
      putPtree(*child,"jointDirX",Vec3d(0,-1,0));
    }
    if(i<nr-1)
      child=addChild(*child,"joint");
  }
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createSpider(tinyxml2::XMLElement& pt,sizeType rootType,scalar bodySz,scalar footLen,scalar footRad,scalar angLegLift,scalar rangeLeg1,scalar rangeLeg2,bool ball)
{
  sizeType nrDOF;
  tinyxml2::XMLElement* body=addRootJoint(pt,rootType);
  put<scalarF>(*body,"geom0D.rad",bodySz);
  for(sizeType i=0; i<4; i++) {
    Mat3 R0=expWGradV<scalar,Vec3>(Vec3(0,M_PI/4+(scalar)i*M_PI/2,0));
    Mat3 R1=expWGradV<scalar,Vec3>(Vec3(0,0,angLegLift));
    Mat3 R2=expWGradV<scalar,Vec3>(Vec3(0,0,angLegLift-M_PI/2));

    //leg1
    nrDOF=ball ? 2: 1;
    tinyxml2::XMLElement* leg1=addChild(*body,"joint");
    put<std::string>(*leg1,"type",ball ? "BALL_JOINT" : "HINGE_JOINT");
    putPtree(*leg1,"limit.lower",Cold::Constant(nrDOF,-rangeLeg1));
    putPtree(*leg1,"limit.upper",Cold::Constant(nrDOF,rangeLeg1));
    //putPtree(*leg1,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
    put<scalarF>(*leg1,"geom1D.lenY",footLen);
    put<scalarF>(*leg1,"geom1D.rad",footRad);
    putPtree<Vec3>(*leg1,"trans",R0*Vec3::UnitX()*bodySz);
    putPtree<Vec3>(*leg1,"jointDirY",R0*Vec3(0,0,1));
    putPtree<Vec3>(*leg1,"jointDirZ",Vec3(0,1,0));
    putPtree<Vec3>(*leg1,"meshDirY",R0*R1*Vec3::UnitX());

    //leg2
    nrDOF=1;
    tinyxml2::XMLElement* leg2=addChild(*leg1,"joint");
    put<std::string>(*leg2,"type",Joint::typeToString(Joint::HINGE_JOINT));
    putPtree(*leg2,"limit.lower",Cold::Constant(nrDOF,-rangeLeg2));
    putPtree(*leg2,"limit.upper",Cold::Constant(nrDOF,rangeLeg2));
    //putPtree(*leg2,"limit.coef",Cold::Constant(nrDOF,1000));  //use global joint limit
    put<scalarF>(*leg2,"geom1D.lenY",footLen);
    put<scalarF>(*leg2,"geom1D.rad",footRad);
    put<scalarF>(*leg2,"transY",1);
    putPtree<Vec3>(*leg2,"jointDirZ",R0*Vec3(0,0,1));
    putPtree<Vec3>(*leg2,"meshDirY",R0*R2*Vec3::UnitX());
  }
  return &pt;
}
tinyxml2::XMLElement* ArticulatedLoader::createArm(tinyxml2::XMLElement& pt,scalar armLen,scalar armRad)
{
  tinyxml2::XMLElement* root=addChild(pt,"joint");
  put<std::string>(*root,"type","FIX_JOINT");
  {
    //the table
    tinyxml2::XMLElement* child=addChild(*root,"joint");
    put<std::string>(*child,"type","FIX_JOINT");
    put<scalarF>(*child,"geom3D.lenX",armRad*3);
    put<scalarF>(*child,"geom3D.lenY",armRad);
    put<scalarF>(*child,"geom3D.lenZ",armRad*3);
    put<scalarF>(*child,"geom3D.rad",0);
    putPtree<Vec3>(*child,"trans",-Vec3::UnitY()*armRad*1.5f);
    putPtree<Vec3>(*child,"meshDirX",Vec3::UnitX());
    putPtree<Vec3>(*child,"meshDirY",Vec3::UnitY());

    tinyxml2::XMLElement* child2=addChild(*root,"joint");
    put<std::string>(*child2,"type","FIX_JOINT");
    put<scalarF>(*child2,"geom3D.lenX",armLen*10);
    put<scalarF>(*child2,"geom3D.lenY",armRad);
    put<scalarF>(*child2,"geom3D.lenZ",armLen*10);
    put<scalarF>(*child2,"geom3D.rad",0);
    putPtree<Vec3>(*child2,"trans",-Vec3::UnitY()*armRad*2.5f);
    putPtree<Vec3>(*child2,"meshDirX",Vec3::UnitX());
    putPtree<Vec3>(*child2,"meshDirY",Vec3::UnitY());
  }
  {
    //the arm
    tinyxml2::XMLElement* child=addChild(*root,"joint");
    put<std::string>(*child,"type","BALL_JOINT");
    putPtree(*child,"limit.lower",Vec2d(-D2R(180.0f),D2R(10.0f)));
    putPtree(*child,"limit.upper",Vec2d(D2R(180.0f),D2R(90.0f)));
    putPtree<Vec3>(*child,"jointDirY",Vec3(0,1,0));
    putPtree<Vec3>(*child,"jointDirZ",Vec3(0,0,1));
    put<scalarF>(*child,"geom1D.lenY",armLen);
    put<scalarF>(*child,"geom1D.rad",armRad);
    putPtree<Vec3>(*child,"meshDirY",Vec3::UnitY());

    tinyxml2::XMLElement* child2=addChild(*child,"joint");
    put<std::string>(*child2,"type","HINGE_JOINT");
    putPtree(*child2,"limit.lower",Cold::Constant(1,D2R(0.0f)));
    putPtree(*child2,"limit.upper",Cold::Constant(1,D2R(150.0f)));
    put<scalarF>(*child2,"transY",1);
    putPtree<Vec3>(*child2,"jointDirZ",Vec3(0,0,1));
    put<scalarF>(*child2,"geom1D.lenY",armLen);
    put<scalarF>(*child2,"geom1D.rad",armRad);
    putPtree<Vec3>(*child2,"meshDirY",Vec3::UnitY());
  }
  return &pt;
}
//build-in bodies as ArticulatedBody
ArticulatedBody ArticulatedLoader::createBird(sizeType root,bool head)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  createBird(*(pt.RootElement()),root,0.25f,0.5f,0.1f,0.6f, 0.3f,0.4f,0.45f, 0.12f,0.1f,0.2f,false,head);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createBipedal(sizeType root,bool withHand)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  tinyxml2::XMLElement* bodyNode=createBipedal(*(pt.RootElement()),root,0.25f,0.5f, 0.4f,0.4f,0.3f, 0.12f,0.1f,0.15f, false);
  if(withHand)
    addBipedalHand(*bodyNode,0.25f,0.5f, 0.4f,0.4f, 0.1f,0.1f);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createChain(sizeType root,scalar rad,sizeType nrLink,sizeType res)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  createChain(*(pt.RootElement()),root,nrLink,0.5f,rad,D2R(120),0,0,0,3,0.0f,0.2f,0.2f);
  put<sizeType>(pt,"globalRes",res);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createSpider(sizeType root,scalar lmt,scalar lmt2)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  scalar footLen=0.2f*sqrt(2.0f)+0.16f;
  ArticulatedLoader::createSpider(*(pt.RootElement()),root,0.25f,footLen,0.08f,D2R(10),D2R(lmt),D2R(lmt2),true);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createArm(scalar armLen,scalar armRad)
{
  tinyxml2::XMLDocument pt;
  pt.InsertEndChild(pt.NewElement("root"));
  createArm(*(pt.RootElement()),armLen,armRad);
  ArticulatedBody body;
  ArticulatedUtils utils(body);
  utils.assemble(*(pt.RootElement()));
  return body;
}
ArticulatedBody ArticulatedLoader::createMesh(const ObjMesh& m)
{
  Joint j;
  j._parent=-1;
  j._depth=0;
  j._typeJoint=Joint::FIX_JOINT;
  j._mimic=-1;
  j._offDOF=0;
  j._offDDT=0;
  j._limits.resize(3,0);
  j._control.resize(0);
  j._damping.resize(0);
  j._trans.setIdentity();
  //mimic
  j._mult=j._offset=0;
  //sphere approx.
  j._name="mesh";
  j._spheres.resize(3,0);
  j._radGeomColl.resize(0);
  j._radSelfColl.resize(0);
  j._color.setZero();
  //mass
  j._mesh.reset(new ObjMeshGeomCell(Mat4::Identity(),m,0,true));
  j.assemble(1);
  //body
  ArticulatedBody b;
  b._joints.push_back(j);
  b._geom.reset(new StaticGeom(3));
  b._geom->addGeomCell(j._mesh);
  b._geom->assemble();
  return b;
}
