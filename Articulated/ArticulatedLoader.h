#ifndef ARTICULATED_LOADER_H
#define ARTICULATED_LOADER_H

#include "ArticulatedBody.h"
#include <experimental/filesystem>
#include <tinyxml2.h>
struct aiScene;
struct aiNode;

PRJ_BEGIN

class ArticulatedLoader
{
public:
  typedef scalarD T;
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  typedef std::unordered_map<std::string,Joint> Links;
  struct JointInfo {
    std::string _child;
    std::string _name;
    std::string _type;
    Mat4 _T;
    Vec3 _axis;
    Vec2 _lmt;
    scalar _ctrl;
    scalar _damp;
    //mimic
    std::string _mimic;
    scalar _multiplier;
    scalar _offset;
  };
  typedef std::unordered_map<std::string,std::vector<JointInfo>> Relations;
  static ObjMesh createMeshFromScene(const aiScene* scene,const aiNode* node);
  static Mat3 RPY2Mat(const Vec3& v);
  static void mergeMesh(ObjMesh& mesh);
  static void makeConvex(tinyxml2::XMLElement& pt);
  static ObjMesh readMesh(const std::experimental::filesystem::v1::path& path,bool adjust);
  static void readLinks(Links& links,const tinyxml2::XMLElement& pt,const std::experimental::filesystem::v1::path& path,bool visualMesh);
  static bool readRelations(Relations& relations,const tinyxml2::XMLElement& pt);
  static std::string findRelationByChild(const std::string& name,const Relations& relations,JointInfo& info);
  static std::string findRelationByName(const std::string& name,const Relations& relations,JointInfo& info);
  static bool buildBody(ArticulatedBody& body,const std::string& name,const Links& links,const Relations& relations);
  static ArticulatedBody readURDF(const std::string& file,bool convex,bool visualMesh);
  static void compareVisualMeshURDF(const std::string& file);
  //visualize
  static void visualizeJointLimitVTK(const ArticulatedBody& body,const std::string& path,Joint::GEOM_TYPE type);
  //built-in bodies
  static sizeType inferDim(sizeType rootType);
  static tinyxml2::XMLElement* addRootJoint(tinyxml2::XMLElement& pt,sizeType rootType);
  static tinyxml2::XMLElement* createBird(tinyxml2::XMLElement& pt,sizeType rootType,scalar bodySz,scalar bodyLen,scalar neckSz,scalar neckLen,scalar footLen1,scalar footLen2,scalar footLen3,scalar footRad1,scalar footRad2,scalar footRad3,bool fixFoot,bool head);
  static tinyxml2::XMLElement* createBipedal(tinyxml2::XMLElement& pt,sizeType rootType,scalar bodySz,scalar bodyLen,scalar footLen1,scalar footLen2,scalar footLen3,scalar footRad1,scalar footRad2,scalar footRad3,bool fixFoot);
  static void addBipedalHand(tinyxml2::XMLElement& pt,scalar bodySz,scalar bodyLen,scalar handLen1,scalar handLen2,scalar handRad1,scalar handRad2);
  static tinyxml2::XMLElement* createChain(tinyxml2::XMLElement& pt,sizeType rootType,sizeType nr,scalar l,scalar rad,scalar rot,scalar rot0,scalar t,scalar t0,sizeType geomDim,scalar yOff=0,scalar ratioX=1,scalar ratioZ=1);
  static tinyxml2::XMLElement* createSpider(tinyxml2::XMLElement& pt,sizeType rootType,scalar bodySz,scalar footLen,scalar footRad,scalar angLegLift,scalar rangeLeg1,scalar rangeLeg2,bool ball);
  static tinyxml2::XMLElement* createArm(tinyxml2::XMLElement& pt,scalar armLen,scalar armRad);
  //build-in bodies as ArticulatedBody
  static ArticulatedBody createBird(sizeType root,bool head=false);
  static ArticulatedBody createBipedal(sizeType root,bool withHand);
  static ArticulatedBody createChain(sizeType root,scalar rad,sizeType nrLink,sizeType res);
  static ArticulatedBody createSpider(sizeType root,scalar lmt=45,scalar lmt2=60);
  static ArticulatedBody createArm(scalar armLen,scalar armRad);
  static ArticulatedBody createMesh(const ObjMesh& m);
};

PRJ_END

#endif
