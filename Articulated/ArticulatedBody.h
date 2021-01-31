#ifndef ARTICULATED_BODY_H
#define ARTICULATED_BODY_H

#include <set>
#include "Joint.h"
#include <tinyxml2.h>

PRJ_BEGIN

struct ArticulatedBody : public SerializableBase
{
  typedef scalarD T;
  DECL_MAP_TYPES_T
  using SerializableBase::read;
  using SerializableBase::write;
  friend class ArticulatedLoader;
  friend class ArticulatedUtils;
  struct TransInfo {
    Mat3X4d _TLast,_X0;
    Vec3d _lenX,_lenY,_lenZ,_origin;
    scalar _rad,_rad2,_rad3,_rho;
  };
  ArticulatedBody();
  ArticulatedBody(const tinyxml2::XMLElement& pt);
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  void beginUpdateGeom(const Mat3Xd& T,std::vector<Mat4,Eigen::aligned_allocator<Mat4>>& tss);
  void endUpdateGeom(const std::vector<Mat4,Eigen::aligned_allocator<Mat4>>& tss);
  const StaticGeom& getGeom() const;
  StaticGeom& getGeom();
  void randomize(sizeType nrLink,bool chain=true);
  ObjMesh writeMesh(const Mat3Xd& T,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask=NULL) const;
  void writeVTK(const Mat3Xd& T,VTKWriter<scalar>& os,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask=NULL,const std::vector<std::vector<unsigned char> >* sphereMask=NULL) const;
  void writeVTK(const Mat3Xd& T,const std::string& path,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask=NULL,const std::vector<std::vector<unsigned char> >* sphereMask=NULL) const;
  void writePov(const Mat3Xd& T,const std::string& path,Joint::GEOM_TYPE type,const std::set<sizeType>* jointMask=NULL) const;
  std::set<sizeType> children(sizeType id,bool direct=false) const;
  sizeType commonRoot(sizeType id,sizeType id2) const;
  sizeType hasJoint(Joint::JOINT_TYPE type) const;
  ArticulatedBody resizeJoints(sizeType nr) const;
  bool isLeaf(sizeType id) const;
  Cold control() const;
  Cold damping() const;
  Cold coefLimit() const;
  Cold lowerLimit() const;
  Cold upperLimit() const;
  Cold lowerLimit(scalarD infty) const;
  Cold upperLimit(scalarD infty) const;
  const Cold& clampLimit(Cold& x) const;
  Cold randomPose(scalarD coef) const;
  Mat3Xd getT(const Cold& x) const;
  void mimic(VecM xMap) const;
  void mimic(MatT& A,Vec& b,Vec& l,Vec& u) const;
  bool movable(sizeType id,sizeType croot=-1) const;
  const Joint& joint(sizeType id) const;
  Joint& joint(sizeType id);
  sizeType jointId(const std::string& name) const;
  sizeType rootJointId() const;
  sizeType depth() const;
  sizeType nrDOF() const;
  sizeType nrDDT() const;
  sizeType nrJ() const;
  void fillChildren();
  void setRootTrans(const Mat3X4d& t);
  void debugBase(const std::string& path,T coef) const;
  void addBase(sizeType dim,const Vec3d& planeNormal);
  void simplify(sizeType nrDebug);
  void eliminateJoint(const std::vector<std::string>& jointName,const Cold& DOF,sizeType nrDebug);
  void scaleMass(scalarD coef);
  scalarD totalMass() const;
protected:
  ALIGN_16 std::vector<Joint> _joints;
  ALIGN_16 std::shared_ptr<StaticGeom> _geom;
};

PRJ_END

#endif
