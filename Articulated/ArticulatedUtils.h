#ifndef ARTICULATED_UTILS_H
#define ARTICULATED_UTILS_H

#include "ArticulatedBody.h"
#include <tinyxml2.h>

PRJ_BEGIN

class ArticulatedUtils
{
public:
  struct TransInfo {
    Mat3X4d _TLast,_X0;
    Vec3d _lenX,_lenY,_lenZ,_origin;
    scalar _rad,_rad2,_rad3,_rho;
  };
  ArticulatedUtils(ArticulatedBody& body);
  void assembleGlobalVars(const tinyxml2::XMLElement& pt);
  void assemble(const tinyxml2::XMLElement& pt);
  void sampleGeom3DBox(Joint& J) const;
  void sampleGeom2Sphere(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const;
  void sampleGeom3Sphere(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const;
  void sampleGeom2D(Joint& J,const Vec3d& base,const Vec3d& lenX,const Vec3d& lenY,scalarD rad,const tinyxml2::XMLElement& pt) const;
  void sampleGeom3D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const;
  void sampleGeom2D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const;
  void sampleGeom1D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const;
  void sampleGeom0D(Joint& J,const TransInfo& I,const tinyxml2::XMLElement& pt) const;
  void assembleJoints(const tinyxml2::XMLElement& pt,std::vector<TransInfo>& infos);
  void assembleJoint(const tinyxml2::XMLElement& pt,sizeType parent,std::vector<TransInfo>& infos);
  void addBase(sizeType dim,const Vec3d& planeNormal);
  void combine(const std::vector<ArticulatedBody>& bodies);
  void simplify(std::function<bool(const Joint&)> canSimplify,const Cold& DOF,sizeType nrDebug);
  void simplify(sizeType nrDebug);
  void eliminateJoint(const std::vector<std::string>& jointName,const Cold& DOF,sizeType nrDebug);
  void scaleMass(scalarD coef);
  scalarD totalMass() const;
private:
  ArticulatedBody& _body;
  //transient data, no included in read/write
  ALIGN_16 scalarD _globalRho;
  ALIGN_16 scalarD _globalLimitCoef;
  ALIGN_16 scalarD _globalControlCoef;
  ALIGN_16 scalarD _globalDampingCoef;
  ALIGN_16 scalarD _globalSampleDist;
  ALIGN_16 scalarD _globalSampleDistLocal;
  ALIGN_16 sizeType _globalRes;
  ALIGN_16 Vec3 _globalColor;
};

PRJ_END

#endif
