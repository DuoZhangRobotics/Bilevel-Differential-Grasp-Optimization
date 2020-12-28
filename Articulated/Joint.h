#ifndef JOINT_H
#define JOINT_H

#include <CommonFile/ObjMesh.h>
#include <CommonFile/StaticPool.h>
#include <Utils/ArticulatedBodyPragma.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

class StaticGeom;
struct StaticGeomCell;
struct ArticulatedBody;
struct Joint : public SerializableBase
{
  typedef scalarD T;
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  friend struct ArticulatedBody;
  friend class ArticulatedLoader;
  friend class ArticulatedUtils;
  using SerializableBase::read;
  using SerializableBase::write;
  enum JOINT_TYPE {
    TRANS_3D    =1<<0,
    TRANS_2D    =1<<1,
    TRANS_1D    =1<<2,
    ROT_3D_XYZ  =1<<3,
    ROT_3D_EXP  =1<<4,
    BALL_JOINT  =1<<5,
    HINGE_JOINT =1<<6,
    FIX_JOINT   =1<<7,
    NR_JOINT_TYPE,
  };
  enum GEOM_TYPE {
    MESH        =1<<0,
    SPHERE_GEOM =1<<1,
    SPHERE_SELF =1<<2,
    AXIS_INFO   =1<<4,
  };
  Joint();
  bool read(std::istream& is,IOData* dat) override;
  bool write(std::ostream& os,IOData* dat) const override;
  std::shared_ptr<SerializableBase> copy() const override;
  std::string type() const override;
  void writeVTK(const Mat3X4d& T,VTKWriter<scalar>& os,std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& css,GEOM_TYPE type,const std::vector<unsigned char>* mask) const;
  ObjMesh getMesh(GEOM_TYPE type,bool render=false) const;
  void setType(const std::string& type);
  static std::string typeToString(JOINT_TYPE type);
  void assemble(scalar rho);
  void transformMesh(const Mat3X4d& T);
  void transformMesh(const Mat3d& R,const Vec3d& X);
  void transformMesh(const Vec3d& X);
  Vec2i CBegEnd() const;
  Vec2i RBegEnd() const;
  sizeType nrDOF() const;
  sizeType nrDDT() const;
  Mat6d getMassC(const Mat3d& R) const;
  Mat6d getMassC() const;
  Mat6d getMass(const Mat3d& R) const;
  Mat6d getMass() const;
  Vec3d getC() const;
  template <sizeType R,sizeType C>
  scalarD getMEntry() const
  {
    if(R==3 && C==3)
      return _M;
    else if(R==3)
      return _MC[C];
    else if(C==3)
      return _MC[R];
    return _MCCT(R,C);
  }
  bool isRotational() const;
  bool isRoot(const ArticulatedBody& body) const;
  Mat3Xd getSpheres(const Mat3X4d& T) const;
  void getSpheres(const Mat3X4d& T,Mat3XTM Sss) const;
  std::shared_ptr<StaticGeomCell> getGeomPtr() const;
  scalarD avgRadGeomColl() const;
  scalarD avgRadSelfColl() const;
  static void loopAllJointTypes(std::function<void(Joint::JOINT_TYPE)> t);
  //joint
  std::vector<sizeType> _children;
  sizeType _parent,_depth,_typeJoint,_mimic;
  sizeType _offDOF,_offDDT;
  Mat3Xd _limits;
  Cold _control;
  Cold _damping;
  Mat3X4d _trans;
  //mimic
  scalarD _mult,_offset;
  //mass
  scalarD _M;
  Vec3d _MC;
  Mat3d _MCCT;
  //sphere approx.
  std::string _name;
  Mat3Xd _spheres;
  Cold _radGeomColl,_radSelfColl;
  Vec3 _color;
private:
  //mesh approx.
  std::shared_ptr<StaticGeomCell> _mesh;
};

PRJ_END

#endif
