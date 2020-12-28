#ifndef GRANULAR_WRENCH_CONSTRUCTOR_H
#define GRANULAR_WRENCH_CONSTRUCTOR_H

#include "EnvWrenchConstructor.h"

PRJ_BEGIN

template <typename T>
class C2GranularWrenchConstructor : public EnvWrenchConstructor<T>
{
public:
  DECL_MAP_TYPES_T
  using EnvWrenchConstructor<T>::_env;
  using EnvWrenchConstructor<T>::operator();
  typedef std::function<T(const Vec2T&,Vec2T* diff)> RBF;
  C2GranularWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n,
                              const std::string& pathInput="inputs.txt",const std::string& pathWeight="weights.txt");
  C2GranularWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res,
                              const std::string& pathInput="inputs.txt",const std::string& pathWeight="weights.txt");
  C2GranularWrenchConstructor(const ArticulatedBody& body,scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z,
                              const std::string& pathInput="inputs.txt",const std::string& pathWeight="weights.txt");
  C2GranularWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD z,
                              const std::string& pathInput="inputs.txt",const std::string& pathWeight="weights.txt");
  C2GranularWrenchConstructor(const ArticulatedBody& body,const Vec4d& plane,
                              const std::string& pathInput="inputs.txt",const std::string& pathWeight="weights.txt");
  C2GranularWrenchConstructor(const ArticulatedBody& body,const std::string& path,bool is2D,scalarD dxMul,
                              const std::string& pathInput="inputs.txt",const std::string& pathWeight="weights.txt");
  C2GranularWrenchConstructor(const ArticulatedBody& body,std::shared_ptr<Environment<T>> env,
                              const std::string& pathInput="inputs.txt",const std::string& pathWeight="weights.txt");
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool grad=true) override;
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,RBF kernel,bool grad);
  void writeEndEffectorVTK(const std::string& path,bool torqueCenter) const;
  void writeVTK(const std::string& path) const;
  static Mat3T rotateVecPlanar(const Vec3T d0,const Vec3T& d1,Mat3T dR0[3]=NULL,Mat3T dR1[3]=NULL);
  static T theta(const Vec3T& d0,const Vec3T& d1,Vec3T* dT0=NULL,Vec3T* dT1=NULL);
  static void debugTheta(sizeType nrIter=10);
  std::vector<EndEffectorBounds> _externalForces;
protected:
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> _torqueCenters;
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> _normalDirs;
  const ArticulatedBody& _body;
  Mat3XT _weight;
  Mat2XT _input;
};

PRJ_END

#endif
