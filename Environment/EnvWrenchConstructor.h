#ifdef ENVIRONMENT_SUPPORT
#ifndef ENV_WRENCH_CONSTRUCTOR_H
#define ENV_WRENCH_CONSTRUCTOR_H

#include "Environment.h"
#include <Articulated/SimplifiedDynamics.h>
#include <Utils/SparseUtils.h>
#include <Utils/Scalar.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T,int DIM>
class PBDDynamicsSequence;
template <typename T>
struct PBDArticulatedGradientInfo;
template <typename T>
struct NEArticulatedGradientInfoMap;
template <typename T>
struct NEArticulatedGradientInfo;
template <typename T>
struct ExternalWrench
{
  DECL_MAP_TYPES_T
  ExternalWrench();
  ExternalWrench& operator*=(T dt) {
    for(sizeType r=0; r<(sizeType)_DBDq.size(); r++)
      _DBDq[r].second*=dt;
    for(sizeType r=0; r<(sizeType)_DBDX.size(); r++)
      for(sizeType c=0; c<(sizeType)_DBDX[r].size(); c++)
        _DBDX[r][c]*=dt;
    _B*=dt;
    return *this;
  }
  template <typename T2>
  ExternalWrench<T2> cast() const {
    T2 tmp;
    ExternalWrench<T2> ret;
    ret._DBDq.resize(_DBDq.size());
    for(sizeType r=0; r<(sizeType)_DBDq.size(); r++)
      ret._DBDq.push_back(std::make_pair(_DBDq[r].first,_DBDq[r].second.unaryExpr([&](T in) {
      std::convert_scalar(in,tmp);
      return tmp;
    })));
    ret._DBDX.resize(_DBDX.size());
    for(sizeType r=0; r<(sizeType)_DBDX.size(); r++) {
      ret._DBDX[r].resize(_DBDX[r].size());
      for(sizeType c=0; c<(sizeType)_DBDX[r].size(); c++)
        ret._DBDX[r][c]=_DBDX[r][c].unaryExpr([&](T in) {
        std::convert_scalar(in,tmp);
        return tmp;
      });
    }
    ret._jointId=_jointId;
    ret._B=_B.unaryExpr([&](T in) {
      std::convert_scalar(in,tmp);
      return tmp;
    });
    ret._w=_w.unaryExpr([&](T in) {
      std::convert_scalar(in,tmp);
      return tmp;
    });
    return ret;
  }
  std::vector<std::pair<sizeType,Mat6XT>> _DBDq;
  std::vector<std::vector<Mat6XT,Eigen::aligned_allocator<Mat6XT>>> _DBDX;
  sizeType _jointId;
  Mat6XT _B;
  Vec _w;
};
template <typename T>
class EnvWrenchConstructor
{
public:
  DECL_MAP_TYPES_T
  EnvWrenchConstructor(std::shared_ptr<Environment<T>> env);
#ifdef TRAJ_OPT_SUPPORT
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,sizeType id,const PBDDynamicsSequence<T,3>& dyn,bool grad=true);
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,sizeType id,const PBDDynamicsSequence<T,2>& dyn,bool grad=true);
#endif
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,const PBDArticulatedGradientInfo<T>& info,bool grad=true);
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,const NEArticulatedGradientInfo<T>& info,bool grad=true);
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool grad=true);
  std::shared_ptr<Environment<T>> _env;
};
template <typename T,typename T2>
class EnvWrenchAdaptor : public EnvWrenchConstructor<T>
{
public:
  DECL_MAP_TYPES_T
  using EnvWrenchConstructor<T>::operator();
  EnvWrenchAdaptor(std::shared_ptr<EnvWrenchConstructor<T2>> env):_env(env) {}
#ifdef TRAJ_OPT_SUPPORT
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,sizeType id,const PBDDynamicsSequence<T,3>& dyn,bool grad=true) override {
    T2 tmp;
    std::vector<ExternalWrench<T2>> externalWrenchT2;
    _env->operator()(externalWrenchT2,[&](sizeType JID)->typename EnvWrenchConstructor<T2>::Mat3X4T{
      return dyn.getJointTrans(id,JID).unaryExpr([&](T in) {
        std::convert_scalar(in,tmp);
        return tmp;
      });
    },grad);
    externalWrench.clear();
    for(auto ew:externalWrenchT2)
      externalWrench.push_back(ew.template cast<T>());
  }
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,sizeType id,const PBDDynamicsSequence<T,2>& dyn,bool grad=true) override {
    T2 tmp;
    std::vector<ExternalWrench<T2>> externalWrenchT2;
    _env->operator()(externalWrenchT2,[&](sizeType JID)->typename EnvWrenchConstructor<T2>::Mat3X4T{
      return dyn.getJointTrans(id,JID).unaryExpr([&](T in) {
        std::convert_scalar(in,tmp);
        return tmp;
      });
    },grad);
    externalWrench.clear();
    for(auto ew:externalWrenchT2)
      externalWrench.push_back(ew.template cast<T>());
  }
#endif
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,const PBDArticulatedGradientInfo<T>& info,bool grad=true) override {
    T2 tmp;
    std::vector<ExternalWrench<T2>> externalWrenchT2;
    _env->operator()(externalWrenchT2,[&](sizeType JID)->typename EnvWrenchConstructor<T2>::Mat3X4T{
      return TRANSI(info._TM,JID).unaryExpr([&](T in) {
        std::convert_scalar(in,tmp);
        return tmp;
      });
    },grad);
    externalWrench.clear();
    for(auto ew:externalWrenchT2)
      externalWrench.push_back(ew.template cast<T>());
  }
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,const NEArticulatedGradientInfo<T>& info,bool grad=true) override {
    T2 tmp;
    std::vector<ExternalWrench<T2>> externalWrenchT2;
    _env->operator()(externalWrenchT2,[&](sizeType JID)->typename EnvWrenchConstructor<T2>::Mat3X4T{
      return info.NEArticulatedGradientInfoMap<T>::getTrans(JID).unaryExpr([&](T in) {
        std::convert_scalar(in,tmp);
        return tmp;
      });
    },grad);
    externalWrench.clear();
    for(auto ew:externalWrenchT2)
      externalWrench.push_back(ew.template cast<T>());
  }
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool grad=true) override {
    T2 tmp;
    std::vector<ExternalWrench<T2>> externalWrenchT2;
    _env->operator()(externalWrenchT2,[&](sizeType JID)->typename EnvWrenchConstructor<T2>::Mat3X4T{
      return tfunc(JID).unaryExpr([&](T in) {
        std::convert_scalar(in,tmp);
        return tmp;
      });
    },grad);
    externalWrench.clear();
    for(auto ew:externalWrenchT2)
      externalWrench.push_back(ew.template cast<T>());
  }
  std::shared_ptr<EnvWrenchConstructor<T2>> _env;
};
template <typename T>
class C0EnvWrenchConstructor : public EnvWrenchConstructor<T>
{
public:
  DECL_MAP_TYPES_T
  using EnvWrenchConstructor<T>::_env;
  using EnvWrenchConstructor<T>::operator();
  C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000);
  C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000);
  C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000);
  C0EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD z,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000);
  C0EnvWrenchConstructor(const ArticulatedBody& body,const Vec4d& plane,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000);
  C0EnvWrenchConstructor(const ArticulatedBody& body,const std::string& path,bool is2D,scalarD dxMul,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000);
  C0EnvWrenchConstructor(const ArticulatedBody& body,std::shared_ptr<Environment<T>> env,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000);
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool grad=true) override;
  static Mat3T rotateVec(const Vec3T a,const Vec3T& b,Mat3T* dRdX=NULL,Mat3T* dRdY=NULL,Mat3T* dRdZ=NULL);
  static void debugRotateVec(sizeType nrIter=10);
  void writeEndEffectorVTK(const std::string& path) const;
  void writeVTK(const std::string& path) const;
  std::vector<EndEffectorBounds> _externalForces;
  const ArticulatedBody& _body;
  Mat3T _frm;
  Mat3XT _B;
};
template <typename T>
class C2EnvWrenchConstructor : public C0EnvWrenchConstructor<T>
{
public:
  DECL_MAP_TYPES_T
  using C0EnvWrenchConstructor<T>::_body;
  using C0EnvWrenchConstructor<T>::_env;
  using C0EnvWrenchConstructor<T>::_frm;
  using C0EnvWrenchConstructor<T>::_B;
  using C0EnvWrenchConstructor<T>::rotateVec;
  using C0EnvWrenchConstructor<T>::writeEndEffectorVTK;
  using C0EnvWrenchConstructor<T>::writeVTK;
  using C0EnvWrenchConstructor<T>::_externalForces;
  using C0EnvWrenchConstructor<T>::operator();
  C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD x0,scalarD z0,scalarD slope,sizeType n,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000,T relax=0);
  C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,std::function<scalarD(scalarD,scalarD)> h,sizeType res,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000,T relax=0);
  C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD wid,scalarD len,scalarD off,sizeType n,scalarD z,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000,T relax=0);
  C2EnvWrenchConstructor(const ArticulatedBody& body,scalarD x,scalarD y,scalarD z,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000,T relax=0);
  C2EnvWrenchConstructor(const ArticulatedBody& body,const Vec4d& plane,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000,T relax=0);
  C2EnvWrenchConstructor(const ArticulatedBody& body,const std::string& path,bool is2D,scalarD dxMul,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000,T relax=0);
  C2EnvWrenchConstructor(const ArticulatedBody& body,std::shared_ptr<Environment<T>> env,const Vec3T& g,sizeType nrDir=6,T mu=0.7f,T coef=1000,T relax=0);
  virtual void operator()(std::vector<ExternalWrench<T>>& externalWrench,std::function<Mat3X4T(sizeType)> tfunc,bool grad=true) override;
  Mat3XT _crossB;
  T _coef,_relax;
};

PRJ_END

#endif
#endif
