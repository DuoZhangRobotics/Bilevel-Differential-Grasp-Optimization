#ifdef ENVIRONMENT_SUPPORT
#ifndef PBD_SIMULATOR_H
#define PBD_SIMULATOR_H

#include "PBDArticulatedGradientInfo.h"
#include <Environment/EnvWrenchConstructor.h>
#include <Optimizer/DSSQPObjective.h>
#include <Utils/DebugGradient.h>
#include <Utils/Options.h>

PRJ_BEGIN

enum PBD_SIMULATOR_MODE {
  PBD_NO_FORCE,
  PBTO,
  NPBD_SQP,
  NPBD_PGM,
  NPBD_ZOPGM,
};
template <typename T>
struct PBDSimulatorTraits;
template <typename T>
class PBDNMDPObjective;
struct PDTarget;
template <typename T>
class PBDSimulator
#ifdef OPTIMIZER_SUPPORT
  : public DSSQPObjective<T>
#endif
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  friend class PBDNMDPObjective<T>;
  friend struct PBDSimulatorTraits<T>;
  typedef ParallelVector<Eigen::Triplet<T,sizeType> > TRIPS;
  PBDSimulator(const ArticulatedBody& body,Options& ops,const Vec3T& g,PBD_SIMULATOR_MODE integrator=PBTO);
  PBDSimulator(const PBDSimulator<T>& other);
  PBDSimulator<T>& operator=(const PBDSimulator<T>& other);
  void setWrenchConstructor(std::shared_ptr<C0EnvWrenchConstructor<T>> wrench);
  void setPDTarget(std::shared_ptr<PDTarget> PDTarget=NULL);
  void writeVTKSeq(const std::string& path,T horizon,T interval,T dt,std::shared_ptr<PDTarget> PDTarget=NULL,const std::map<std::string,std::set<sizeType>>* jointMask=NULL);
  void writeVTK(const Vec& q,const std::string& path,const std::map<std::string,std::set<sizeType>>* jointMask=NULL) const;
  void setState(const Vec& qM,const Vec& qMM);
  bool step(T dt,T& lastDt,const Vec* tau0,Vec& qqMqMM,T* betaInit=NULL);
  bool step(T dt,T& lastDt,VecCM tau0,const VecM qqMqMM,T* betaInit=NULL);
  bool eval(const PBDArticulatedGradientInfo<T>& x,const PBDArticulatedGradientInfo<T>& xM,const PBDArticulatedGradientInfo<T>& xMM,const Vec& ctrl,T dt,T lastDt,T* E,Vec* G,TRIPS* H,sizeType off) const;
  void debug(T dt,sizeType nrIter,T scale=1000);
  void reset(Options& ops);
  void setMode(PBD_SIMULATOR_MODE m);
  PBD_SIMULATOR_MODE getMode() const;
private:
  Vec stepPBTO(T dt,T lastDt,VecCM tau0,VecCM qMqMM);
  Vec stepNMDPSQP(T dt,T lastDt,VecCM tau0,VecCM qMqMM);
  Vec stepNMDPPGM(T dt,T lastDt,VecCM tau0,VecCM qMqMM,T* betaInit=NULL);
  bool stepNMDPPGMAdaptive(T dt,T& lastDt,VecCM tau0,VecM qMqMM,T* betaInit=NULL);
  Vec manifoldProjection(T dt,T lastDt,VecCM tau0,VecCM init,Vec* w,T* alphaInit=NULL);
  //helper
  T E(T dt,T lastDt,VecCM tau0) const;
  T E(T dt,T lastDt,VecCM tau0,
      const PBDArticulatedGradientInfo<T>& q,
      const PBDArticulatedGradientInfo<T>& qM,
      const PBDArticulatedGradientInfo<T>& qMM) const;
  T K(T dt,VecM DKDTheta,MatTM DDKDDTheta) const;
  Vec G(T dt,T lastDt,VecCM tau0,MatTM DGDTheta,Vec* w,MatTM DGDw);
  Vec G(T dt,T lastDt,VecCM tau0,MatTM DGDTheta,Vec* w,MatTM DGDw,
        const PBDArticulatedGradientInfo<T>& q,
        const PBDArticulatedGradientInfo<T>& qM,
        const PBDArticulatedGradientInfo<T>& qMM,
        std::vector<ExternalWrench<T>>& externalWrench,sizeType& iterG) const;
  Vec3T forcePBTO(T dt,T phi0,const Vec3T qqM[2],Mat3T* DFDPos,Mat3T* DFDPosN) const;
  Vec3T forcePBTO(T dt,const EndEffectorBounds& EE,Mat3T* DFDPos,Mat3T* DFDPosN) const;
  void DFDTheta(MatTM DGDTheta,const ExternalWrench<T>& ew,const EndEffectorBounds& ee,const SMat& dPos,const Vec& w) const;
  SMat DPos(const EndEffectorBounds& EE) const;
  sizeType nW() const;
  //data
  std::shared_ptr<PDTarget> _PDTarget;
  std::shared_ptr<C0EnvWrenchConstructor<T>> _wrench;
  std::vector<ExternalWrench<T>> _externalWrench;
  const ArticulatedBody& _body;
  PBD_SIMULATOR_MODE _mode;
  T _alphaPBTO,_betaPBTO;
  sizeType _maxIter;
  T _tolG,_tolLS,_tolK,_minDt;
  T _betaMin,_alphaBetaInc;
  bool _callback,_useKEE;
  Mat3XT _JRCF;
  //transient data
  PBDArticulatedGradientInfo<T> _qqMqMM[3];
  sizeType _iter,_iterG;
};
template <>
struct PBDSimulatorTraits<double>
{
  typedef double T;
  DECL_MAP_TYPES_T
  static void registerOptions(Options& ops);
  static void initOptions(PBDSimulator<double>& sol);
};
template <>
struct PBDSimulatorTraits<__float128>
{
  typedef __float128 T;
  DECL_MAP_TYPES_T
  static void registerOptions(Options& ops);
  static void initOptions(PBDSimulator<__float128>& sol);
};
template <>
struct PBDSimulatorTraits<mpfr::mpreal>
{
  typedef mpfr::mpreal T;
  DECL_MAP_TYPES_T
  static void registerOptions(Options& ops);
  static void initOptions(PBDSimulator<mpfr::mpreal>& sol);
};
#ifdef OPTIMIZER_SUPPORT
template <typename T>
class PBDNMDPObjective : public DSSQPObjective<T>
{
public:
  DECL_MAP_TYPES_T
  DECL_MAP_FUNCS
  PBDNMDPObjective(T dt,T lastDt,VecCM tau0,PBDSimulator<T>& sim);
  //objective
  virtual Vec lb() const;
  virtual Vec ub() const;
  virtual Vec gl() const;
  virtual Vec gu() const;
  virtual Vec init() const;
  //constraint
  virtual int operator()(const Vec& x,Vec& fvec,STrips* fjac);
  //objective
  virtual T operator()(const Vec& x,Vec* fgrad);
  //problem size
  virtual int inputs() const;
  virtual int values() const;
protected:
  T _dt,_lastDt;
  VecCM _tau0;
  PBDSimulator<T>& _sim;
};
#endif

PRJ_END

#endif
#endif
