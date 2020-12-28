#ifdef OPTIMIZER_SUPPORT
#ifdef ENVIRONMENT_SUPPORT
#ifdef DEFORMABLE_SUPPORT
#include "Simulator.h"
#include <Deformable/FEMSystem.h>
#include <Articulated/JointFunc.h>
#include <Articulated/PBDSimulator.h>
#include <Utils/DebugGradient.h>

USE_PRJ_NAMESPACE

template <typename T>
Simulator<T>::Simulator():_dt(0.01f),_lastDt(0.01f) {}
template <typename T>
sizeType Simulator<T>::addFEMObject(std::shared_ptr<FEMSystem<T>> sys,const Vec& xInit)
{
  SimulatorObject obj;
  obj._sysFEM=sys;
  obj._xFEM[0].reset(*sys,mapV(xInit));
  obj._xFEM[2]=obj._xFEM[1]=obj._xFEM[0];
  _oss.push_back(obj);
  return (sizeType)_oss.size()-1;
}
template <typename T>
sizeType Simulator<T>::addPBDObject(std::shared_ptr<ArticulatedBody> body,const Vec3T& g,const Vec& xInit)
{
  SimulatorObject obj;
  obj._body=body;
  Options ops;
  obj._sysPBD.reset(new PBDSimulator<T>(*body,ops,g,PBTO));
  obj._xPBD[0].reset(*body,mapV(xInit));
  obj._xPBD[2]=obj._xPBD[1]=obj._xPBD[0];
  _oss.push_back(obj);
  return (sizeType)_oss.size()-1;
}
template <typename T>
void Simulator<T>::writeVTK(bool embed,const std::string& path) const
{
  ObjMesh mesh;
  writeObj(embed,mesh);
  mesh.writeVTK(path,true);
}
template <typename T>
void Simulator<T>::writeObj(bool embed,ObjMesh& mesh) const
{
  for(sizeType i=0; i<(sizeType)_oss.size(); i++) {
    ObjMesh m;
    if(_oss[i]._sysFEM) {
      _oss[i]._sysFEM->writeObj(embed,_oss[i]._xFEM[0],m);
    } else if(_oss[i]._sysPBD) {
      m=_oss[i]._body->writeMesh(_oss[i]._xPBD[0]._TM.unaryExpr([&](const T& in) {
        return (scalarD)std::to_double(in);
      }),Joint::MESH);
    } else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
    mesh.addMesh(m);
  }
}
template <typename T>
bool Simulator<T>::solve(const Vec& ctrl,T GThres,T alphaThres,sizeType itMax,bool direct,bool callback)
{
  Vec G,d;
  bool valid;
  T alpha=1.0f,E,ENew;
  std::vector<SimulatorState> state,stateNew;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<scalarD,0,sizeType>,Eigen::Lower|Eigen::Upper,Eigen::DiagonalPreconditioner<scalarD>> solverI;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<scalarD,0,sizeType>> solverD;
  solverI.setTolerance(1e-4f);
  for(sizeType it=0; alpha>alphaThres && it<itMax; it++) {
    Eigen::SparseMatrix<T,0,sizeType> H;
    valid=buildSystem(ctrl,&E,&G,&H);
    if(!valid) {
      if(callback) {
        std::cout << "solve iter=" << it << " E=" << E << " GVio=" << G.unaryExpr([&](const T& in) {
          return std::abs(in);
        }).maxCoeff() << " failed(invalid system)" << std::endl;
      }
      return false;
    }
    //termination
    if(G.unaryExpr([&](const T& in) {
    return std::abs(in);
    }).maxCoeff()<GThres) {
      if(callback) {
        std::cout << "solve iter=" << it << " E=" << E << " GVio=" << G.unaryExpr([&](const T& in) {
          return std::abs(in);
        }).maxCoeff() << " succeeded" << std::endl;
      }
      break;
    }
    //solve
    if(direct)
      solverD.compute(H.template cast<scalarD>());
    else solverI.compute(H.template cast<scalarD>());
    if((direct?solverD.info():solverI.info())!=Eigen::Success) {
      if(callback) {
        std::cout << "solve iter=" << it << " E=" << E << " GVio=" << G.unaryExpr([&](const T& in) {
          return std::abs(in);
        }).maxCoeff() << " alpha=" << alpha << " failed(factorization)" << std::endl;
      }
      break;
    }
    if(direct)
      d=-solverD.solve(G.template cast<scalarD>()).template cast<T>();
    else d=-solverI.solve(G.template cast<scalarD>()).template cast<T>();
    //line search
    saveState(state);
    while(alpha>alphaThres) {
      updateState(stateNew=state,d*alpha);
      loadState(stateNew);
      valid=buildSystem(ctrl,&ENew,NULL,NULL);
      //update line search
      if(valid && std::isfinite(ENew) && ENew<E) {
        alpha=std::min<T>(alpha*1.5f,1.f);
        if(callback) {
          std::cout << "solve iter=" << it << " E=" << E << " alpha=" << alpha << std::endl;
        }
        break;
      } else alpha*=0.5f;
    }
    if(alpha<alphaThres) {
      if(callback) {
        std::cout << "solve iter=" << it << " E=" << E << " alpha=" << alpha << " failed(small alpha)" << std::endl;
      }
    }
  }
  return G.unaryExpr([&](const T& in) {
    return std::abs(in);
  }).maxCoeff()<GThres;
}
template <typename T>
void Simulator<T>::attach(sizeType oidA,sizeType oidB)
{
  if(oidA<0||oidA>=(sizeType)_oss.size()) {
    WARNINGV("invalid object id %d",oidA)
    return;
  }
  if(oidB<0||oidB>=(sizeType)_oss.size()) {
    WARNINGV("invalid object id %d",oidB)
    return;
  }
  if(oidA==oidB) {
    WARNINGV("cannot attach object %d to itself",oidA)
    return;
  }
  std::unordered_set<sizeType> children=attachedChildren(oidA);
  if(children.find(oidB)!=children.end()) {
    WARNINGV("object %d is already attached to object %d",oidB,oidA)
    return;
  }
  children=attachedChildren(oidB);
  if(children.find(oidA)!=children.end()) {
    WARNINGV("object %d is already attached to object %d, cannot reverse the attachment",oidA,oidB)
    return;
  }
  Mat3X4T TA=getTrans(oidA),invTA;
  Mat3X4T TB=getTrans(oidB),invTATB;
  INV(invTA,TA)
  APPLY_TRANS(invTATB,invTA,TB)
  _oss[oidA]._attachedObjects[oidB]=invTATB;
}
template <typename T>
void Simulator<T>::setTrans(sizeType oid,sizeType id,const Mat3X4T& t)
{
  if(oid<0||oid>=(sizeType)_oss.size()) {
    WARNINGV("Invalid object id %d",oid)
    return;
  }
  SimulatorObject& o=_oss[oid];
  if(o._sysFEM) {
    o._xFEM[id].setT(*(o._sysFEM),t);
  } else if(o._sysPBD) {
    Vec x=o._xPBD[id]._xM;
    o._body->setRootTrans(t.unaryExpr([&](const T& in) {
      return (scalarD)std::to_double(in);
    }));
    o._xPBD[id].reset(*(o._body),x);
  } else {
    ASSERT_MSG(false,"Unknown dynamics type")
  }
  Mat3X4T tc;
  for(const std::pair<sizeType,Mat3X4T>& c:_oss[oid]._attachedObjects) {
    APPLY_TRANS(tc,t,c.second)
    setTrans(c.first,id,tc);
  }
}
template <typename T>
typename Simulator<T>::Mat3X4T Simulator<T>::getTrans(sizeType oid,sizeType id) const
{
  if(oid<0||oid>=(sizeType)_oss.size()) {
    WARNINGV("invalid object id %d",oid)
    return Mat3X4T::Identity();
  }
  const SimulatorObject& o=_oss[oid];
  if(o._sysFEM)
    return o._xFEM[id].getT();
  else if(o._sysPBD)
    return TRANSI(o._xPBD[id]._TM,0);
  else {
    ASSERT_MSG(false,"Unknown dynamics type")
  }
}
template <typename T>
void Simulator<T>::advance()
{
  std::vector<SimulatorState> s;
  saveState(s,1);
  loadState(s,2);
  saveState(s,0);
  loadState(s,1);
}
template <typename T>
void Simulator<T>::debug(sizeType nrIter,T scale)
{
  DEFINE_NUMERIC_DELTA_T(T)
  std::vector<SimulatorState> s0,s;
  saveState(s0);
  for(sizeType i=0; i<nrIter; i++) {
    T E,E2;
    Vec dx=Vec::Random(nrDOF());
    Vec G,G2,ctrl=Vec::Random(nrControl());
    Eigen::SparseMatrix<T,0,sizeType> H;
    //evaluate at s0+dx*scale
    updateState(s=s0,dx*scale);
    loadState(s);
    if(!buildSystem(ctrl,&E,&G,&H,false)) {
      i--;
      continue;
    }
    //make sure dx=0 for non-optimizable DOFs
    dx.setRandom();
    for(sizeType i=0,off=0; i<(sizeType)_oss.size(); i++)
      if(_oss[i]._sysFEM) {
        dx.segment(off+_oss[i]._sysFEM->nrDOFOptimizable(),_oss[i]._sysFEM->nrDOF()-_oss[i]._sysFEM->nrDOFOptimizable()).setZero();
        off+=_oss[i]._sysFEM->nrDOF();
      } else if(_oss[i]._sysPBD)
        off+=_oss[i]._body->nrDOF();
      else {
        ASSERT_MSG(false,"Unknown dynamics type")
      }
    //evaluate at perturbed point
    updateState(s,dx*DELTA);
    loadState(s);
    buildSystem(ctrl,&E2,&G2,NULL,false);
    DEBUG_GRADIENT("E",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("G",std::sqrt((H*dx).squaredNorm()),std::sqrt((H*dx-(G2-G)/DELTA).squaredNorm()))
  }
}
template <typename T>
sizeType Simulator<T>::nrControl() const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM)
      ret+=_oss[i]._sysFEM->nrControl();
    else if(_oss[i]._sysPBD)
      ret+=_oss[i]._body->nrDOF();
    else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
  return ret;
}
template <typename T>
sizeType Simulator<T>::nrDOF() const
{
  sizeType ret=0;
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM)
      ret+=_oss[i]._sysFEM->nrDOF();
    else if(_oss[i]._sysPBD)
      ret+=_oss[i]._body->nrDOF();
    else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
  return ret;
}
template <typename T>
const T& Simulator<T>::lastDt() const
{
  return _lastDt;
}
template <typename T>
T& Simulator<T>::lastDt()
{
  return _lastDt;
}
template <typename T>
const T& Simulator<T>::dt() const
{
  return _dt;
}
template <typename T>
T& Simulator<T>::dt()
{
  return _dt;
}
//helper
template <typename T>
std::unordered_set<sizeType> Simulator<T>::attachedChildren(sizeType root) const
{
  std::unordered_set<sizeType> ret;
  ret.insert(root);
  for(const std::pair<sizeType,Mat3X4T>& c:_oss[root]._attachedObjects) {
    std::unordered_set<sizeType> other=attachedChildren(c.first);
    ret.insert(other.begin(),other.end());
  }
  return ret;
}
template <typename T>
void Simulator<T>::saveState(std::vector<SimulatorState>& s,sizeType id) const
{
  s.resize(_oss.size());
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM)
      s[i]._xFEM=_oss[i]._xFEM[id];
    else if(_oss[i]._sysPBD)
      s[i]._xPBD=_oss[i]._xPBD[id];
    else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
}
template <typename T>
void Simulator<T>::loadState(const std::vector<SimulatorState>& s,sizeType id)
{
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM)
      _oss[i]._xFEM[id]=s[i]._xFEM;
    else if(_oss[i]._sysPBD)
      _oss[i]._xPBD[id]=s[i]._xPBD;
    else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
}
template <typename T>
void Simulator<T>::updateState(std::vector<SimulatorState>& s,const Vec& x) const
{
  sizeType ret=0;
  ASSERT(_oss.size()==s.size())
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM) {
      VecCM xM(x.data()+ret,_oss[i]._sysFEM->nrDOF());
      s[i]._xFEM.reset(*(_oss[i]._sysFEM),s[i]._xFEM._x+xM);
      ret+=_oss[i]._sysFEM->nrDOF();
    } else if(_oss[i]._sysPBD) {
      VecCM xM(x.data()+ret,_oss[i]._body->nrDOF());
      s[i]._xPBD.reset(*(_oss[i]._body),s[i]._xPBD._xM+xM);
      ret+=_oss[i]._body->nrDOF();
    } else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
}
template <typename T>
bool Simulator<T>::buildSystem(const Vec& ctrl,T* E,Vec* G,Eigen::SparseMatrix<T,0,sizeType>* H,bool pad) const
{
  sizeType nDOF=nrDOF();
  if(E)
    *E=0;
  if(G)
    G->setZero(nDOF);
  TRIPS trips;
  for(sizeType i=0,off=0,offC=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM) {
      VecCM ctrlI(ctrl.data()+offC,_oss[i]._sysFEM->nrControl());
      if(!_oss[i]._sysFEM->eval(_oss[i]._xFEM[0],_oss[i]._xFEM[1],_oss[i]._xFEM[2],ctrlI,_dt,_lastDt,E,G,H?&trips:NULL,off))
        return false;
      //make sure matrix is SPD
      if(pad)
        for(sizeType j=_oss[i]._sysFEM->nrDOFOptimizable(); j<_oss[i]._sysFEM->nrDOF(); j++)
          trips.push_back(Eigen::Triplet<T,sizeType>(off+j,off+j,1));
      off+=_oss[i]._sysFEM->nrDOF();
      offC+=_oss[i]._sysFEM->nrControl();
    } else {
      VecCM ctrlI(ctrl.data()+offC,_oss[i]._body->nrDOF());
      if(!_oss[i]._sysPBD->eval(_oss[i]._xPBD[0],_oss[i]._xPBD[1],_oss[i]._xPBD[2],ctrlI,_dt,_lastDt,E,G,H?&trips:NULL,off))
        return false;
      off+=_oss[i]._body->nrDOF();
      offC+=_oss[i]._body->nrDOF();
    }
  if(H) {
    H->resize(nDOF,nDOF);
    H->setFromTriplets(trips.begin(),trips.end());
  }
  return true;
}
//instance
PRJ_BEGIN
template class Simulator<double>;
#ifdef ALL_TYPES
template class Simulator<__float128>;
template class Simulator<mpfr::mpreal>;
#endif
PRJ_END
#endif
#endif
#endif
