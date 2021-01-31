#ifdef OPTIMIZER_SUPPORT
#ifdef ENVIRONMENT_SUPPORT
#ifdef DEFORMABLE_SUPPORT
#include "Simulator.h"
#include <Deformable/FEMSystem.h>
#include <Articulated/JointFunc.h>
#include <Articulated/PBDSimulator.h>
#include <Optimizer/QCQPSolverQPOASES.h>
//#include <Optimizer/QCQPSolverMosek.h>
//#include <Optimizer/QCQPSolverGurobi.h>
#include "CollisionConstraint.h"
#include <Utils/DebugGradient.h>

USE_PRJ_NAMESPACE

template <typename T>
void callbackSimulator(const std::string& callbackPrefix,sizeType it,T E,T alpha,T reg,T dMax,const std::string& info="")
{
  if(callbackPrefix.empty())
    return;
  if(it>=0)
    std::cout << "iter=" << it;
  else std::cout << "frm=" << callbackPrefix;
  if(it>=0)
    std::cout << " E=" << E;
  std::cout << " alpha=" << alpha << " reg=" << reg;
  if(it>=0)
    std::cout << " dMax=" << dMax;
  if(!info.empty())
    std::cout << " " << info;
  std::cout << std::endl;
}
//SimulatorObject
template <typename T>
sizeType SimulatorObject<T>::nrControl() const
{
  return _sysFEM?_sysFEM->nrControl():_body->nrDOF();
}
template <typename T>
sizeType SimulatorObject<T>::nrDOF() const
{
  return _sysFEM?_sysFEM->nrDOF():_body->nrDOF();
}
//Simulator
template <typename T>
Simulator<T>::Simulator():_dt(0.01f),_lastDt(0.01f),_dtMin(0.001f),_dtMax(0.01f),_dist(0.1f),_dThres(1e-4f),_rThres(1),_aThres(1e-20f),_itMax(1000),_fDir(6) {}
template <typename T>
typename Simulator<T>::Vec Simulator<T>::getState(sizeType id) const
{
  Vec ret=Vec::Zero(nrDOF());
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM)
      ret.segment(_oss[i]._offDOF,_oss[i].nrDOF())=_oss[i]._xFEM[id]._x.segment(0,_oss[i].nrDOF());
    else if(_oss[i]._sysPBD)
      ret.segment(_oss[i]._offDOF,_oss[i].nrDOF())=_oss[i]._xPBD[id]._xM.segment(0,_oss[i].nrDOF());
    else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
  return ret;
}
template <typename T>
void Simulator<T>::setState(const Vec& x,sizeType id)
{
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM)
      _oss[i]._xFEM[id].reset(*(_oss[i]._sysFEM),x.segment(_oss[i]._offDOF,_oss[i].nrDOF()));
    else if(_oss[i]._sysPBD)
      _oss[i]._xPBD[id].reset(*(_oss[i]._body),x.segment(_oss[i]._offDOF,_oss[i].nrDOF()));
    else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
}
template <typename T>
sizeType Simulator<T>::addFEMObject(std::shared_ptr<FEMSystem<T>> sys,const Vec& xInit)
{
  SimulatorObject<T> obj;
  obj._sysFEM=sys;
  obj._xFEM[0].reset(*sys,mapV(xInit));
  obj._xFEM[2]=obj._xFEM[1]=obj._xFEM[0];
  obj._offControl=_oss.empty()?0:_oss.back()._offControl+_oss.back().nrControl();
  obj._offDOF=_oss.empty()?0:_oss.back()._offDOF+_oss.back().nrDOF();
  _css.push_back(std::shared_ptr<CollisionNode<T>>(new CollisionNode<T>(obj,_oss.size(),-1)));
  _sap.add(_css.back()->getBB());
  _oss.push_back(obj);
  _lb=concat<Vec>(_lb,Vec::Constant(obj.nrDOF(),-DSSQPObjective<T>::infty()));
  _ub=concat<Vec>(_ub,Vec::Constant(obj.nrDOF(), DSSQPObjective<T>::infty()));
  return (sizeType)_oss.size()-1;
}
template <typename T>
sizeType Simulator<T>::addPBDObject(std::shared_ptr<ArticulatedBody> body,const Vec3T& g,const Vec& xInit)
{
  SimulatorObject<T> obj;
  obj._body=body;
  Options ops;
  obj._sysPBD.reset(new PBDSimulator<T>(*body,ops,g,PBTO));
  obj._xPBD[0].reset(*body,mapV(xInit));
  obj._xPBD[2]=obj._xPBD[1]=obj._xPBD[0];
  obj._offControl=_oss.empty()?0:_oss.back()._offControl+_oss.back().nrControl();
  obj._offDOF=_oss.empty()?0:_oss.back()._offDOF+_oss.back().nrDOF();
  for(sizeType i=0; i<obj._body->nrJ(); i++)
    if(obj._body->joint(i)._M>0) {
      _css.push_back(std::shared_ptr<CollisionNode<T>>(new CollisionNode<T>(obj,_oss.size(),i)));
      _sap.add(_css.back()->getBB());
    }
  _oss.push_back(obj);
  _lb=concat<Vec>(_lb,body->lowerLimit(DSSQPObjective<scalarD>::infty()).template cast<T>());
  _ub=concat<Vec>(_ub,body->upperLimit(DSSQPObjective<scalarD>::infty()).template cast<T>());
  return (sizeType)_oss.size()-1;
}
template <typename T>
void Simulator<T>::writeConstraintsVTK(const std::string& path) const
{
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  for(sizeType i=0; i<(sizeType)_constraints.size(); i++) {
    vss.push_back(_constraints[i]._pL.template cast<scalar>());
    vss.push_back(_constraints[i]._pR.template cast<scalar>());
  }
  VTKWriter<scalar> os("constraints",path,true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(VTKWriter<scalar>::IteratorIndex<Vec3i>(0,2,0),
                 VTKWriter<scalar>::IteratorIndex<Vec3i>(vss.size()/2,2,0),
                 VTKWriter<scalar>::LINE);
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
bool Simulator<T>::step(const Vec& ctrl,const std::string& callbackPrefix)
{
  if(_dtMin==_dtMax)
    return stepInner(ctrl,callbackPrefix);
  else return stepInnerAdaptive(_dtMax,ctrl,callbackPrefix);
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
  SimulatorObject<T>& o=_oss[oid];
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
  const SimulatorObject<T>& o=_oss[oid];
  if(o._sysFEM)
    return o._xFEM[id].getT();
  else if(o._sysPBD)
    return TRANSI(o._xPBD[id]._TM,0);
  else {
    ASSERT_MSG(false,"Unknown dynamics type")
  }
}
template <typename T>
void Simulator<T>::debugCollisionConstraint(sizeType N,sizeType nrIter,T scale)
{
  _constraints.clear();
  while((sizeType)_constraints.size()<N) {
    CollisionConstraint<T> c;
    c._L=_css[RandEngine::randI(0,_css.size()-1)];
    c._R=_css[RandEngine::randI(0,_css.size()-1)];
    if(c._L==c._R)
      continue;
    c._pL.setRandom();
    c._pR.setRandom();
    c._bary.setRandom();
    if(RandEngine::randR01()>0.5f) {
      c._isP2T=true;
      c._vid[0]=RandEngine::randI(0,c._L->nrV()-1);
      c._vid[1]=RandEngine::randI(0,c._R->nrV()-1);
      c._vid[2]=RandEngine::randI(0,c._R->nrV()-1);
      c._vid[3]=RandEngine::randI(0,c._R->nrV()-1);
    } else {
      c._isP2T=false;
      c._vid[0]=RandEngine::randI(0,c._L->nrV()-1);
      c._vid[1]=RandEngine::randI(0,c._L->nrV()-1);
      c._vid[2]=RandEngine::randI(0,c._R->nrV()-1);
      c._vid[3]=RandEngine::randI(0,c._R->nrV()-1);
    }
    _constraints.push_back(c);
  }

  DEFINE_NUMERIC_DELTA_T(T)
  std::vector<SimulatorState<T>> s0,s;
  saveState(s0);
  Vec c,c2;
  SMat J;
  for(sizeType i=0; i<nrIter; i++) {
    Vec dx=Vec::Random(nrDOF());
    //evaluate at s0+dx*scale
    updateState(s=s0,dx*scale);
    loadState(s);
    c=getCollisionConstraintJacobian(&J);
    //evaluate at perturbed point
    updateState(s,dx*DELTA);
    loadState(s);
    c2=getCollisionConstraintJacobian(NULL);
    DEBUG_GRADIENT("J",std::sqrt((J*dx).squaredNorm()),std::sqrt((J*dx-(c2-c)/DELTA).squaredNorm()))
  }
}
template <typename T>
void Simulator<T>::debug(sizeType nrIter,T scale)
{
  DEFINE_NUMERIC_DELTA_T(T)
  std::vector<SimulatorState<T>> s0,s;
  saveState(s0);
  for(sizeType i=0; i<nrIter; i++) {
    T E,E2;
    Vec dx=Vec::Random(nrDOF());
    Vec G,G2,ctrl=Vec::Random(nrControl());
    SMat H;
    //evaluate at s0+dx*scale
    updateState(s=s0,dx*scale);
    loadState(s);
    if(!buildSystem(ctrl,&E,&G,&H,NULL,false)) {
      i--;
      continue;
    }
    //make sure dx=0 for non-optimizable DOFs
    dx.setRandom();
    for(sizeType i=0; i<(sizeType)_oss.size(); i++)
      if(_oss[i]._sysFEM) {
        sizeType nrDOF=_oss[i]._sysFEM->nrDOF(),nrDOFO=_oss[i]._sysFEM->nrDOFOptimizable();
        dx.segment(_oss[i]._offDOF+nrDOFO,nrDOF-nrDOFO).setZero();
      }
    //evaluate at perturbed point
    updateState(s,dx*DELTA);
    loadState(s);
    buildSystem(ctrl,&E2,&G2,NULL,NULL,false);
    DEBUG_GRADIENT("E",G.dot(dx),G.dot(dx)-(E2-E)/DELTA)
    DEBUG_GRADIENT("G",std::sqrt((H*dx).squaredNorm()),std::sqrt((H*dx-(G2-G)/DELTA).squaredNorm()))
  }
}
//parameter
template <typename T>
sizeType Simulator<T>::nrControl() const
{
  return _oss.empty()?0:_oss.back()._offControl+_oss.back().nrControl();
}
template <typename T>
sizeType Simulator<T>::nrDOF() const
{
  return _oss.empty()?0:_oss.back()._offDOF+_oss.back().nrDOF();
}
template <typename T>
const T& Simulator<T>::dtMin() const
{
  return _dtMin;
}
template <typename T>
T& Simulator<T>::dtMin()
{
  return _dtMin;
}
template <typename T>
const T& Simulator<T>::dtMax() const
{
  return _dtMax;
}
template <typename T>
T& Simulator<T>::dtMax()
{
  return _dtMax;
}
template <typename T>
const T& Simulator<T>::dist() const
{
  return _dist;
}
template <typename T>
T& Simulator<T>::dist()
{
  return _dist;
}
template <typename T>
const T& Simulator<T>::dThres() const
{
  return _dThres;
}
template <typename T>
T& Simulator<T>::dThres()
{
  return _dThres;
}
template <typename T>
const T& Simulator<T>::rThres() const
{
  return _rThres;
}
template <typename T>
T& Simulator<T>::rThres()
{
  return _rThres;
}
template <typename T>
const T& Simulator<T>::aThres() const
{
  return _aThres;
}
template <typename T>
T& Simulator<T>::aThres()
{
  return _aThres;
}
template <typename T>
const sizeType& Simulator<T>::itMax() const
{
  return _itMax;
}
template <typename T>
sizeType& Simulator<T>::itMax()
{
  return _itMax;
}
//helper
template <typename T>
void Simulator<T>::advance()
{
  std::vector<SimulatorState<T>> s;
  saveState(s,1);
  loadState(s,2);
  saveState(s,0);
  loadState(s,1);
  _lastDt=_dt;
}
template <typename T>
bool Simulator<T>::stepInnerAdaptive(T dt,const Vec& ctrl,const std::string& callbackPrefix)
{
  if(dt<_dtMin)
    return false;

  //try larger step
  _dt=dt;
  std::string prefix=callbackPrefix;
  if(!prefix.empty())
    prefix+="[lastDt="+std::to_string(_lastDt)+",dt="+std::to_string(dt)+"]";
  if(stepInner(ctrl,prefix))
    return true;

  //try subdivide step
  _dt=dt/2;
  if(!stepInnerAdaptive(dt/2,ctrl,prefix))
    return false;
  if(!stepInnerAdaptive(dt/2,ctrl,prefix))
    return false;
  return true;
}
template <typename T>
bool Simulator<T>::stepInner(const Vec& ctrl,const std::string& callbackPrefix)
{
  Vec G,C,d,lbC;
  bool valid,succ=true;
  T regMin=1e-5f,reg=regMin;
  T E=0,ENew,alpha=1.0f,dMax=0;
  std::vector<SimulatorState<T>> state,stateNew;
  SMat H,HD,J,diagonalPerturb;
  QCQPSolverQPOASES<T> sol;
  sizeType it;

  //initialize: copy state[1] to state[0]
  saveState(state,1);
  loadState(state,0);

  //detect collision constraint
  Vec x=getState(),lb=_lb-x,ub=_ub-x;
  getCollisionConstraint();
  C=getCollisionConstraintJacobian(&J);
  callbackSimulator(callbackPrefix,-1,E,alpha,reg,dMax,"#collision-constraint="+std::to_string(_constraints.size()));

  //main loop
  for(it=0; alpha>_aThres && it<_itMax; it++) {
    valid=buildSystem(ctrl,&E,&G,&H,diagonalPerturb.size()==0?&diagonalPerturb:NULL);
    if(!valid) {
      callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax,"failed(invalid system)");
      succ=false;
      break;
    }
    //solve
    lbC=Vec::Constant(C.size(),_dist/2)-C;
    while(reg<_rThres) {
      HD=H+diagonalPerturb*reg;
      typename QCQPSolver<T>::QCQP_RETURN_CODE ret=sol.solveQP(d=Vec::Zero(nrDOF()),HD,G,J.rows()>0?&J:NULL,&lb,&ub,J.rows()>0?&lbC:NULL,NULL,std::vector<Coli,Eigen::aligned_allocator<Coli>>(),false);
      if(ret==QCQPSolver<T>::SOLVED) {
        reg=std::max(reg*0.5f,regMin);
        break;
      } else if(ret==QCQPSolver<T>::NOT_POSITIVE_DEFINITE)
        reg*=10;
      else if(ret==QCQPSolver<T>::INFEASIBLE) {
        callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax,"failed(QP infeasible)");
        succ=false;
        break;
      } else {
        reg*=10;
        //callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax,"failed(QP unknown error)");
        //succ=false;
        //break;
      }
    }
    if(!succ)
      break;
    if(reg>=_rThres) {
      callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax,"failed(QP too large regularization)");
      succ=false;
      break;
    }
    //termination
    dMax=d.unaryExpr([&](const T& in) {
      return std::abs(in);
    }).maxCoeff();
    if(dMax<_dThres) {
      callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax,"succeeded");
      break;
    }
    //line search
    saveState(state);
    while(alpha>_aThres) {
      updateState(stateNew=state,d*alpha);
      loadState(stateNew);
      if(it==0)   //since constraint has changed, always take the first step
        break;
      valid=buildSystem(ctrl,&ENew,NULL,NULL,NULL);
      //update line search
      if(valid && std::isfinite(ENew) && ENew<E) {
        alpha=std::min<T>(alpha*1.5f,1.f);
        callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax);
        break;
      } else alpha*=0.5f;
    }
    if(alpha<_aThres) {
      callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax,"failed(small alpha)");
      //for debug
      /*std::cout << "d: " << d.unaryExpr([&](T in) {
        return (scalarD)std::to_double(in);
      }).transpose() << std::endl;
      std::cout << "lbC: " << lbC.unaryExpr([&](T in) {
        return (scalarD)std::to_double(in);
      }).transpose() << std::endl;
      std::cout << "lb: " << lb.unaryExpr([&](T in) {
        return (scalarD)std::to_double(in);
      }).transpose() << std::endl;
      std::cout << "ub: " << ub.unaryExpr([&](T in) {
        return (scalarD)std::to_double(in);
      }).transpose() << std::endl;
      typename QCQPSolver<T>::QCQP_RETURN_CODE ret=sol.solveQP(d=Vec::Zero(nrDOF()),HD,G,J.rows()>0?&J:NULL,&lb,&ub,J.rows()>0?&lbC:NULL,NULL,std::vector<Coli,Eigen::aligned_allocator<Coli>>(),false);
      std::cout << std::to_double(d.dot(G+HD*d/2)) << std::endl;
      debug(10,0);*/
      succ=false;
      break;
    }
    //update constraint
    C+=J*d*alpha;
    lb-=d*alpha;
    ub-=d*alpha;
  }

  //make sure no collision after solve
  if(succ && checkCollision()) {
    callbackSimulator(callbackPrefix,it,E,alpha,reg,dMax,"failed(collision after solve)");
    succ=false;
  }
  if(succ)
    advance();
  return succ;
}
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
void Simulator<T>::saveState(std::vector<SimulatorState<T>>& s,sizeType id) const
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
void Simulator<T>::loadState(const std::vector<SimulatorState<T>>& s,sizeType id)
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
void Simulator<T>::updateState(std::vector<SimulatorState<T>>& s,const Vec& x) const
{
  ASSERT(_oss.size()==s.size())
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM) {
      VecCM xM(x.data()+_oss[i]._offDOF,_oss[i]._sysFEM->nrDOF());
      s[i]._xFEM.reset(*(_oss[i]._sysFEM),s[i]._xFEM._x+xM);
    } else if(_oss[i]._sysPBD) {
      VecCM xM(x.data()+_oss[i]._offDOF,_oss[i]._body->nrDOF());
      s[i]._xPBD.reset(*(_oss[i]._body),s[i]._xPBD._xM+xM);
    } else {
      ASSERT_MSG(false,"Unknown dynamics type")
    }
}
template <typename T>
bool Simulator<T>::buildSystem(const Vec& ctrl,T* E,Vec* G,SMat* H,SMat* diagonalPerturb,bool pad) const
{
  sizeType nDOF=nrDOF();
  if(E)
    *E=0;
  if(G)
    G->setZero(nDOF);
  TRIPS trips,tripsD;
  for(sizeType i=0; i<(sizeType)_oss.size(); i++)
    if(_oss[i]._sysFEM) {
      VecCM ctrlI(ctrl.data()+_oss[i]._offControl,_oss[i]._sysFEM->nrControl());
      if(!_oss[i]._sysFEM->eval(_oss[i]._xFEM[0],_oss[i]._xFEM[1],_oss[i]._xFEM[2],ctrlI,_dt,_lastDt,E,G,H?&trips:NULL,_oss[i]._offDOF))
        return false;
      //make sure matrix is SPD
      if(pad)
        for(sizeType j=_oss[i]._sysFEM->nrDOFOptimizable(); j<_oss[i]._sysFEM->nrDOF(); j++)
          trips.push_back(Eigen::Triplet<T,sizeType>(_oss[i]._offDOF+j,_oss[i]._offDOF+j,1));
    } else {
      VecCM ctrlI(ctrl.data()+_oss[i]._offControl,_oss[i]._body->nrDOF());
      if(!_oss[i]._sysPBD->eval(_oss[i]._xPBD[0],_oss[i]._xPBD[1],_oss[i]._xPBD[2],ctrlI,_dt,_lastDt,E,G,H?&trips:NULL,_oss[i]._offDOF))
        return false;
    }
  if(H) {
    H->resize(nDOF,nDOF);
    H->setFromTriplets(trips.begin(),trips.end());
  }
  if(diagonalPerturb) {
    diagonalPerturb->resize(nDOF,nDOF);
    sizeType oid=0;
    std::unordered_map<sizeType,T> oDiag;
    for(sizeType k=0; k<H->outerSize(); ++k)
      for(typename SMat::InnerIterator it(*H,k); it; ++it)
        if(it.row()==it.col()) {
          while(it.row()>_oss[oid]._offDOF+_oss[oid].nrDOF())
            oid++;
          oDiag[oid]=std::max<T>(oDiag[oid],std::abs(it.value()));
        }
    for(oid=0; oid<(sizeType)_oss.size(); oid++)
      for(sizeType d=0; d<_oss[oid].nrDOF(); d++)
        tripsD.push_back(STrip(d+_oss[oid]._offDOF,d+_oss[oid]._offDOF,oDiag[oid]));
    diagonalPerturb->setFromTriplets(tripsD.begin(),tripsD.end());
  }
  return true;
}
template <typename T>
typename Simulator<T>::Vec Simulator<T>::getCollisionConstraintJacobian(SMat* J)
{
  STrips trips;
  Vec ret=Vec::Zero(_constraints.size());
  for(sizeType i=0; i<(sizeType)_constraints.size(); i++) {
    CollisionConstraint<T>& c=_constraints[i];
    c.buildFrame(_fDir);
    ret[i]=c.assembleSystem(_oss[c._L->objId()],_oss[c._R->objId()],i,J?&trips:NULL);
  }
  if(J) {
    J->resize(_constraints.size(),nrDOF());
    J->setFromTriplets(trips.begin(),trips.end());
  }
  return ret;
}
template <typename T>
void Simulator<T>::getCollisionConstraint(sizeType id)
{
  //update collision
  for(sizeType i=0; i<(sizeType)_css.size(); i++) {
    sizeType jid=_css[i]->objId();
    _css[i]->updateBVH(_oss[jid],id);
    _sap.update(_css[i]->getBB().enlarge(std::to_double(_dist)),i);
  }
  //get constraints
  _constraints.clear();
  _sap.intersect([&](sizeType i,sizeType j) {
    if(filterCollision(i,j))
      return;
    CollisionNode<T>::addConstraint(_constraints,_css[i],_css[j],_dist);
  });
}
template <typename T>
bool Simulator<T>::checkCollision(sizeType id)
{
  //update collision
  for(sizeType i=0; i<(sizeType)_css.size(); i++) {
    sizeType jid=_css[i]->objId();
    _css[i]->updateBVH(_oss[jid],id);
    _sap.update(_css[i]->getBB(),i);
  }
  bool hasColl=false;
  _sap.intersect([&](sizeType i,sizeType j) {
    if(hasColl)
      return;
    if(filterCollision(i,j))
      return;
    else if(CollisionNode<T>::checkCollision(_css[i],_css[j]))
      hasColl=true;
  });
  return hasColl;
}
template <typename T>
bool Simulator<T>::filterCollision(sizeType i,sizeType j) const
{
  if(_css[i]->objId()==_css[j]->objId() && _css[i]->jointId()>=0 && _css[j]->jointId()>=0) {
    const Joint& JI=_oss[_css[i]->objId()]._body->joint(_css[i]->jointId());
    if(JI._parent==_css[j]->jointId())
      return true;
    const Joint& JJ=_oss[_css[j]->objId()]._body->joint(_css[j]->jointId());
    if(JJ._parent==_css[i]->jointId())
      return true;
  }
  return false;
}
//instance
PRJ_BEGIN
template class SimulatorObject<double>;
template class SimulatorState<double>;
template class Simulator<double>;

#ifdef ALL_TYPES
template class SimulatorObject<__float128>;
template class SimulatorState<__float128>;
template class Simulator<__float128>;

template class SimulatorObject<mpfr::mpreal>;
template class SimulatorState<mpfr::mpreal>;
template class Simulator<mpfr::mpreal>;
#endif
PRJ_END
#endif
#endif
#endif
