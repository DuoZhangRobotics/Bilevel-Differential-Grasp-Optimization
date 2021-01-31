#include <Utils/Scalar.h>
#include "DSSQPObjective.h"
#include <CommonFile/Timing.h>
#include <Utils/DebugGradient.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

//KTRObjective
template <typename T>
T DSSQPObjective<T>::infty()
{
  return 1E6;
}
template <typename T>
typename DSSQPObjective<T>::Vec DSSQPObjective<T>::lb() const
{
  return Vec::Constant(inputs(), -infty());
}
template <typename T>
typename DSSQPObjective<T>::Vec DSSQPObjective<T>::ub() const
{
  return Vec::Constant(inputs(), infty());
}
template <typename T>
typename DSSQPObjective<T>::Vec DSSQPObjective<T>::gl() const
{
  return Vec::Constant(values(), 0);
}
template <typename T>
typename DSSQPObjective<T>::Vec DSSQPObjective<T>::gu() const
{
  return Vec::Constant(values(), infty());
}
template <typename T>
typename DSSQPObjective<T>::Vec DSSQPObjective<T>::init() const
{
  return Vec::Zero(inputs(), 0);
}
template <typename T>
int DSSQPObjective<T>::nnzJ()
{
  Vec fvec;
  SMat fjac;
  Vec ret = Vec::Zero(inputs());
  operator()(ret, fvec, &fjac);
  //count nnz
  int off = 0;
  for (sizeType k = 0; k < fjac.outerSize(); ++k)
    for (typename SMat::InnerIterator it(fjac, k); it; ++it, off++);
  return off;
}
//constraint
template <typename T>
int DSSQPObjective<T>::operator()(const Vec& x,Vec& fvec,DMat* fjac) {
  int ret=operator()(x,fvec,fjac?&_tmpFjacs:NULL);
  if(fjac)
    *fjac=_tmpFjacs.toDense();
  return ret;
}
template <typename T>
int DSSQPObjective<T>::operator()(const Vec& x,Vec& fvec,SMat* fjac) {
  if(fjac)  //to remember number of entries
    _tmp.clear();
  int ret=operator()(x,fvec,fjac?&_tmp:NULL);
  if(fjac) {
    fjac->resize(values(),inputs());
    if(!_tmp.getVector().empty())
      fjac->setFromTriplets(_tmp.begin(),_tmp.end());
  }
  return ret;
}
template <typename T>
int DSSQPObjective<T>::operator()(const Vec&,Vec&,STrips*) {
  ASSERT_MSG(false,"Not Implemented: Sparse Constraint Function!")
  return -1;
}
//objective
template <typename T>
T DSSQPObjective<T>::operator()(const Vec& x,Vec* fgrad) {
  return operator()(x,fgrad,(STrips*)NULL);
}
template <typename T>
T DSSQPObjective<T>::operator()(const Vec& x,Vec* fgrad,DMat* fhess) {
  T ret=operator()(x,fgrad,fhess?&_tmpFjacs:NULL);
  if(fhess)
    *fhess=_tmpFjacs.toDense();
  return ret;
}
template <typename T>
T DSSQPObjective<T>::operator()(const Vec& x,Vec* fgrad,SMat* fhess) {
  if(fhess)  //to remember number of entries
    _tmp.clear();
  T ret=operator()(x,fgrad,fhess?&_tmp:NULL);
  if(fhess) {
    fhess->resize(inputs(),inputs());
    if(!_tmp.getVector().empty())
      fhess->setFromTriplets(_tmp.begin(),_tmp.end());
  }
  return ret;
}
template <typename T>
T DSSQPObjective<T>::operator()(const Vec&,Vec*,STrips*) {
  ASSERT_MSG(false,"Not Implemented: Sparse Objective Function!")
  return -1;
}
//problem size
template <typename T>
int DSSQPObjective<T>::inputs() const {
  return 0;
}
template <typename T>
int DSSQPObjective<T>::values() const {
  return 0;
}
//KTRObjectiveComponent
template <typename T>
DSSQPObjectiveComponent<T>::DSSQPObjectiveComponent(DSSQPObjectiveCompound<T>& obj,const std::string& name,bool force):_vars(&(obj.vars())),_offset(0)
{
  _name=name;
  if(force) {
    if(obj.components().find(_name)!=obj.components().end()) {
      sizeType dupId=0;
      do {
        dupId++;
        _name=name+"[Duplicate"+std::to_string(dupId)+"]";
      } while(obj.components().find(_name)!=obj.components().end());
    }
  }
}
template <typename T>
DSSQPObjectiveComponent<T>::~DSSQPObjectiveComponent() {}
template <typename T>
int DSSQPObjectiveComponent<T>::inputs() const
{
  return _inputs;
}
template <typename T>
int DSSQPObjectiveComponent<T>::values() const
{
  ASSERT(_gl.size()==_gu.size())
  return (sizeType)_gl.size();
}
//whether modifying objective function expression is allowed
template <typename T>
void DSSQPObjectiveComponent<T>::setUpdateCache(const Vec&,bool) {}
//constraint
template <typename T>
int DSSQPObjectiveComponent<T>::operator()(const Vec&,Vec&,STrips*)
{
  return 0;
}
//objective
template <typename T>
T DSSQPObjectiveComponent<T>::operator()(const Vec&,Vec*)
{
  //default to pure constraints
  return 0;
}
template <typename T>
bool DSSQPObjectiveComponent<T>::debug(sizeType inputs,sizeType nrTrial,T thres,Vec* x0,bool hess)
{
  sizeType offset=_offset;
  _inputs=inputs;
  _offset=0;
  DEFINE_NUMERIC_DELTA_T(T)
  for(sizeType i=0; i<nrTrial; i++)
    for(sizeType vid=thres!=0?0:-1; vid<sizeType(thres!=0?inputs:0); vid++) {
      T FX,FX2;
      SMat fjac,fhess;
      Vec x=makeValid(Vec::Random(inputs)),grad,grad2;
      if(x0)
        x.segment(0,x0->size())=*x0;
      Vec delta=thres>0?Vec(Vec::Unit(inputs,vid)):Vec(Vec::Random(inputs));
      //objective
      grad.setZero(inputs);
      grad2.setZero(inputs);
      if(hess) {
        setUpdateCache(x,true);
        FX=DSSQPObjective<T>::operator()(x,&grad,&fhess);
        setUpdateCache(x+delta*DELTA,false);
        FX2=DSSQPObjective<T>::operator()(x+delta*DELTA,&grad2,(SMat*)NULL);
        fhess=fhess.toDense().block(0,0,inputs,inputs).sparseView();
      } else {
        setUpdateCache(x,true);
        FX=DSSQPObjective<T>::operator()(x,&grad);
        setUpdateCache(x+delta*DELTA,false);
        FX2=DSSQPObjective<T>::operator()(x+delta*DELTA,&grad2);
      }
      if(thres!=0) {
        std::string entryStr=_name+":grad["+(vid<0?"-":std::to_string(vid))+"]";
        T ref=grad.dot(delta),err=grad.dot(delta)-(FX2-FX)/DELTA;
        if(ref==0)
          continue;
        T thresRel=thres>0?thres:std::max<T>(-thres,-thres*std::abs(ref));
        if(std::abs(err)>thresRel) {
          DEBUG_GRADIENT(entryStr,ref,err)
          return false;
        }
        if(hess) {
          std::string entryStr=_name+":hess["+(vid<0?"-":std::to_string(vid))+"]";
          ref=std::sqrt((hess*delta).squaredNorm()),err=std::sqrt((hess*delta-(grad2-grad)/DELTA).squaredNorm());
          if(ref==0)
            continue;
          T thresRel=thres>0?thres:std::max<T>(-thres,-thres*std::abs(ref));
          if(std::abs(err)>thresRel) {
            DEBUG_GRADIENT(entryStr,ref,err)
            return false;
          }
        }
      } else {
        std::string entryStr=_name+":grad["+(vid<0?"-":std::to_string(vid))+"]";
        T ref=grad.dot(delta),err=grad.dot(delta)-(FX2-FX)/DELTA;
        if(ref==0)
          continue;
        DEBUG_GRADIENT(entryStr,ref,err)
        if(hess) {
          std::string entryStr=_name+":hess["+(vid<0?"-":std::to_string(vid))+"]";
          ref=std::sqrt((fhess*delta).squaredNorm()),err=std::sqrt((fhess*delta-(grad2-grad)/DELTA).squaredNorm());
          if(ref==0)
            continue;
          DEBUG_GRADIENT(entryStr,ref,err)
        }
      }
    }
  for(sizeType i=0; i<nrTrial; i++)
    for(sizeType vid=thres!=0?0:-1; vid<sizeType(thres!=0?inputs:0); vid++) {
      SMat fjac;
      Vec x=makeValid(Vec::Random(inputs)),grad,grad2;
      if(x0)
        x.segment(0,x0->size())=*x0;
      Vec delta=thres>0?Vec(Vec::Unit(inputs,vid)):Vec(Vec::Random(inputs));
      //constraint
      Vec fvec=Vec::Zero(values());
      Vec fvec2=Vec::Zero(values());
      fjac.resize(values(),inputs);
      setUpdateCache(x,true);
      DSSQPObjective<T>::operator()(x,fvec,(SMat*)&fjac);
      setUpdateCache(x+delta*DELTA,false);
      DSSQPObjective<T>::operator()(x+delta*DELTA,fvec2,(SMat*)NULL);
      if(thres!=0) {
        for(sizeType r=0; r<fvec.size(); r++) {
          std::string entryStr=_name+":fJac("+std::to_string(r)+","+(vid<0?"-":std::to_string(vid))+")";
          T ref=(fjac*delta)[r],err=(fjac*delta-(fvec2-fvec)/DELTA)[r];
          if(ref==0 && fvec2==fvec)
            continue;
          T thresRel=thres>0?thres:std::max<T>(-thres,-thres*std::abs(ref));
          if(std::abs(err)>thresRel) {
            DEBUG_GRADIENT(entryStr,ref,err)
            return false;
          }
        }
      } else {
        std::string entryStr=_name+":fJac(-,"+(vid<0?"-":std::to_string(vid))+")";
        T ref=std::sqrt((fjac*delta).squaredNorm()),err=std::sqrt((fjac*delta-(fvec2-fvec)/DELTA).squaredNorm());
        if(ref==0 && fvec2==fvec)
          continue;
        DEBUG_GRADIENT(entryStr,ref,err)
      }
    }
  _offset=offset;
  return true;
}
template <typename T>
bool DSSQPObjectiveComponent<T>::debug(sizeType nrTrial,T thres,Vec* x0,bool hess)
{
  sizeType inputs=0;
  for(const std::pair<std::string,DSSQPVariable<T>>& v:*_vars)
    inputs=std::max(inputs,v.second._id+1);
  return debug(inputs,nrTrial,thres,x0,hess);
}
template <typename T>
typename DSSQPObjectiveComponent<T>::Vec DSSQPObjectiveComponent<T>::makeValid(const Vec& x) const
{
  return x;
}
template <typename T>
bool DSSQPObjectiveComponent<T>::debugGradientConditional(const std::string& entryStr,T ref,T err,T thres) const
{
  if(ref==0 && err==0)
    return true;
  std::string entryStrName=_name+":"+entryStr;
  DEFINE_NUMERIC_DELTA_T(T)
  if(thres!=0) {
    T thresRel=thres>0?thres:std::max<T>(-thres,-thres*std::abs(ref));
    if(std::abs(err)>thresRel) {
      DEBUG_GRADIENT(entryStrName,ref,err)
      return false;
    }
  } else {
    DEBUG_GRADIENT(entryStrName,ref,err)
  }
  return true;
}
template <typename T>
bool DSSQPObjectiveComponent<T>::debugGradientConditionalVID(const std::string& entryStr,sizeType vid,T ref,T err,T thres) const
{
  if(vid==-1)
    return debugGradientConditional(entryStr,ref,err,thres);
  else return debugGradientConditional(entryStr+"("+std::to_string(vid)+")",ref,err,thres);
}
template <typename T>
void DSSQPObjectiveComponent<T>::setOffset(sizeType offset)
{
  _offset=offset;
}
//KTRObjectiveCompound
template <typename T>
typename DSSQPObjectiveCompound<T>::Vec DSSQPObjectiveCompound<T>::lb() const
{
  Vec ret=Vec::Zero(_vars.size());
  for(const std::pair<std::string,DSSQPVariable<T>>& v:_vars)
    ret[v.second._id]=v.second._l;
  return ret;
}
template <typename T>
typename DSSQPObjectiveCompound<T>::Vec DSSQPObjectiveCompound<T>::ub() const
{
  Vec ret=Vec::Zero(_vars.size());
  for(const std::pair<std::string,DSSQPVariable<T>>& v:_vars)
    ret[v.second._id]=v.second._u;
  return ret;
}
template <typename T>
typename DSSQPObjectiveCompound<T>::Vec DSSQPObjectiveCompound<T>::gl() const
{
  Vec ret=Vec::Zero(values());
  for(const std::pair<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>& v:_components)
    ret.segment(v.second->_offset,v.second->values())=Eigen::Map<const Vec>(&(v.second->_gl[0]),v.second->_gl.size());
  return ret;
}
template <typename T>
typename DSSQPObjectiveCompound<T>::Vec DSSQPObjectiveCompound<T>::gu() const
{
  Vec ret=Vec::Zero(values());
  for(const std::pair<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>& v:_components)
    ret.segment(v.second->_offset,v.second->values())=Eigen::Map<const Vec>(&(v.second->_gu[0]),v.second->_gu.size());
  return ret;
}
template <typename T>
typename DSSQPObjectiveCompound<T>::Vec DSSQPObjectiveCompound<T>::init() const
{
  Vec ret=Vec::Zero(_vars.size());
  for(const std::pair<std::string,DSSQPVariable<T>>& v:_vars)
    ret[v.second._id]=v.second._init;
  return ret;
}
template <typename T>
int DSSQPObjectiveCompound<T>::inputs() const
{
  return _vars.size();
}
template <typename T>
int DSSQPObjectiveCompound<T>::values() const
{
  int ret=0;
  for(const std::pair<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>& v:_components)
    ret+=v.second->values();
  return ret;
}
//constraint
template <typename T>
int DSSQPObjectiveCompound<T>::operator()(const Vec& x,Vec& fvec,STrips* fjac)
{
  fvec.setZero(values());
  if(fjac)
    fjac->clear();
  for(const std::pair<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>& v:_components)
    v.second->operator()(x,fvec,fjac);
  return 0;
}
template <typename T>
const std::vector<Coli,Eigen::aligned_allocator<Coli>>& DSSQPObjectiveCompound<T>::getQCones() const
{
  return _QCones;
}
template <typename T>
void DSSQPObjectiveCompound<T>::addQCone(const Coli& vss)
{
  for(sizeType i=0; i<vss.size(); i++) {
    ASSERT_MSGV(vss[i]>=0 && vss[i]<(sizeType)vars().size(),"Invalid variable id: 0<=%d<%d",vss[i],vars().size())
  }
  _QCones.push_back(vss);
}
//objective
template <typename T>
T DSSQPObjectiveCompound<T>::operator()(const Vec& x,Vec* fgrad)
{
  T FX=0;
  if(fgrad)
    fgrad->setZero(inputs());
  for(const std::pair<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>& v:_components)
    FX+=v.second->operator()(x,fgrad);
  return FX;
}
template <typename T>
const DSSQPVariable<T>& DSSQPObjectiveCompound<T>::addVar(const std::string& name,T l,T u,VARIABLE_OP op)
{
  typename VARMAP::iterator it=_vars.find(name);
  if(op==MUST_EXIST) {
    ASSERT_MSGV(it!=_vars.end(),"Variable: %s does not exist!",name.c_str())
  } else {
    if(op==MUST_NEW) {
      ASSERT_MSGV(it==_vars.end(),"Variable: %s already exists!",name.c_str())
    }
    if(it==_vars.end()) {
      _vars[name]=DSSQPVariable<T>();
      _vars[name]._id=(sizeType)_vars.size()-1;
      _varsInv[(sizeType)_vars.size()-1]=name;
      it=_vars.find(name);
    }
  }
  it->second._l=l;
  it->second._u=u;
  return it->second;
}
template <typename T>
const DSSQPVariable<T>& DSSQPObjectiveCompound<T>::addVar(const std::string& name,VARIABLE_OP op)
{
  typename VARMAP::iterator it=_vars.find(name);
  if(op==MUST_EXIST) {
    ASSERT_MSGV(it!=_vars.end(),"Variable: %s does not exist!",name.c_str())
  } else {
    if(op==MUST_NEW) {
      ASSERT_MSGV(it==_vars.end(),"Variable: %s already exists!",name.c_str())
    }
    if(it==_vars.end()) {
      _vars[name]=DSSQPVariable<T>();
      _vars[name]._id=(sizeType)_vars.size()-1;
      _varsInv[(sizeType)_vars.size()-1]=name;
      it=_vars.find(name);
    }
  }
  return it->second;
}
template <typename T>
const DSSQPVariable<T>& DSSQPObjectiveCompound<T>::addVar(sizeType id) const
{
  ASSERT(_varsInv.find(id)!=_varsInv.end())
  return _vars.find(_varsInv.find(id)->second)->second;
}
template <typename T>
DSSQPVariable<T>& DSSQPObjectiveCompound<T>::addVar(sizeType id)
{
  ASSERT(_varsInv.find(id)!=_varsInv.end())
  return _vars.find(_varsInv.find(id)->second)->second;
}
template <typename T>
const typename DSSQPObjectiveCompound<T>::VARMAP& DSSQPObjectiveCompound<T>::vars() const
{
  return _vars;
}
template <typename T>
void DSSQPObjectiveCompound<T>::setVarInit(const std::string& name,T init)
{
  typename VARMAP::iterator it=_vars.find(name);
  ASSERT_MSGV(it!=_vars.end(),"Variable: %s does not exist!",name.c_str())
  it->second._init=init;
}
template <typename T>
void DSSQPObjectiveCompound<T>::setVarInit(sizeType id,T init)
{
  ASSERT(_varsInv.find(id)!=_varsInv.end())
  _vars.find(_varsInv.find(id)->second)->second._init=init;
}
template <typename T>
void DSSQPObjectiveCompound<T>::checkViolation(const Vec* at)
{
  Vec fvec,grad;
  Vec l=gl();
  Vec u=gu();
  STrips fjac;
  sizeType index;
  Vec x=at?*at:init();
  operator()(x,fvec,&fjac);
  for(const std::pair<std::string,std::shared_ptr<DSSQPObjectiveComponent<T>>>& v:_components) {
    T maxVioBelow,maxVioAbove;
    Vec grad=Vec::Zero(inputs()),vioMin,vioMax;
    T obj=v.second->operator()(x,&grad);
    if(v.second->values()==0) {
      maxVioBelow=maxVioAbove=0;
    } else {
      vioMin=(fvec-l).segment(v.second->_offset,v.second->values());
      maxVioBelow=std::min<T>(vioMin.minCoeff(&index),0);
      vioMax=(fvec-u).segment(v.second->_offset,v.second->values());
      maxVioAbove=std::max<T>(vioMax.maxCoeff(&index),0);
    }
    std::ostringstream oss;
    oss << "Component " << v.second->_name << ": violation=(below: " << maxVioBelow << " above: " << maxVioAbove << ") objective=" << obj << "!";
    INFO(oss.str().c_str())
  }
}
template <typename T>
const typename DSSQPObjectiveCompound<T>::CONSMAP& DSSQPObjectiveCompound<T>::components() const
{
  return _components;
}
template <typename T>
void DSSQPObjectiveCompound<T>::addComponent(std::shared_ptr<DSSQPObjectiveComponent<T>> c)
{
  ASSERT_MSGV(_components.find(c->_name)==_components.end(),"Constraint %s already exists!",c->_name.c_str())
  c->setOffset(values());
  _components[c->_name]=c;
}
template <typename T>
bool DSSQPObjectiveCompound<T>::debug(const std::string& str,sizeType nrTrial,T thres)
{
  ASSERT_MSGV(_components.find(str)!=_components.end(),"Cannot find component: %s!",str.c_str())
  std::shared_ptr<DSSQPObjectiveComponent<T>> comp=_components.find(str)->second;
  sizeType offset=comp->_offset;
  comp->setOffset(0);
  bool ret=comp->debug(inputs(),nrTrial,thres);
  comp->setOffset(offset);
  return ret;
}
//instance
PRJ_BEGIN
#define INSTANTIATE(T)    \
template struct DSSQPObjective<T>;  \
template struct DSSQPObjectiveComponent<T>; \
template struct DSSQPObjectiveCompound<T>;
INSTANTIATE(double)
#ifdef ALL_TYPES
INSTANTIATE(__float128)
INSTANTIATE(mpfr::mpreal)
#endif
PRJ_END
