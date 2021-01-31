#include "ContactIndexConstraint.h"
#include "GraspPlanner.h"
#include <Environment/Environment.h>
#include <Articulated/PBDArticulatedGradientInfo.h>

USE_PRJ_NAMESPACE

template <typename T>
bool ContactIndexConstraint<T>::ContactIndexMap::operator<(const ContactIndexMap& other) const
{
  if(_oid<other._oid)
    return true;
  else if(_oid>other._oid)
    return false;

  if(_oidOther<other._oidOther)
    return true;
  else if(_oidOther>other._oidOther)
    return false;

  if(_pid<other._pid)
    return true;
  else if(_pid>other._pid)
    return false;

  return false;
}
template <typename T>
bool ContactIndexConstraint<T>::ContactIndexMap::operator!=(const ContactIndexMap& other) const {
  return *this<other || other<*this;
}
template <typename T>
bool ContactIndexConstraint<T>::ContactIndexMap::operator==(const ContactIndexMap& other) const {
  return !(*this<other) && !(other<*this);
}
template <typename T>
ContactIndexConstraint<T>::ContactIndexConstraint(DSSQPObjectiveCompound<T>& obj,const PBDArticulatedGradientInfo<T>& info,const GraspPlanner<T>& planner,const PointCloudObject<T>& object,const Vec3T& g,const T& tau)
  :ArticulatedObjective<T>(obj,"ContactIndexConstraint(tau="+std::to_string(tau)+",g[0]="+std::to_string(g[0])+",g[1]="+std::to_string(g[1])+",g[2]="+std::to_string(g[2])+")",info,planner,object),_obj(obj),_g(g),_tau(tau)
{
  for(sizeType i=0; i<_planner.body().nrDOF(); i++) {
    DSSQPObjectiveComponent<T>::_gl.push_back(0);
    DSSQPObjectiveComponent<T>::_gu.push_back(0);
  }
}
template <typename T>
void ContactIndexConstraint<T>::addContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid,T mu,sizeType nDir)
{
  ContactIndexMap cm;
  cm._oid=oid;
  cm._oidOther=oidOther;
  cm._pid=pid;

  Vec3T normal;
  ASSERT_MSGV(cm._oid>=-1 && cm._oid<_planner.body().nrDOF()/6,"Invalid object id: -1<=%d<%d",cm._oid,_planner.body().nrDOF()/6)
  if(oid>=0) {
    cm._oid=cm._oid*2+2;  //remap to linkid
    ASSERT_MSGV(cm._pid>=0 && cm._pid<_planner.pnss()[cm._oid].first.cols(),"Invalid point id: 0<=%d<%d",cm._pid,_planner.pnss()[cm._oid].first.cols())
    cm._objPos=_planner.pnss()[cm._oid].first.col(cm._pid);
    normal=_planner.pnss()[cm._oid].second.col(cm._pid);
  } else {
    ASSERT_MSGV(cm._pid>=0 && cm._pid<_object.pss().cols(),"Invalid point id: 0<=%d<%d",cm._pid,_object.pss().cols())
    cm._objPos=_object.pss().col(cm._pid);
    normal=_object.nss().col(cm._pid);
  }
  ASSERT_MSGV(cm._oidOther>=0 && cm._oidOther<_planner.body().nrDOF()/6,"Invalid object id: 0<=%d<%d",cm._oidOther,_planner.body().nrDOF()/6)
  cm._oidOther=cm._oidOther*2+2;  //remap to linkid

  Mat3T frm;
  sizeType id=-1;
  frm.col(0)=-normal;
  frm.col(0)/=std::sqrt(frm.col(0).squaredNorm());
  for(sizeType i=0; i<3; i++)
    if(id==-1 || std::abs(frm.col(0)[i])<std::abs(frm.col(0)[id]))
      id=i;
  frm.col(1)=Vec3T::Unit(id).cross(frm.col(0));
  frm.col(1)/=std::sqrt(frm.col(1).squaredNorm());
  frm.col(2)=frm.col(0).cross(frm.col(1));

  if(nDir<=3) {
    cm._off.resize(3);
    cm._objDir.resize(3,3);
    cm._objDir=frm;
    cm._objDir.col(1)*=mu;
    cm._objDir.col(2)*=mu;
    cm._off[0]=_obj.addVar("ContactO"+std::to_string(cm._oid)+
                           "OO"+std::to_string(cm._oidOther)+
                           "P"+std::to_string(cm._pid)+
                           "D"+std::to_string(0),0,
                           DSSQPObjectiveCompound<T>::infty(),
                           DSSQPObjectiveCompound<T>::NEW_OR_EXIST)._id;
    cm._off[1]=_obj.addVar("ContactO"+std::to_string(cm._oid)+
                           "OO"+std::to_string(cm._oidOther)+
                           "P"+std::to_string(cm._pid)+
                           "D"+std::to_string(1),
                           -DSSQPObjectiveCompound<T>::infty(),
                           DSSQPObjectiveCompound<T>::infty(),
                           DSSQPObjectiveCompound<T>::NEW_OR_EXIST)._id;
    cm._off[2]=_obj.addVar("ContactO"+std::to_string(cm._oid)+
                           "OO"+std::to_string(cm._oidOther)+
                           "P"+std::to_string(cm._pid)+
                           "D"+std::to_string(2),
                           -DSSQPObjectiveCompound<T>::infty(),
                           DSSQPObjectiveCompound<T>::infty(),
                           DSSQPObjectiveCompound<T>::NEW_OR_EXIST)._id;
    _obj.addQCone(cm._off);
  } else {
    cm._off.resize(nDir);
    cm._objDir.resize(3,nDir);
    for(sizeType i=0; i<nDir; i++) {
      T theta=M_PI*2*i/nDir;
      cm._off[i]=_obj.addVar("ContactO"+std::to_string(cm._oid)+
                             "OO"+std::to_string(cm._oidOther)+
                             "P"+std::to_string(cm._pid)+
                             "D"+std::to_string(i),0,
                             DSSQPObjectiveCompound<T>::infty(),
                             DSSQPObjectiveCompound<T>::NEW_OR_EXIST)._id;
      Vec3T nt(1,std::cos(theta)*mu,std::sin(theta)*mu);
      cm._objDir.col(i)=frm*nt;
    }
  }

  typename std::vector<ContactIndexMap>::iterator it=std::lower_bound(_cmss.begin(),_cmss.end(),cm);
  if(it==_cmss.end())
    _cmss.insert(it,cm);
  else if(*it!=cm)
    _cmss.insert(it,cm);
  else *it=cm;

  for(sizeType i=0; i<2; i++) {
    DSSQPObjectiveComponent<T>::_gl.push_back(0);
    DSSQPObjectiveComponent<T>::_gu.push_back(DSSQPObjective<T>::infty());
  }
}
template <typename T>
bool ContactIndexConstraint<T>::existContactIndexMap(sizeType oid,sizeType oidOther,sizeType pid) const
{
  ContactIndexMap cm;
  cm._oid=oid;
  cm._oidOther=oidOther;
  cm._pid=pid;

  ASSERT_MSGV(cm._oid>=-1 && cm._oid<_planner.body().nrDOF()/6,"Invalid object id: -1<=%d<%d",cm._oid,_planner.body().nrDOF()/6)
  if(oid>=0) {
    cm._oid=cm._oid*2+2;  //remap to linkid
    ASSERT_MSGV(cm._pid>=0 && cm._pid<_planner.pnss()[cm._oid].first.cols(),"Invalid point id: 0<=%d<%d",cm._pid,_planner.pnss()[cm._oid].first.cols())
  } else {
    ASSERT_MSGV(cm._pid>=0 && cm._pid<_object.pss().cols(),"Invalid point id: 0<=%d<%d",cm._pid,_object.pss().cols())
  }
  ASSERT_MSGV(cm._oidOther>=0 && cm._oidOther<_planner.body().nrDOF()/6,"Invalid object id: 0<=%d<%d",cm._oidOther,_planner.body().nrDOF()/6)
  cm._oidOther=cm._oidOther*2+2;  //remap to linkid

  typename std::vector<ContactIndexMap>::const_iterator it=std::lower_bound(_cmss.begin(),_cmss.end(),cm);
  return it!=_cmss.end() && *it==cm;
}
template <typename T>
void ContactIndexConstraint<T>::clearContactIndexMap()
{
  _cmss.clear();
}
template <typename T>
void ContactIndexConstraint<T>::writeContactVTK(const Vec& x,const PBDArticulatedGradientInfo<T>& info,const std::string& path) const
{
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>> pss,fss;
  for(const ContactIndexMap& cm:_cmss) {
    Vec coef=Vec::Zero(cm._off.size());
    for(sizeType i=0; i<coef.size(); i++)
      coef[i]=x[cm._off[i]];

    Vec3T f=cm._objDir*coef;
    Vec3T p=cm._objPos;
    if(cm._oid>=0)
      p=ROTI(info._TM,cm._oid)*p+CTRI(info._TM,cm._oid);

    const Environment<T>& env=_planner.env(cm._oidOther);
    Vec3T pOther=ROTI(info._TM,cm._oidOther).transpose()*(p-CTRI(info._TM,cm._oidOther)),normal;
    T phi=env.phi(pOther,&normal);
    pOther-=normal*phi;
    pOther=ROTI(info._TM,cm._oidOther)*pOther+CTRI(info._TM,cm._oidOther);

    sizeType off=(sizeType)vss.size();
    vss.push_back(p.unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
    vss.push_back((p+f).unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
    vss.push_back(pOther.unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
    vss.push_back((pOther-f).unaryExpr([&](const T& in) {
      return (scalar)std::to_double(in);
    }));
    pss.push_back(Vec3i::Constant(off+0));
    pss.push_back(Vec3i::Constant(off+2));
    fss.push_back(Vec3i(0,1,0)+Vec3i::Constant(off));
    fss.push_back(Vec3i(2,3,0)+Vec3i::Constant(off));
    fss.push_back(Vec3i(0,2,0)+Vec3i::Constant(off));
  }
  VTKWriter<scalar> os("force",path,true);
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(pss.begin(),pss.end(),VTKWriter<scalar>::POINT);
  os.appendCells(fss.begin(),fss.end(),VTKWriter<scalar>::LINE);
}
//constraints
template <typename T>
int ContactIndexConstraint<T>::operator()(const Vec& x,Vec& fvec,STrips* fjac)
{
  sizeType nDOF=_planner.body().nrDOF();
  fvec.segment(_offset,values()).setZero();
  for(sizeType i=0,off=_offset; i<nDOF/6; i++,off+=6)
    fvec.template segment<3>(off)=_planner.body().joint(i*2+2)._M*_g;
  for(sizeType i=0,off=_offset+nDOF; i<(sizeType)_cmss.size(); i++,off+=2) {
    const ContactIndexMap& cm=_cmss[i];
    const Environment<T>& env=_planner.env(cm._oidOther);
    Vec3T normal,force,forceO,forceG,objDirO,pG;
    if(cm._oid>=0)
      pG=ROTI(_info._TM,cm._oid)*cm._objPos+CTRI(_info._TM,cm._oid);
    else pG=cm._objPos;
    Vec3T pGL=ROTI(_info._TM,cm._oidOther).transpose()*(pG-CTRI(_info._TM,cm._oidOther));
    T dist=env.phi(pGL,fjac?&normal:NULL),xSum=0;

    force.setZero();
    for(sizeType vid=0; vid<cm._off.size(); vid++) {
      force+=cm._objDir.col(vid)*x[cm._off[vid]];
      xSum+=x[cm._off[vid]];
    }
    fvec[off+0]=dist;
    fvec[off+1]=_tau-xSum*dist;
#define FORCE_OFF _offset+(cm._oid-2)*3
#define TORQUE_OFF FORCE_OFF+3
#define FORCE_OTHER_OFF _offset+(cm._oidOther-2)*3
#define TORQUE_OTHER_OFF FORCE_OTHER_OFF+3
    if(cm._oid>=0) {
      forceG=ROTI(_info._TM,cm._oid)*force;
      fvec.template segment<3>(FORCE_OFF)+=force;
      fvec.template segment<3>(TORQUE_OFF)+=cm._objPos.cross(force);
    } else forceG=force;
    forceO=ROTI(_info._TM,cm._oidOther).transpose()*forceG;
    fvec.template segment<3>(FORCE_OTHER_OFF)-=forceO;
    fvec.template segment<3>(TORQUE_OTHER_OFF)-=pGL.cross(forceO);

    if(fjac) {
      Vec jac=Vec::Zero(nDOF);
      Mat3XT G=Mat3XT::Zero(3,_planner.body().nrJ()*4);
      if(cm._oid>=0) {
        ROTI(G,cm._oid)=(ROTI(_info._TM,cm._oidOther)*normal)*cm._objPos.transpose();
        CTRI(G,cm._oid)=ROTI(_info._TM,cm._oidOther)*normal;
      }
      ROTI(G,cm._oidOther)=(pG-CTRI(_info._TM,cm._oidOther))*normal.transpose();
      CTRI(G,cm._oidOther)=-ROTI(_info._TM,cm._oidOther)*normal;

      _info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(G),ArticulatedObjective<T>::mapV(jac));
      //jacobian of distance W.R.T. pose
      addBlock(*fjac,off+0,0,jac.transpose());
      addBlock(*fjac,off+1,0,-jac.transpose()*xSum);
      for(sizeType vid=0; vid<cm._off.size(); vid++) {
        //jacobian of complementarity W.R.T. dir
        fjac->push_back(STrip(off+1,cm._off[vid],-dist));
        //jacobian for force W.R.T. dir
        if(cm._oid>=0) {
          addBlock(*fjac,FORCE_OFF,cm._off[vid],cm._objDir.col(vid));
          addBlock(*fjac,TORQUE_OFF,cm._off[vid],cm._objPos.cross(cm._objDir.col(vid)));
          objDirO=ROTI(_info._TM,cm._oidOther).transpose()*(ROTI(_info._TM,cm._oid)*cm._objDir.col(vid));
        } else objDirO=ROTI(_info._TM,cm._oidOther).transpose()*cm._objDir.col(vid);
        addBlock(*fjac,FORCE_OTHER_OFF,cm._off[vid],-objDirO);
        addBlock(*fjac,TORQUE_OTHER_OFF,cm._off[vid],pGL.cross(-objDirO));
      }
      //jacobian of force W.R.T. pose
      for(sizeType d=0; d<3; d++) {
        G.setZero();
        jac.setZero();
        if(cm._oid>=0)
          ROTI(G,cm._oid)=-ROTI(_info._TM,cm._oidOther).col(d)*force.transpose();
        ROTI(G,cm._oidOther)=-forceG*Vec3T::Unit(d).transpose();
        _info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(G),ArticulatedObjective<T>::mapV(jac));
        addBlock(*fjac,FORCE_OTHER_OFF+d,0,jac.transpose());

        G.setZero();
        jac.setZero();
        if(cm._oid>=0) {
          ROTI(G,cm._oid)=-(pG-CTRI(_info._TM,cm._oidOther)).cross(ROTI(_info._TM,cm._oidOther).col(d))*force.transpose();
          ROTI(G,cm._oid)+=forceG.cross(ROTI(_info._TM,cm._oidOther).col(d))*cm._objPos.transpose();
          CTRI(G,cm._oid)=forceG.cross(ROTI(_info._TM,cm._oidOther).col(d));
        }
        ROTI(G,cm._oidOther)=(pG-CTRI(_info._TM,cm._oidOther)).cross(forceG)*Vec3T::Unit(d).transpose();
        CTRI(G,cm._oidOther)=-forceG.cross(ROTI(_info._TM,cm._oidOther).col(d));
        _info.DTG(_planner.body(),ArticulatedObjective<T>::mapM(G),ArticulatedObjective<T>::mapV(jac));
        addBlock(*fjac,TORQUE_OTHER_OFF+d,0,-jac.transpose());
      }
    }
#undef FORCE_OFF
#undef TORQUE_OFF
#undef FORCE_OTHER_OFF
#undef TORQUE_OTHER_OFF
  }
  return 0;
}
template <typename T>
int ContactIndexConstraint<T>::values() const
{
  return _planner.body().nrDOF()+(sizeType)_cmss.size()*2;
}
//instance
PRJ_BEGIN
template class ContactIndexConstraint<double>;
#ifdef ALL_TYPES
template class ContactIndexConstraint<__float128>;
template class ContactIndexConstraint<mpfr::mpreal>;
#endif
PRJ_END
