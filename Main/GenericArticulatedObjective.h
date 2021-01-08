#ifndef GENERIC_ARTICULATED_OBJECTIVE_H
#define GENERIC_ARTICULATED_OBJECTIVE_H

#include <Articulated/ArticulatedBody.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <Utils/DebugGradient.h>

PRJ_BEGIN

template <typename T>
class GenericArticulatedObjective
{
public:
  DECL_MAP_FUNCS
  DECL_MAP_TYPES_T
  GenericArticulatedObjective(const ArticulatedBody& body):_body(body) {}
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g=NULL,MatT* h=NULL)=0;
  virtual int operator()(const Vec& x,T& e,Vec* g=NULL,MatT* h=NULL)
  {
    PBDArticulatedGradientInfo<T> info(_body,x.segment(0,_body.nrDOF()));
    Mat3XT G,tmpG;
    MatT H,tmpH;
    if(g)
      G.setZero(3,_body.nrJ()*4);
    if(h)
      H.setZero(12,_body.nrJ()*12);
    int ret=operator()(info,e,g?&G:NULL,h?&H:NULL);
    if(ret<0)
      return ret;
    if(g)
      info.DTG(_body,mapM(tmpG=G),mapV(*g));
    if(h)
      info.toolAB(_body,mapCM(tmpH=H),mapM(tmpG=G),mapM(h));
    return ret;
  }
  void debug(const Vec& x)
  {
    sizeType nDOF=_body.nrDOF();
    DEFINE_NUMERIC_DELTA_T(T)
    T e=0,e2=0;
    MatT h=MatT::Zero(nDOF,nDOF);
    Vec g=Vec::Zero(nDOF),g2=Vec::Zero(nDOF);
    Vec dx=Vec::Random(nDOF);
    if(operator()(x,e,&g,&h)<0) {
      std::ostringstream oss;
      for(sizeType i=0; i<x.size(); i++)
        oss << x[i] << " ";
      WARNINGV("Invalid debug position (x=%s)",oss.str().c_str())
      return;
    }
    operator()(x+dx*DELTA,e2,&g2);
    DEBUG_GRADIENT("GenericArticulatedObjective-G",g.dot(dx),g.dot(dx)-(e2-e)/DELTA)
    DEBUG_GRADIENT("GenericArticulatedObjective-H",std::sqrt((h*dx).squaredNorm()),std::sqrt((h*dx-(g2-g)/DELTA).squaredNorm()))
  }
protected:
  const ArticulatedBody& _body;
};
template <typename T>
class ExampleArticulatedObjective : public GenericArticulatedObjective<T>
{
public:
  DECL_MAP_TYPES_T
  using GenericArticulatedObjective<T>::_body;
  ExampleArticulatedObjective(const ArticulatedBody& body,sizeType N):GenericArticulatedObjective<T>(body)
  {
    for(sizeType i=0; i<N; i++) {
      _vLocal.push_back(Vec3T::Random());
      _H.push_back(Mat3T::Random());
      _H.back()=_H.back()*_H.back().transpose();
      _idLink.push_back(RandEngine::randI(0,_body.nrJ()-1));
    }
  }
  virtual int operator()(const PBDArticulatedGradientInfo<T>& info,T& e,Mat3XT* g=NULL,MatT* h=NULL)
  {
    for(sizeType i=0; i<(sizeType)_vLocal.size(); i++) {
      Vec3T vGlobal;
      Mat3X4T t=TRANSI(info._TM,_idLink[i]);
      vGlobal=ROT(t)*_vLocal[i]+CTR(t);
      e+=vGlobal.dot(_H[i]*vGlobal)/2;
      if(g) {
        ROTI(*g,_idLink[i])+=(_H[i]*vGlobal)*_vLocal[i].transpose();
        CTRI(*g,_idLink[i])+=_H[i]*vGlobal;
      }
      if(h) {
        Vec4T vLocalHomo(_vLocal[i][0],_vLocal[i][1],_vLocal[i][2],1);
        Eigen::Map<Eigen::Matrix<T,12,12>> hBlk(&(h->coeffRef(0,_idLink[i]*12)));
        for(sizeType r=0; r<4; r++)
          for(sizeType c=0; c<4; c++)
            hBlk.template block<3,3>(r*3,c*3)+=_H[i]*vLocalHomo[r]*vLocalHomo[c];
      }
    }
    return 0;
  }
protected:
  std::vector<Vec3T,Eigen::aligned_allocator<Vec3T>> _vLocal;
  std::vector<Mat3T,Eigen::aligned_allocator<Mat3T>> _H;
  std::vector<sizeType> _idLink;
};

PRJ_END

#endif
