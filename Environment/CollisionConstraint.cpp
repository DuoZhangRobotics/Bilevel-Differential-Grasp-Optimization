#ifdef OPTIMIZER_SUPPORT
#ifdef ENVIRONMENT_SUPPORT
#ifdef DEFORMABLE_SUPPORT
#include <Utils/Scalar.h>
#include "CollisionConstraint.h"
#include "Simulator.h"
#include <Deformable/FEMMesh.h>
#include <Deformable/FEMCell.h>
#include <Deformable/FEMSystem.h>
#include <Deformable/FEMGradientInfo.h>
#include <Articulated/PBDArticulatedGradientInfo.h>
#include <CommonFile/geom/BVHBuilder.h>

USE_PRJ_NAMESPACE

//CollisionConstraint
template <typename T>
void CollisionConstraint<T>::buildFrame(sizeType n)
{
  Vec3T t1,t2;
  _frame.resize(3,n+1);
  _frame.col(0)=(_pR-_pL).template cast<T>();
  _frame.col(0)/=std::sqrt(_frame.col(0).squaredNorm());

  sizeType id;
  _frame.col(0).unaryExpr([&](const T& val) {
    return std::abs(val);
  }).minCoeff(&id);
  t2.setUnit(id);
  t1=_frame.col(0).cross(t2);
  t1/=std::sqrt(t1.squaredNorm());
  t2=_frame.col(0).cross(t1);
  t2/=std::sqrt(t2.squaredNorm());

  for(sizeType i=1; i<=n; i++) {
    T theta=M_PI*2*(i-1)/n;
    _frame.col(i)=t1*std::cos(theta)+t2*std::sin(theta);
  }
}
template <typename T>
T CollisionConstraint<T>::assembleSystem(const SimulatorObject<T>& objL,const SimulatorObject<T>& objR,sizeType row,STrips* J) const
{
  T val=0;
  if(_isP2T) {
    //P2T
    val+=_L->assembleSystem(objL,row,J,_vid[0],-_bary[0]*_frame.col(0));
    val+=_R->assembleSystem(objR,row,J,_vid[1], _bary[1]*_frame.col(0));
    val+=_R->assembleSystem(objR,row,J,_vid[2], _bary[2]*_frame.col(0));
    val+=_R->assembleSystem(objR,row,J,_vid[3], _bary[3]*_frame.col(0));
  } else {
    //E2E
    val+=_L->assembleSystem(objL,row,J,_vid[0],-_bary[0]*_frame.col(0));
    val+=_L->assembleSystem(objL,row,J,_vid[1],-_bary[1]*_frame.col(0));
    val+=_R->assembleSystem(objR,row,J,_vid[2], _bary[2]*_frame.col(0));
    val+=_R->assembleSystem(objR,row,J,_vid[3], _bary[3]*_frame.col(0));
  }
  return val;
}
//CollisionNode
template <typename T>
CollisionNode<T>::CollisionNode(const SimulatorObject<T>& obj,sizeType objId,sizeType jointId)
{
  _objId=objId;
  ObjMesh mesh;
  if(obj._sysFEM) {
    _jointId=-1;
    ASSERT_MSG(obj._sysFEM->mesh().nrB()==1,"We only support FEMMesh with one body")
    const FEMBody<T>& body=obj._sysFEM->mesh().getB(0);
    ASSERT_MSG(body.nrEM()==1,"We only support FEMBody with one embedded mesh")
    body.writeObjEmbedded(NULL,mesh,0);
    _vlss.resize(mesh.getV().size());
    for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++)
      _vlss[i]=mesh.getV(i).template cast<scalarD>();
  } else {
    _jointId=jointId;
    ASSERT(obj._body)
    obj._body->joint(jointId).getGeomPtr()->getMesh(mesh);
    _vlss.resize(mesh.getV().size());
    for(sizeType i=0; i<(sizeType)mesh.getV().size(); i++)
      _vlss[i]=mesh.getV(i).template cast<scalarD>();
  }

  ObjMesh::EdgeMap eMap;
  mesh.buildEdge(eMap);
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>> ess;
  for(auto E:eMap._ess)
    ess.push_back(Vec2i(E.first.first,E.first.second));
  buildBVH(mesh.getI(),ess);
}
template <typename T>
void CollisionNode<T>::updateBVH(const SimulatorObject<T>& obj,sizeType id)
{
  _vgss=_vlss;
  if(_jointId==-1) {
    const FEMBody<T>& body=obj._sysFEM->mesh().getB(0);
    for(sizeType i=0; i<(sizeType)_vgss.size(); i++) {
      Vec3T v=obj._xFEM[id].getVert(*(obj._sysFEM),body.getEM(0)._emvss[i]);
      _vgss[i]=v.unaryExpr([&](const T& val) {
        return (scalarD)std::to_double(val);
      });
    }
  } else {
    _vgss=_vlss;
    Mat3X4d TM=TRANSI(obj._xPBD[id]._TM,_jointId).unaryExpr([&](const T& val) {
      return (scalarD)std::to_double(val);
    });
    for(sizeType i=0; i<(sizeType)_vgss.size(); i++)
      _vgss[i]=ROT(TM)*_vlss[i]+CTR(TM);
  }
  updateBVH();
}
//check collision
template <typename T>
bool CollisionNode<T>::checkCollisionE2T(std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R)
{
  scalarD tmpS;
  std::stack<std::pair<sizeType,sizeType>> ss;
  ss.push(std::make_pair((sizeType)L->_bvhE.size()-1,(sizeType)R->_bvhT.size()-1));
  while(!ss.empty()) {
    sizeType pid=ss.top().first;
    sizeType tid=ss.top().second;
    const Node<Vec2i,BBOX>& eNode=L->_bvhE[pid];
    const Node<Vec3i,BBOX>& tNode=R->_bvhT[tid];
    ss.pop();
    if(!eNode._bb.intersect(tNode._bb))
      continue;
    else if(eNode._cell[0]>=0 && tNode._cell[0]>=0) {
      const LineSegTpl<scalarD> e(L->_vgss[eNode._cell[0]],L->_vgss[eNode._cell[1]]);
      const TriangleTpl<scalarD> t(R->_vgss[tNode._cell[0]],R->_vgss[tNode._cell[1]],R->_vgss[tNode._cell[2]]);
      if(t.intersect(e,tmpS,false))
        return true;
    } else if(eNode._cell[0]>=0) {
      ss.push(std::make_pair(pid,tNode._l));
      ss.push(std::make_pair(pid,tNode._r));
    } else if(tNode._cell[0]>=0) {
      ss.push(std::make_pair(eNode._l,tid));
      ss.push(std::make_pair(eNode._r,tid));
    } else {
      ss.push(std::make_pair(eNode._l,tNode._l));
      ss.push(std::make_pair(eNode._r,tNode._l));
      ss.push(std::make_pair(eNode._l,tNode._r));
      ss.push(std::make_pair(eNode._r,tNode._r));
    }
  }
  return false;
}
template <typename T>
bool CollisionNode<T>::checkCollision(std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R)
{
  return checkCollisionE2T(L,R)||checkCollisionE2T(R,L);
}
//add constraint
template <typename T>
void CollisionNode<T>::addConstraintP2T(std::vector<CollisionConstraint<T>>& cons,std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R,T dist)
{
  if(L->_bvhV.empty())
    return;
  if(R->_bvhT.empty())
    return;

  CollisionConstraint<T> c;
  c._L=L;
  c._R=R;
  c._isP2T=true;

  std::stack<std::pair<sizeType,sizeType>> ss;
  ss.push(std::make_pair((sizeType)L->_bvhV.size()-1,(sizeType)R->_bvhT.size()-1));
  while(!ss.empty()) {
    sizeType pid=ss.top().first;
    sizeType tid=ss.top().second;
    const Node<sizeType,BBOX>& pNode=L->_bvhV[pid];
    const Node<Vec3i,BBOX>& tNode=R->_bvhT[tid];
    ss.pop();
    if(!pNode._bb.enlarge(std::to_double(dist)).intersect(tNode._bb))
      continue;
    else if(pNode._cell>=0 && tNode._cell[0]>=0) {
      const Vec3d& p=L->_vgss[pNode._cell];
      const TriangleTpl<scalarD> t(R->_vgss[tNode._cell[0]],R->_vgss[tNode._cell[1]],R->_vgss[tNode._cell[2]]);
      //constraint
      Vec3d cp,b;
      scalarD sqrDist=ScalarUtil<scalarD>::scalar_max();
      t.calcPointDist(p,sqrDist,cp,b);
      if(sqrDist<dist*dist) {
        c._vid=Vec4i(pNode._cell,tNode._cell[0],tNode._cell[1],tNode._cell[2]);
        c._bary=Vec4d(1,b[0],b[1],b[2]);
        c._pL=p;
        c._pR=cp;
        cons.push_back(c);
      }
    } else if(pNode._cell>=0) {
      ss.push(std::make_pair(pid,tNode._l));
      ss.push(std::make_pair(pid,tNode._r));
    } else if(tNode._cell[0]>=0) {
      ss.push(std::make_pair(pNode._l,tid));
      ss.push(std::make_pair(pNode._r,tid));
    } else {
      ss.push(std::make_pair(pNode._l,tNode._l));
      ss.push(std::make_pair(pNode._r,tNode._l));
      ss.push(std::make_pair(pNode._l,tNode._r));
      ss.push(std::make_pair(pNode._r,tNode._r));
    }
  }
}
template <typename T>
void CollisionNode<T>::addConstraintE2E(std::vector<CollisionConstraint<T>>& cons,std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R,T dist)
{
  if(L->_bvhE.empty())
    return;
  if(R->_bvhE.empty())
    return;

  CollisionConstraint<T> c;
  c._L=L;
  c._R=R;
  c._isP2T=false;

  std::stack<std::pair<sizeType,sizeType>> ss;
  ss.push(std::make_pair((sizeType)L->_bvhE.size()-1,(sizeType)R->_bvhE.size()-1));
  while(!ss.empty()) {
    sizeType e1id=ss.top().first;
    sizeType e2id=ss.top().second;
    const Node<Vec2i,BBOX>& e1Node=L->_bvhE[e1id];
    const Node<Vec2i,BBOX>& e2Node=R->_bvhE[e2id];
    ss.pop();
    if(!e1Node._bb.enlarge(std::to_double(dist)).intersect(e2Node._bb))
      continue;
    else if(e1Node._cell[0]>=0 && e2Node._cell[0]>=0) {
      const LineSegTpl<scalarD> e1(L->_vgss[e1Node._cell[0]],L->_vgss[e1Node._cell[1]]);
      const LineSegTpl<scalarD> e2(R->_vgss[e2Node._cell[0]],R->_vgss[e2Node._cell[1]]);
      //constraint
      scalarD a,b;
      scalarD sqrDist=ScalarUtil<scalarD>::scalar_max();
      e1.calcLineDist(e2,sqrDist,a,b);
      if(sqrDist<dist*dist) {
        c._vid=Vec4i(e1Node._cell[0],e1Node._cell[1],e2Node._cell[0],e2Node._cell[1]);
        c._bary=Vec4d(1-a,a,1-b,b);
        c._pL=e1._x*(1-a)+e1._y*a;
        c._pR=e2._x*(1-b)+e2._y*b;
        cons.push_back(c);
      }
    } else if(e1Node._cell[0]>=0) {
      ss.push(std::make_pair(e1id,e2Node._l));
      ss.push(std::make_pair(e1id,e2Node._r));
    } else if(e2Node._cell[0]>=0) {
      ss.push(std::make_pair(e1Node._l,e2id));
      ss.push(std::make_pair(e1Node._r,e2id));
    } else {
      ss.push(std::make_pair(e1Node._l,e2Node._l));
      ss.push(std::make_pair(e1Node._r,e2Node._l));
      ss.push(std::make_pair(e1Node._l,e2Node._r));
      ss.push(std::make_pair(e1Node._r,e2Node._r));
    }
  }
}
template <typename T>
void CollisionNode<T>::addConstraint(std::vector<CollisionConstraint<T>>& cons,std::shared_ptr<CollisionNode<T>> L,std::shared_ptr<CollisionNode<T>> R,T dist)
{
  addConstraintP2T(cons,L,R,dist);
  addConstraintP2T(cons,R,L,dist);
  addConstraintE2E(cons,L,R,dist);
  addConstraintE2E(cons,R,L,dist);
}
template <typename T>
T CollisionNode<T>::assembleSystem(const SimulatorObject<T>& obj,sizeType row,STrips* J,sizeType id,const Vec3T& coef) const
{
  if(_jointId>=0) {
    //this is an articulated object
    const PBDArticulatedGradientInfo<T>& info=obj._xPBD[0];
    if(J) {
      info.JRCSparse(*(obj._body),_jointId,[&](sizeType col,const Vec3T& w) {
        J->push_back(STrip(row,obj._offDOF+col,w.cross(ROTI(info._TM,_jointId)*_vlss[id].template cast<T>()).dot(coef)));
      },[&](sizeType col,const Vec3T& w) {
        J->push_back(STrip(row,obj._offDOF+col,w.dot(coef)));
      });
    }
    return (ROTI(info._TM,_jointId)*_vlss[id].template cast<T>()+CTRI(info._TM,_jointId)).dot(coef);
  } else {
    //this is an deformable object
    const FEMGradientInfo<T>& info=obj._xFEM[0];
    const FEMInterp<T>& I=obj._sysFEM->mesh().getB(0).getEM(0)._emvss[id];
    if(J) {
      info.JSparse(*(obj._sysFEM),I,[&](sizeType col,const Vec3T& w) {
        J->push_back(STrip(row,obj._offDOF+col,w.dot(coef)));
      });
    }
    return info.getVert(*(obj._sysFEM),I).dot(coef);
  }
}
template <typename T>
BBox<scalarD> CollisionNode<T>::getBB() const
{
  return BBox<scalarD>(_bvhV.back()._bb.minCorner(),_bvhV.back()._bb.maxCorner());
}
template <typename T>
sizeType CollisionNode<T>::objId() const
{
  return _objId;
}
template <typename T>
sizeType CollisionNode<T>::jointId() const
{
  return _jointId;
}
template <typename T>
sizeType CollisionNode<T>::nrV() const
{
  return (sizeType)_vlss.size();
}
//helper
template <typename T>
void CollisionNode<T>::updateBVH()
{
  //bvhV
  for(sizeType i=0; i<((sizeType)_bvhV.size()+1)/2; i++) {
    _bvhV[i]._bb.reset();
    _bvhV[i]._bb.setUnion(_vgss[_bvhV[i]._cell]);
  }
  BVHQuery<sizeType,BBOX>(_bvhV,3,-1).updateBVH();
  //bvhE
  for(sizeType i=0; i<((sizeType)_bvhE.size()+1)/2; i++) {
    _bvhE[i]._bb.reset();
    for(sizeType d=0; d<2; d++)
      _bvhE[i]._bb.setUnion(_vgss[_bvhE[i]._cell[d]]);
  }
  BVHQuery<Vec2i,BBOX>(_bvhE,3,Vec2i(-1,-1)).updateBVH();
  //bvhT
  for(sizeType i=0; i<((sizeType)_bvhT.size()+1)/2; i++) {
    _bvhT[i]._bb.reset();
    for(sizeType d=0; d<3; d++)
      _bvhT[i]._bb.setUnion(_vgss[_bvhT[i]._cell[d]]);
  }
  BVHQuery<Vec3i,BBOX>(_bvhT,3,Vec3i(-1,-1,-1)).updateBVH();
}
template <typename T>
void CollisionNode<T>::buildBVH(const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>>& iss,
                                const std::vector<Vec2i,Eigen::aligned_allocator<Vec2i>>& ess)
{
  //bvhV
  _bvhV.assign(_vlss.size(),Node<sizeType,BBOX>());
  for(sizeType i=0; i<(sizeType)_vlss.size(); i++) {
    _bvhV[i]._nrCell=1;
    _bvhV[i]._cell=i;
    _bvhV[i]._bb.reset();
    _bvhV[i]._bb.setUnion(_vlss[i]);
  }
  COMMON::buildBVH<sizeType,BBOX>(_bvhV,3,-1);
  //bvhE
  _bvhE.assign(ess.size(),Node<Vec2i,BBOX>());
  for(sizeType i=0; i<(sizeType)ess.size(); i++) {
    _bvhE[i]._nrCell=1;
    _bvhE[i]._cell=ess[i];
    for(sizeType d=0; d<2; d++)
      _bvhE[i]._bb.setUnion(_vlss[ess[i][d]]);
  }
  COMMON::buildBVH<Vec2i,BBOX>(_bvhE,3,Vec2i(-1,-1));
  //bvhT
  _bvhT.assign(iss.size(),Node<Vec3i,BBOX>());
  for(sizeType i=0; i<(sizeType)iss.size(); i++) {
    _bvhT[i]._nrCell=1;
    _bvhT[i]._cell=iss[i];
    _bvhT[i]._bb.reset();
    for(sizeType d=0; d<3; d++)
      _bvhT[i]._bb.setUnion(_vlss[iss[i][d]]);
  }
  COMMON::buildBVH<Vec3i,BBOX>(_bvhT,3,Vec3i(-1,-1,-1));
}
//instance
PRJ_BEGIN
template class CollisionConstraint<double>;
template class CollisionNode<double>;

#ifdef ALL_TYPES
template class CollisionConstraint<__float128>;
template class CollisionNode<__float128>;

template class CollisionConstraint<mpfr::mpreal>;
template class CollisionNode<mpfr::mpreal>;
#endif
PRJ_END
#endif
#endif
#endif
