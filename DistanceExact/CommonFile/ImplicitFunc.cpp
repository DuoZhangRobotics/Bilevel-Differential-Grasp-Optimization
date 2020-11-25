#include "ImplicitFunc.h"
#include "GridOp.h"
#include "geom/StaticGeomCell.h"

USE_PRJ_NAMESPACE

//ImplicitFuncCSG
ImplicitFuncCSG::ImplicitFuncCSG(OP_TYPE op):_a((ImplicitFuncCSG*)NULL),_b((ImplicitFuncCSG*)NULL),_alpha(0.8f),_op(op) {}
ImplicitFuncCSG::ImplicitFuncCSG(const BBox<scalar,3>& bb,OP_TYPE op):_alpha(0.8f),_op(op)
{
  BBox<scalar,2> bb2(bb._minC.block(0,0,2,1),bb._maxC.block(0,0,2,1));
  std::shared_ptr<ImplicitFuncCSG> axis01(new ImplicitFuncCSG(bb2,op));
  std::shared_ptr<ImplicitFuncCSG> axis2(new ImplicitFuncCSG(op));
  if(op == INTERSECT) {
    axis2->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._minC.z()),Vec3(0.0f,0.0f,-1.0f)));
    axis2->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._maxC.z()),Vec3(0.0f,0.0f, 1.0f)));
  } else {
    axis2->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._minC.z()),Vec3(0.0f,0.0f, 1.0f)));
    axis2->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,0.0f,bb._maxC.z()),Vec3(0.0f,0.0f,-1.0f)));
  }
  _a=axis01;
  _b=axis2;
}
ImplicitFuncCSG::ImplicitFuncCSG(const BBox<scalar,2>& bb,OP_TYPE op):_alpha(0.8f)
{
  std::shared_ptr<ImplicitFuncCSG> axis0(new ImplicitFuncCSG(op));
  std::shared_ptr<ImplicitFuncCSG> axis1(new ImplicitFuncCSG(op));
  if(op == INTERSECT) {
    axis0->_a.reset(new ImplicitFuncPlane(Vec3(bb._minC.x(),0.0f,0.0f),Vec3(-1.0f,0.0f,0.0f)));
    axis0->_b.reset(new ImplicitFuncPlane(Vec3(bb._maxC.x(),0.0f,0.0f),Vec3( 1.0f,0.0f,0.0f)));
    axis1->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._minC.y(),0.0f),Vec3(0.0f,-1.0f,0.0f)));
    axis1->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._maxC.y(),0.0f),Vec3(0.0f, 1.0f,0.0f)));
  } else {
    axis0->_a.reset(new ImplicitFuncPlane(Vec3(bb._minC.x(),0.0f,0.0f),Vec3( 1.0f,0.0f,0.0f)));
    axis0->_b.reset(new ImplicitFuncPlane(Vec3(bb._maxC.x(),0.0f,0.0f),Vec3(-1.0f,0.0f,0.0f)));
    axis1->_a.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._minC.y(),0.0f),Vec3(0.0f, 1.0f,0.0f)));
    axis1->_b.reset(new ImplicitFuncPlane(Vec3(0.0f,bb._maxC.y(),0.0f),Vec3(0.0f,-1.0f,0.0f)));
  }
  _a=axis0;
  _b=axis1;
}
scalar ImplicitFuncCSG::operator()(const Vec3& pos) const
{
  if(!_a && !_b)
    return 1.0f;
  else if(!_b)
    return (*_a)(pos);
  else if(!_a) {
    scalar vb=(*_b)(pos);
    if(_op == SUBTRACT)
      vb*=-1.0f;
    return vb;
  } else {
    scalar va=(*_a)(pos);
    scalar vb=(*_b)(pos);
    if(_op == SUBTRACT)
      vb*=-1.0f;

    scalar sgn=_op == UNION ? -1.0f : 1.0f;
    return (va+vb+sgn*std::sqrt(std::abs(va*va+vb*vb-2.0f*va*vb*_alpha)))/(1.0f+_alpha);
  }
}
void ImplicitFuncCSG::setAlpha(const scalar& alpha)
{
  _alpha=alpha;
  ImplicitFuncCSG* a=dynamic_cast<ImplicitFuncCSG*>(_a.get());
  ImplicitFuncCSG* b=dynamic_cast<ImplicitFuncCSG*>(_b.get());
  if(a)a->setAlpha(alpha);
  if(b)b->setAlpha(alpha);
}
//ImplicitFuncPlane
ImplicitFuncPlane::ImplicitFuncPlane() {}
ImplicitFuncPlane::ImplicitFuncPlane(const Vec3& x0,const Vec3& n):_p(x0,n) {}
ImplicitFuncPlane::ImplicitFuncPlane(const Vec3& a,const Vec3& b,const Vec3& c):_p(a,b,c) {}
scalar ImplicitFuncPlane::operator()(const Vec3& pos) const
{
  return _p.side(pos);
}
//ImplicitFuncGridRef
scalar ImplicitFuncGridRef::operator()(const Vec3& pos) const
{
  return _lsRef.sampleSafe(pos);
}
//ImplicitFuncReinit
ImplicitFuncReinit::ImplicitFuncReinit(scalar cellSz,const ImplicitFunc<scalar>& inner)
{
  BBox<scalar> bb=inner.getBB();
  bb.enlarged(cellSz*3,bb.getExtent()[2] == 0 ? 2 : 3);
  Vec3i nrCell=ceilV(Vec3(bb.getExtent()/cellSz));
  _ls.reset(nrCell,bb,0.0f);

  GridOp<scalar,scalar>::copyFromImplictFunc(_ls,inner);
  GridOp<scalar,scalar>::reinitialize(_ls);
}
ImplicitFuncReinit::ImplicitFuncReinit(scalar cellSz,ImplicitFunc<scalar>& inner,bool verbose)
{
  BBox<scalar> bb=inner.getBB();
  bb.enlarged(cellSz*3,bb.getExtent()[2] == 0 ? 2 : 3);
  Vec3i nrCell=ceilV(Vec3(bb.getExtent()/cellSz));
  _ls.reset(nrCell,bb,0.0f);

  GridOp<scalar,scalar>::copyFromImplictFuncCached(_ls,inner);
  GridOp<scalar,scalar>::reinitialize(_ls);
}
ImplicitFuncReinit::ImplicitFuncReinit(const Grid<scalar,scalar> &tpl,const ImplicitFunc<scalar>& inner)
{
  _ls.makeSameGeometry(tpl);
  GridOp<scalar,scalar>::copyFromImplictFunc(_ls,inner);
  GridOp<scalar,scalar>::reinitialize(_ls);
}
scalar ImplicitFuncReinit::operator()(const Vec3& pos) const
{
  return _ls.sampleSafe(pos);
}
//ImplicitFuncMesh
ImplicitFuncMeshRef::ImplicitFuncMeshRef(const ObjMeshGeomCell& mesh,scalar eps)
  :_mesh(mesh),_eps(eps),_ext(mesh.getBB().getExtent()*2) {}
void ImplicitFuncMeshRef::beginSampleSet(const ScalarField& toBeSampled)
{
  _tag.reset(new TagField);
  _tag->makeSameGeometry(toBeSampled);
  _tag->init(0);
  OMP_PARALLEL_FOR_
  for(sizeType y=0; y<_tag->getNrPoint()[1]; y++)
    for(sizeType z=0; z<_tag->getNrPoint()[2]; z++)
      search(Vec3i(0,y,z),0);
  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<_tag->getNrPoint()[0]; x++)
    for(sizeType z=0; z<_tag->getNrPoint()[2]; z++)
      search(Vec3i(x,0,z),1);
  if(_tag->getDim() == 3) {
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<_tag->getNrPoint()[0]; x++)
      for(sizeType y=0; y<_tag->getNrPoint()[1]; y++)
        search(Vec3i(x,y,0),2);
  }
}
void ImplicitFuncMeshRef::endSampleSet(const ScalarField& toBeSampled)
{
  if(_tag)_tag.reset((TagField*)NULL);
}
scalar ImplicitFuncMeshRef::operator()(const Vec3& pos) const
{
  Vec3 n;
  bool inside=_mesh.closest(pos,n);
  scalar c=n.norm();

  //if tag is built beforehand
  if(_tag) {
    if(_tag->sampleSafe(pos) > 0.5f)
      return c;
    else return c*(inside ? -1 : 1);
  }

  //ASSERT_MSG(false,"Not implemented, used cached version!")
  /*//brute force visibility test for enhanced in-out testing
  else if(_mesh.rayQuery(pos,-Vec3::Unit(0)*_ext[0]) == 1)
    return c;
  else if(_mesh.rayQuery(pos, Vec3::Unit(0)*_ext[0]) == 1)
    return c;
  else if(_mesh.rayQuery(pos,-Vec3::Unit(1)*_ext[1]) == 1)
    return c;
  else if(_mesh.rayQuery(pos, Vec3::Unit(1)*_ext[1]) == 1)
    return c;
  else if(_mesh.dim() == 3) {
    if(_mesh.rayQuery(pos,-Vec3::Unit(2)*_ext[2]) == 1)
      return c;
    else if(_mesh.rayQuery(pos, Vec3::Unit(2)*_ext[2]) == 1)
      return c;
  }*/
  return c*(inside ? -1 : 1);
}
BBox<scalar> ImplicitFuncMeshRef::getBB() const
{
  return _mesh.getBB().enlargeEps(_eps);
}
void ImplicitFuncMeshRef::search(Vec3i base,sizeType dir)
{
  //forward
  Vec3i currI=base;
  Vec3i stepI=Vec3i::Unit(dir);
  Vec3 currP=_tag->getPt(base);
  Vec3 stepP=Vec3::Unit(dir)*_tag->getCellSize()[dir];
  Vec3 ray=Vec3::Unit(dir)*_ext[dir];
  while(currI[dir]<_tag->getNrPoint()[dir]) {
    scalar frac=_mesh.rayQuery(currP,ray);
    scalar lmt=currP[dir]+_ext[dir]*frac;
    while(currP[dir] <= lmt && currI[dir]<_tag->getNrPoint()[dir]) {
      if(frac == 1)
        _tag->get(currI)=1;
      currP+=stepP;
      currI+=stepI;
    }
  }

  //reverse
  stepI*=-1;
  stepP*=-1;
  ray*=-1;
  currP+=stepP;
  currI+=stepI;
  while(currI[dir]>=0) {
    scalar frac=_mesh.rayQuery(currP,ray);
    scalar lmt=currP[dir]-_ext[dir]*frac;
    while(currP[dir] >= lmt && currI[dir]>=0) {
      if(frac == 1)
        _tag->get(currI)=1;
      currP+=stepP;
      currI+=stepI;
    }
  }
}
//ImplicitFuncRosy
ImplicitFuncRosy::ImplicitFuncRosy(const Vec3& origin,const Vec3& X,const scalar& step,const scalar& coef)
  :_step(step),_coef(coef),_origin(origin),_X(X) {}
scalar ImplicitFuncRosy::operator()(const Vec3& pos) const
{
  Vec3 rel=pos-_origin;
  return _coef*dist(Vec2(rel.dot(_X),(rel-rel.dot(_X)*_X).norm()));
}
scalar ImplicitFuncRosy::dist(const Vec2& p) const
{
  //scalar minX=p.x();
  scalar dist=ScalarUtil<scalar>::scalar_max();
  for(scalar currX=p.x();; currX-=_step) {
    scalar newDist=(Vec2(currX,y(currX))-p).norm();
    if(newDist < dist) {
      dist=newDist;
      //minX=currX;
    } else break;
  }
  for(scalar currX=p.x();; currX+=_step) {
    scalar newDist=(Vec2(currX,y(currX))-p).norm();
    if(newDist < dist) {
      dist=newDist;
      //minX=currX;
    } else break;
  }
  return (y(p.x()) < p.y()) ? dist : -dist;
}
