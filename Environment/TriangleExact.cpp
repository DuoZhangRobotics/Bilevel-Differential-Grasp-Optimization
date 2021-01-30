#include "TriangleExact.h"
#include "MPQZIO.h"
#include <CommonFile/IO.h>
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

TriangleExact::TriangleExact():_vid(-1,-1,-1),_eNId(-1,-1,-1) {}
TriangleExact::TriangleExact(const TriangleTpl<scalar>& t):_vid(-1,-1,-1),_eNId(-1,-1,-1)
{
  _v[0]=t._a.cast<double>().cast<T>();
  _v[1]=t._b.cast<double>().cast<T>();
  _v[2]=t._c.cast<double>().cast<T>();
  reset();
}
TriangleExact::TriangleExact(const PT& a,const PT& b,const PT& c):_vid(-1,-1,-1),_eNId(-1,-1,-1)
{ 
  _v[0]=a;
  _v[1]=b;
  _v[2]=c;
  reset();
}
bool TriangleExact::read(std::istream& is,IOData*)
{
  readBinaryData(_vid,is);
  readBinaryData(_eNId,is);
  for(int d=0; d<3; d++) {
    readBinaryData(_v[d],is);
    readBinaryData(_n[d],is);
    readBinaryData(_nO[d],is);
  }
  readBinaryData(_nOSqr,is);
  readBinaryData(_nt,is);
  readBinaryData(_nTnN,is);
  readBinaryData(_invM,is);
  return is.good();
}
bool TriangleExact::write(std::ostream& os,IOData*) const
{
  writeBinaryData(_vid,os);
  writeBinaryData(_eNId,os);
  for(int d=0; d<3; d++) {
    writeBinaryData(_v[d],os);
    writeBinaryData(_n[d],os);
    writeBinaryData(_nO[d],os);
  }
  writeBinaryData(_nOSqr,os);
  writeBinaryData(_nt,os);
  writeBinaryData(_nTnN,os);
  writeBinaryData(_invM,os);
  return os.good();
}
std::shared_ptr<SerializableBase> TriangleExact::copy() const
{
  return std::shared_ptr<SerializableBase>(new TriangleExact);
}
std::string TriangleExact::type() const
{
  return typeid(TriangleExact).name();
}
const TriangleExact::PT& TriangleExact::v(sizeType i) const
{
  return _v[i];
}
const TriangleExact::PT& TriangleExact::normal() const
{
  return _nt;
}
TriangleExact::PT TriangleExact::bary(const PT& p) const
{
  MAT3X2 LHS;
  LHS.col(0)=_v[1]-_v[0];
  LHS.col(1)=_v[2]-_v[0];
  MAT2 LHS2=LHS.transpose()*LHS;
  PT2 ab=LHS2.inverse()*(LHS.transpose()*(p-_v[0]));
  return PT(1-ab.sum(),ab[0],ab[1]);
}
TriangleExact::PT TriangleExact::interp(const PT& bary) const
{
  return _v[0]*bary[0]+_v[1]*bary[1]+_v[2]*bary[2];
}
TriangleExact::PT2 TriangleExact::project(const PT& d) const
{
  T a=_v[0].dot(d);
  T b=_v[1].dot(d);
  T c=_v[2].dot(d);
  return PT2(std::min(a,std::min(b,c)),std::max(a,std::max(b,c)));
}
bool TriangleExact::intersect(const BBoxExact& bb) const
{
  std::function<bool(const PT2&,const PT2&)> testSegment=[&](const PT2& a,const PT2& b) {
    return a.x() > b.y() || b.x() > a.y();
  };
  std::function<bool(const PT&,const TriangleExact&,const BBoxExact&)> sep=[&](const PT& d,const TriangleExact& A,const BBoxExact& B) {
    return testSegment(A.project(d),B.project(d));
  };
  PT d=(_v[1]-_v[0]).cross(_v[2]-_v[0]);
  if(testSegment(bb.project(d),PT2::Constant(_v[0].dot(d))))
    return false;
  for(sizeType i=0; i<3; i++)
    if(sep(PT::Unit(i),*this,bb))
      return false;
  for(sizeType i=0; i<3; i++)
    for(sizeType j=0; j<3; j++)
      if(sep(PT::Unit(i).cross(_v[j]-_v[(j+1)%3]),*this,bb))
        return false;
  return true;
}
TriangleExact::T TriangleExact::dirGrad(const Vec2i& evid,const PT& D) const
{
  for(sizeType d=0; d<3; d++)
    if(_vid[d]!=evid[0] && _vid[d]!=evid[1]) {
      ASSERT((_vid[(d+1)%3]==evid[0] && _vid[(d+2)%3]==evid[1]) || (_vid[(d+2)%3]==evid[0] && _vid[(d+1)%3]==evid[1]))
      return -_n[d].dot(D);
    }
  ASSERT_MSGV(false,"Cannot find vertex opposite edge (%d,%d)",evid[0],evid[1])
  return 0;
}
std::pair<sizeType,TriangleExact::PT> TriangleExact::moveFromEdge(const Vec2i& evid,const PT& p0,PT D) const
{
  T t[2];
  sizeType vid=-1;
  for(sizeType d=0; d<3; d++)
    if(_vid[d]!=evid[0] && _vid[d]!=evid[1]) {
      ASSERT((_vid[(d+1)%3]==evid[0] && _vid[(d+2)%3]==evid[1]) || (_vid[(d+2)%3]==evid[0] && _vid[(d+1)%3]==evid[1]))
      vid=d;
      break;
    }
  D-=_nTnN*D;
  for(sizeType off=0; off<2; off++) {
    const PT& tp0=_v[(vid+off+2)%3];
    const PT& td=_nO[(vid+off+1)%3];
    PT2 RHS2;
    MAT2 LHS2;
    MAT3X2 LHS;
    LHS.col(0)=D;
    LHS.col(1)=-td;
    LHS2=LHS.transpose()*LHS;
    RHS2=LHS.transpose()*(tp0-p0);
    //invert
    T detLHS2=LHS2.determinant();
    if(detLHS2==0)
      t[off]=2;
    else {
      std::swap(LHS2(0,0),LHS2(1,1));
      LHS2(0,1)=-LHS2(0,1);
      LHS2(1,0)=-LHS2(1,0);
      LHS2/=detLHS2;
      t[off]=(LHS2*RHS2)[0];
    }
  }
  ASSERT(t[0]>0 || t[1]>0)
  if(t[0]<0) {
    if(t[1]<1)
      return std::make_pair((vid+2)%3,p0+D*t[1]);
    else return std::make_pair(-1,p0+D);
  } else if(t[1]<0) {
    if(t[0]<1)
      return std::make_pair((vid+1)%3,p0+D*t[0]);
    else return std::make_pair(-1,p0+D);
  } else if(t[0]>1 && t[1]>1)
    return std::make_pair(-1,p0+D);
  else if(t[0]<t[1])
    return std::make_pair((vid+1)%3,p0+D*t[0]);
  else if(t[1]<t[0])
    return std::make_pair((vid+2)%3,p0+D*t[1]);
  else
    return std::make_pair(vid,p0+D*t[0]);
}
void TriangleExact::calcPointDist(const PT& pt,T& sqrDistance,PT& cp,PT& b,Vec2i& feat) const
{
  feat=Vec2i(-1,-1);
  bool outside=false;
  for(int d=0; d<3; d++) {
    if((pt-_v[(d+1)%3]).dot(_n[d])>=0) {
      outside=true;
      T t=(pt-_v[(d+1)%3]).dot(_nO[d])/_nOSqr[d];
      if(t>0&&t<1) {
        //Edge's dirichlet region
        b.setZero();
        b[(d+1)%3]=1-t;
        b[(d+2)%3]=t;
        feat=Vec2i((d+1)%3,(d+2)%3);
        cp=b[0]*_v[0]+b[1]*_v[1]+b[2]*_v[2];
        sqrDistance=(pt-cp).squaredNorm();
        break;
      }
    }
  }
  if(!outside) {
    //Triangle's dirichlet region
    PT2 RHS;
    RHS[0]=(pt-_v[0]).dot(_nO[1]);
    RHS[1]=(pt-_v[0]).dot(_nO[2]);
    b.segment<2>(1)=_invM*RHS;
    b[0]=1-b[1]-b[2];
    feat=Vec2i(-1,-1);
    cp=b[0]*_v[0]+b[1]*_v[1]+b[2]*_v[2];
    sqrDistance=(pt-cp).squaredNorm();
  } else if(feat[0]==-1) {
    //Vertex's dirichlet region
    for(int d=0; d<3; d++) {
      T sqrDistanceNew=(pt-_v[d]).squaredNorm();
      if(d==0||sqrDistanceNew<sqrDistance) {
        b=PT::Unit(d);
        feat=Vec2i(d,-1);
        cp=_v[d];
        sqrDistance=sqrDistanceNew;
      }
    }
  }
}
void TriangleExact::writeVTK(const std::string& path,const PTss* pt,const PTss* cp) const
{
  VTKWriter<scalar> os("TriangleExact",path,true);
  std::vector<Vec3,Eigen::aligned_allocator<Vec3>> vss;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>> tss,pss,lss;
  vss.push_back(castRational<Vec3,PT>(_v[0]));
  vss.push_back(castRational<Vec3,PT>(_v[1]));
  vss.push_back(castRational<Vec3,PT>(_v[2]));
  tss.push_back(Vec3i(0,1,2));
  pss.push_back(Vec3i::Constant(0));
  pss.push_back(Vec3i::Constant(1));
  pss.push_back(Vec3i::Constant(2));
  lss.push_back(Vec3i(0,1,-1));
  lss.push_back(Vec3i(1,2,-1));
  lss.push_back(Vec3i(2,0,-1));

  if(pt && cp) {
    for(sizeType i=0; i<(sizeType)pt->size(); i++) {
      vss.push_back(castRational<Vec3,PT>(pt->at(i)));
      vss.push_back(castRational<Vec3,PT>(cp->at(i)));
      pss.push_back(Vec3i::Constant(3+i*2));
      pss.push_back(Vec3i::Constant(4+i*2));
      lss.push_back(Vec3i(3+i*2,4+i*2,-1));
    }
  }
  os.appendPoints(vss.begin(),vss.end());
  os.appendCells(tss.begin(),tss.end(),VTKWriter<scalar>::TRIANGLE);
  os.appendCells(pss.begin(),pss.end(),VTKWriter<scalar>::POINT);
  os.appendCells(lss.begin(),lss.end(),VTKWriter<scalar>::LINE);
}
void TriangleExact::debugPointDist(const std::string& path,sizeType nrIter)
{
#define DIST 10
  recreate(path);
  for(sizeType i=0; i<nrIter; i++) {
    TriangleTpl<scalar> t(Vec3::Random(),Vec3::Random(),Vec3::Random());
    TriangleExact et(t);

    T sqrDist;
    Vec2i feat;
    PTss pt,cp;
    PT l0=Vec3::Random().cast<double>().cast<T>()*DIST;
    PT l1=Vec3::Random().cast<double>().cast<T>()*DIST,b;
    for(sizeType i=0; i<DIST*10; i++) {
      T alpha=T(i)/T(DIST*10);
      pt.push_back(l0*(1-alpha)+l1*alpha);
      cp.push_back(PT());
      et.calcPointDist(pt.back(),sqrDist,cp.back(),b,feat);
    }
    et.writeVTK(path+"/frm"+std::to_string(i)+".vtk",&pt,&cp);
  }
}
void TriangleExact::reset()
{
  for(int d=0; d<3; d++) {
    _nO[d]=_v[(d+2)%3]-_v[(d+1)%3];

    _nOSqr[d]=_nO[d].dot(_nO[d]);

    T t=(_v[d]-_v[(d+1)%3]).dot(_nO[d])/_nOSqr[d];

    _n[d]=_v[(d+1)%3]+_nO[d]*t-_v[d];

  }
  _nt=(_v[1]-_v[0]).cross(_v[2]-_v[0]);

  _nTnN=(_nt*_nt.transpose())/_nt.squaredNorm();

  MAT2 M;
  M(0,0)= _nO[1].dot(_nO[2]);
  M(1,0)= _nO[2].dot(_nO[2]);
  M(0,1)=-_nO[1].dot(_nO[1]);
  M(1,1)=-_nO[2].dot(_nO[1]);

  _invM=M.inverse();

}
