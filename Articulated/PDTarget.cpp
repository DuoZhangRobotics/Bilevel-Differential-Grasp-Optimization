#include "PDTarget.h"
#include "PBDArticulatedGradientInfo.h"
#include <Utils/Utils.h>

USE_PRJ_NAMESPACE

PDTarget::PDTarget() {}
PDTarget::PDTarget(const Vec& PCoef,const Vec& DCoef,const std::vector<std::tuple<scalarD,Vec,Vec>>& spline)
  :_PCoef(PCoef),_DCoef(DCoef),_spline(spline),_id(0),_t(0)
{
  ASSERT(PCoef.size()==std::get<1>(spline[0]).size())
  ASSERT(DCoef.size()==std::get<1>(spline[0]).size())
  reset();
}
PDTarget::PDTarget(const Vec& PCoef,const Vec& DCoef,const Vec& s)
  :_PCoef(PCoef),_DCoef(DCoef),_s(s),_t(0)
{
  ASSERT(PCoef.size()==s.size()/2)
  ASSERT(DCoef.size()==s.size()/2)
  reset();
}
bool PDTarget::read(std::istream& is,IOData* dat)
{
  readBinaryData(_PCoef,is,dat);
  readBinaryData(_DCoef,is,dat);
  readBinaryData(_spline,is,dat);
  readBinaryData(_id,is,dat);
  readBinaryData(_s,is,dat);
  readBinaryData(_t,is,dat);
  return is.good();
}
bool PDTarget::write(std::ostream& os,IOData* dat) const
{
  writeBinaryData(_PCoef,os,dat);
  writeBinaryData(_DCoef,os,dat);
  writeBinaryData(_spline,os,dat);
  writeBinaryData(_id,os,dat);
  writeBinaryData(_s,os,dat);
  writeBinaryData(_t,os,dat);
  return os.good();
}
std::shared_ptr<SerializableBase> PDTarget::copy() const
{
  return std::shared_ptr<SerializableBase>(new PDTarget());
}
std::string PDTarget::type() const
{
  return typeid(PDTarget).name();
}
void PDTarget::writeVTKSeq(const ArticulatedBody& body,const std::string& path,T dt)
{
  reset();
  recreate(path);
  sizeType nrF=(sizeType)std::ceil(std::get<0>(_spline.back())/dt);
  for(sizeType id=0; id<nrF; id++) {
    INFOV("Writing PDTarget %ld/%ld",id,nrF)
    PBDArticulatedGradientInfo<T> info(body,s().segment(0,body.nrDOF()));
    body.writeVTK(info._TM,path+"/frm"+std::to_string(id)+".vtk",Joint::MESH);
    advance(dt);
  }
  reset();
}
void PDTarget::reset() {
  _t=0;
  if(_spline.size()>0)
    advanceSpline();
}
void PDTarget::advance(T dt) {
  _t+=dt;
  if(_spline.size()>0)
    advanceSpline();
}
PDTarget::Vec PDTarget::s() const {
  if(_spline.size()>0)
    return splineInterp();
  else if(_s.size()>0)
    return _s;
  else return Vec::Zero(0);
}
void PDTarget::setS(const Vec& s) {
  _s=s;
}
PDTarget::Vec PDTarget::splineInterp() const {
  T a=_t-std::get<0>(_spline[_id]);
  a/=(std::get<0>(_spline[_id+1])-std::get<0>(_spline[_id]));
  T a2=a*a,a3=a2*a;
  const Vec& p0=std::get<1>(_spline[_id]);
  const Vec& m0=std::get<2>(_spline[_id]);
  const Vec& p1=std::get<1>(_spline[_id+1]);
  const Vec& m1=std::get<2>(_spline[_id+1]);
  Vec pos=p0*(2*a3-3*a2+1)+m0*(a3-2*a2+a)+p1*(-2*a3+3*a2)+m1*(a3-a2);
  Vec vel=p0*(6*a2-6*a)+m0*(3*a2-4*a+1)+p1*(-6*a2+6*a)+m1*(3*a2-2*a);
  return concat(pos,vel);
}
void PDTarget::advanceSpline() {
  T horizon=std::get<0>(_spline.back());
  while(_t>horizon)
    _t-=horizon;
  while(_id<(sizeType)_spline.size()-2 && _t>std::get<0>(_spline[_id]))
    _id++;
  while(_id>0 && _t<std::get<0>(_spline[_id]))
    _id--;
}
