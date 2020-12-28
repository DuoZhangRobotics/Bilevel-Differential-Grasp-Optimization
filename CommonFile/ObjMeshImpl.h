#ifdef OBJMESH

OBJMESH::Edge::Edge():_nor(0.0f,0.0f,0.0f) {}
bool OBJMESH::EdgeMap::LSS::operator()(const std::pair<int,int>& a,const std::pair<int,int>& b) const
{
  return (a.first < b.first) || (a.first == b.first && a.second < b.second);
}
int OBJMESH::PovTexture::nrMAT() const
{
  return (int)_mat.size();
}
const int& OBJMESH::PovTexture::FMAT(int i) const
{
  return _fmat[i];
}
const std::string& OBJMESH::PovTexture::MAT(int i) const
{
  return _mat[i];
}
int OBJMESH::PovTexture::nrUV() const
{
  return (int)_uv.size();
}
const Vec3i& OBJMESH::PovTexture::FUV(int i) const
{
  return _fuv[i];
}
const OBJMESH::PT2& OBJMESH::PovTexture::UV(int i) const
{
  return _uv[i];
}

//the OBJMESH class
OBJMESH::OBJMESH():_id(0),_pos(0.0f,0.0f,0.0f),_scale(1.0f),_ctrOff(0.0f,0.0f,0.0f)
{
  _trans=MAT3::Identity();
  _dim=3;
}
OBJMESH::OBJMESH(const int& id):_id(id),_pos(0.0f,0.0f,0.0f),_scale(1.0f),_ctrOff(0.0f,0.0f,0.0f)
{
  _trans=MAT3::Identity();
  _dim=3;
}
void OBJMESH::setDim(int dim)
{
  _dim=dim;
}
OBJMESH::~OBJMESH() {}
bool OBJMESH::read(std::istream &is,bool move,bool scale,int assertFaceType,int assertTriangle)
{
  _dim=3;
  _ctrOff=PT3(0.0f,0.0f,0.0f);

  std::vector<PT3,Eigen::aligned_allocator<PT3> > tmpVss;
  _vss.swap(tmpVss);

  std::vector<PT3,Eigen::aligned_allocator<PT3> > tmpNss;
  _nss.swap(tmpNss);

  std::vector<PT3,Eigen::aligned_allocator<PT3> > tmpFNss;
  _fnss.swap(tmpFNss);

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > tmpIss;
  _iss.swap(tmpIss);

  std::vector<int> tmpIssg;
  _issg.swap(tmpIssg);

  std::map<std::string,int> tmpGss;
  _gss.swap(tmpGss);

  std::map<int,std::string> tmpIgss;
  _igss.swap(tmpIgss);

  _tex=PovTexture();

  char buf[4096];
  std::string gName;

  PT3 v;
  Eigen::Matrix<int,3,1> a,b;
  Eigen::Matrix<int,3,1> at,bt;
  Eigen::Matrix<int,3,1> an,bn;
  int gIndex=-1;
  int currIndex=-1;

  while(!is.eof() && is.good()) {
    int faceType=-1;
    int isTriangle=1;
    int hasG=1;

    bool match=false;
    is.getline(buf,std::streamsize(4096),'\n');
    switch(buf[0]) {
    case 'g':
#ifdef _MSC_VER
      hasG=sscanf_s(buf,"g %s",buf,4096);
#else
      hasG=sscanf(buf,"g %s",buf);
#endif
      if(hasG < 1)
        buf[0]='\0';
      gName=buf;
      if(_gss.find(gName) != _gss.end()) {
        currIndex=_gss[gName];
      } else {
        gIndex++;
        currIndex=gIndex;
        _gss[gName]=currIndex;
        _igss[currIndex]=gName;
      }
      break;
    case 'v':
      switch(buf[1]) {
      case ' ':
        WriteObjVertex<T>::write(buf+2,v.x(),v.y(),v.z());
        _vss.push_back(v);
        break;
      case 'n':
        WriteObjVertex<T>::write(buf+3,v.x(),v.y(),v.z());
        _nss.push_back(v);
        break;
      case 't':
        WriteObjVertex<T>::write(buf+3,v.x(),v.y(),v.z());
        _tex._uv.push_back(v.segment<2>(0));
        break;
      }
      break;
    case 'f':
      //0
      if(!match)
#ifdef _MSC_VER
        switch(sscanf_s(buf,"f %d %d %d %d",&(a.x()),&(a.y()),&(a.z()),&(b.x())))
#else
        switch(sscanf(buf,"f %d %d %d %d",&(a.x()),&(a.y()),&(a.z()),&(b.x())))
#endif
        {
        case 3:
          if(_tex.nrUV() > 0) {
            WARNING("Skip face without texture!")
            match=true;
            break;
          }
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          match=true;
          faceType=0;
          break;
        case 4:
          if(_tex.nrUV() > 0) {
            WARNING("Skip face without texture!")
            match=true;
            break;
          }
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          match=true;
          faceType=0;
          isTriangle=0;
          break;
        }
      //1
      if(!match)
#ifdef _MSC_VER
        switch(sscanf_s(buf,"f %d/%d %d/%d %d/%d %d/%d",&(a.x()),&(at.x()), &(a.y()),&(at.y()), &(a.z()),&(at.z()), &(b.x()),&(bt.x())))
#else
        switch(sscanf(buf,"f %d/%d %d/%d %d/%d %d/%d",&(a.x()),&(at.x()), &(a.y()),&(at.y()), &(a.z()),&(at.z()), &(b.x()),&(bt.x())))
#endif
        {
        case 6:
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _tex._fuv.push_back(Vec3i(at.x(),at.y(),at.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          match=true;
          faceType=1;
          break;
        case 8:
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _tex._fuv.push_back(Vec3i(at.x(),at.y(),at.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
          _tex._fuv.push_back(Vec3i(at.x(),at.z(),bt.x())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          match=true;
          faceType=1;
          isTriangle=0;
          break;
        }
      //2
      if(!match)
#ifdef _MSC_VER
        switch(sscanf_s(buf,"f %d//%d %d//%d %d//%d %d//%d",&(a.x()),&(an.x()), &(a.y()),&(an.y()), &(a.z()),&(an.z()), &(b.x()),&(bn.x())))
#else
        switch(sscanf(buf,"f %d//%d %d//%d %d//%d %d//%d",&(a.x()),&(an.x()), &(a.y()),&(an.y()), &(a.z()),&(an.z()), &(b.x()),&(bn.x())))
#endif
        {
        case 6:
          if(_tex.nrUV() > 0) {
            WARNING("Skip face without texture!")
            match=true;
            break;
          }
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _fnss.push_back(_nss[an.x()-1]);
          _fnss.push_back(_nss[an.y()-1]);
          _fnss.push_back(_nss[an.z()-1]);
          match=true;
          faceType=2;
          break;
        case 8:
          if(_tex.nrUV() > 0) {
            WARNING("Skip face without texture!")
            match=true;
            break;
          }
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _fnss.push_back(_nss[an.x()-1]);
          _fnss.push_back(_nss[an.y()-1]);
          _fnss.push_back(_nss[an.z()-1]);
          _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _fnss.push_back(_nss[an.x()-1]);
          _fnss.push_back(_nss[an.z()-1]);
          _fnss.push_back(_nss[bn.x()-1]);
          match=true;
          faceType=2;
          isTriangle=0;
          break;
        }
      //3
      if(!match)
#ifdef _MSC_VER
        switch(sscanf_s(buf,"f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",&(a.x()),&(at.x()),&(an.x()), &(a.y()),&(at.y()),&(an.y()), &(a.z()),&(at.z()),&(an.z()), &(b.x()),&(bt.x()),&(bn.x())))
#else
        switch(sscanf(buf,"f %d/%d/%d %d/%d/%d %d/%d/%d %d/%d/%d",&(a.x()),&(at.x()),&(an.x()), &(a.y()),&(at.y()),&(an.y()), &(a.z()),&(at.z()),&(an.z()), &(b.x()),&(bt.x()),&(bn.x())))
#endif
        {
        case 9:
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _tex._fuv.push_back(Vec3i(at.x(),at.y(),at.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _fnss.push_back(_nss[an.x()-1]);
          _fnss.push_back(_nss[an.y()-1]);
          _fnss.push_back(_nss[an.z()-1]);
          match=true;
          faceType=3;
          break;
        case 12:
          _iss.push_back(Vec3i(a.x(),a.y(),a.z())-Vec3i(1,1,1));
          _tex._fuv.push_back(Vec3i(at.x(),at.y(),at.z())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _fnss.push_back(_nss[an.x()-1]);
          _fnss.push_back(_nss[an.y()-1]);
          _fnss.push_back(_nss[an.z()-1]);
          _iss.push_back(Vec3i(a.x(),a.z(),b.x())-Vec3i(1,1,1));
          _tex._fuv.push_back(Vec3i(at.x(),at.z(),bt.x())-Vec3i(1,1,1));
          _issg.push_back(currIndex);
          _fnss.push_back(_nss[an.x()-1]);
          _fnss.push_back(_nss[an.z()-1]);
          _fnss.push_back(_nss[bn.x()-1]);
          match=true;
          faceType=3;
          isTriangle=0;
          break;
        }
      //4
      if(!match)
        return false;

      if(assertFaceType != -1 && assertFaceType != faceType)
        return false;
      if(assertTriangle != -1 && isTriangle != assertTriangle)
        return false;
    }
  }

  if(_tex.nrUV() > 0) {
    ASSERT_MSG(_tex._fuv.size() == _iss.size(),"Face UV size mismatch!")
  }

  if(_vss.empty() || _iss.empty()) {
    WARNING("No Vertex In Mesh")
    applyTrans();
    return true;
  }

  for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
    const Vec3i& tri=_iss[i];
    if(tri.x() < (int)0 || tri.x() >= (int)_vss.size()) {
      WARNING("Triangle Index Out Of Bounds")
      return false;
    }
    if(tri.y() < (int)0 || tri.y() >= (int)_vss.size()) {
      WARNING("Triangle Index Out Of Bounds")
      return false;
    }
    if(tri.z() < (int)0 || tri.z() >= (int)_vss.size()) {
      WARNING("Triangle Index Out Of Bounds")
      return false;
    }
  }

  //BBox for vssInput
  BBox<T> box;
  for(int i=0; i<(int)_vss.size(); i++)
    box.setUnion(_vss[i]);

  //if(!(compG(box._minC,PT3(-1E3f,-1E3f,-1E3f)) && compL(box._maxC,PT3(1E3,1E3,1E3))))
  //{
  //	WARNING("Too Big A Mesh")
  //	return false;
  //}

  _trans=MAT3::Identity();
  if(scale) {
    PT3 extent=box.getExtent();
    _scale=1.0f/std::min(extent.x(),std::min(extent.y(),extent.z()));
    for(int i=0,sz=(int)_vss.size(); i<sz; i++)
      _vss[i]*=_scale;
  }
  _scale=1.0f;

  _pos=PT3(0.0f,0.0f,0.0f);
  if(move) {
    BBox<T> bb;
    for(int i=0; i<(int)_vss.size(); i++)
      bb.setUnion(_vss[i]);

    for(int i=0,sz=(int)_vss.size(); i<sz; i++)
      _vss[i]-=bb._minC;
  }

  if(scale || move)
    applyTrans();
  else if(_nss.size() != _vss.size())
    smooth();
  return true;
}
bool OBJMESH::readBinary(std::istream& is)
{
  if(!is.good())
    return false;
  if(!readBinaryData(_dim,is).good())
    return false;

  std::vector<PT3,Eigen::aligned_allocator<PT3> >& vss=getV();
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss=getI();
  std::vector<int>& issg=getIG();
  std::vector<PT3,Eigen::aligned_allocator<PT3> >& nss=getN();
  std::vector<PT3,Eigen::aligned_allocator<PT3> >& fnss=getFN();
  std::vector<PT3,Eigen::aligned_allocator<PT3> >& tnss=getTN();

  int szV;
  int szT;
  int szTG;
  int szG;

  if(!readBinaryData(szV,is).good())
    return false;
  if(!readBinaryData(szT,is).good())
    return false;
  if(!readBinaryData(szTG,is).good())
    return false;
  if(!readBinaryData(szG,is).good())
    return false;

  vss.resize(szV);
  nss.resize(szV);
  iss.resize(szT);
  issg.resize(szTG);

  if(!readBinaryData(vss,is).good())
    return false;
  if(!readBinaryData(nss,is).good())
    return false;
  if(!readBinaryData(fnss,is).good())
    return false;
  if(!readBinaryData(tnss,is).good())
    return false;
  if(!readBinaryData(iss,is).good())
    return false;
  if(!readBinaryData(issg,is).good())
    return false;

  for(int ig=0; ig<szG; ig++) {
    std::string name;
    int id;

    readBinaryData(name,is);
    readBinaryData(id,is);
    _gss[name]=id;
    _igss[id]=name;

    if(!is.good())
      return false;
  }

  if(!readBinaryData(_ctrOff,is).good())
    return false;

  return true;
}
bool OBJMESH::write(std::ostream &os) const
{
  if(!os.good())
    return false;
  if(_dim == 2)
    return false;

  char buf[4096];
  for(int i=0,sz=(int)_vss.size(); i<sz; i++) {
    const PT3D vd=_vss[i].cast<double>();
#ifdef _MSC_VER
    sprintf_s(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#else
    sprintf(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#endif
    os << buf;
    if(!os.good())
      return false;
  }

  if(_nss.size() == _vss.size()) {
    for(int i=0,sz=(int)_nss.size(); i<sz; i++) {
      const PT3D vn=_nss[i].cast<double>();
#ifdef _MSC_VER
      sprintf_s(buf,"vn %f %f %f\n",vn.x(),vn.y(),vn.z());
#else
      sprintf(buf,"vn %f %f %f\n",vn.x(),vn.y(),vn.z());
#endif
      os << buf;
      if(!os.good())
        return false;
    }
  }

  if(_issg.empty()) {
    for(int i=0,sz=(int)_iss.size(); i<sz; i++) {
#ifdef _MSC_VER
      if(_nss.size() == _vss.size())
        sprintf_s(buf,"f %I64d//%I64d %I64d//%I64d %I64d//%I64d\n",
                  (int64_t)_iss[i].x()+1,(int64_t)_iss[i].x()+1,
                  (int64_t)_iss[i].y()+1,(int64_t)_iss[i].y()+1,
                  (int64_t)_iss[i].z()+1,(int64_t)_iss[i].z()+1);
      else
        sprintf_s(buf,"f %I64d %I64d %I64d\n",
                  (int64_t)_iss[i].x()+1,
                  (int64_t)_iss[i].y()+1,
                  (int64_t)_iss[i].z()+1);
#else
      if(_nss.size() == _vss.size())
        sprintf(buf,"f %ld//%ld %ld//%ld %ld//%ld\n",
                (int64_t)_iss[i].x()+1,(int64_t)_iss[i].x()+1,
                (int64_t)_iss[i].y()+1,(int64_t)_iss[i].y()+1,
                (int64_t)_iss[i].z()+1,(int64_t)_iss[i].z()+1);
      else
        sprintf(buf,"f %ld %ld %ld\n",
                (int64_t)_iss[i].x()+1,
                (int64_t)_iss[i].y()+1,
                (int64_t)_iss[i].z()+1);
#endif
      os << buf;
    }
  } else {
    std::map<int,std::string> groupOs;
    for(int i=0,sz=(int)_iss.size(); i<sz; i++) {
      if(groupOs.find(_issg[i]) == groupOs.end()) {
        groupOs[_issg[i]]=std::string();
        if(_issg[i] != -1) {
          std::string currGroupName=_igss.find(_issg[i])->second;
          std::string& str=groupOs[_issg[i]];
          str+="g ";
          str+=currGroupName;
          str+="\n";
        }
      }

#ifdef _MSC_VER
      if(_nss.size() == _vss.size())
        sprintf_s(buf,"f %I64d//%I64d %I64d//%I64d %I64d//%I64d\n",
                  (int64_t)_iss[i].x()+1,(int64_t)_iss[i].x()+1,
                  (int64_t)_iss[i].y()+1,(int64_t)_iss[i].y()+1,
                  (int64_t)_iss[i].z()+1,(int64_t)_iss[i].z()+1);
      else
        sprintf_s(buf,"f %I64d %I64d %I64d\n",
                  (int64_t)_iss[i].x()+1,
                  (int64_t)_iss[i].y()+1,
                  (int64_t)_iss[i].z()+1);
#else
      if(_nss.size() == _vss.size())
        sprintf(buf,"f %ld//%ld %ld//%ld %ld//%ld\n",
                (int64_t)_iss[i].x()+1,(int64_t)_iss[i].x()+1,
                (int64_t)_iss[i].y()+1,(int64_t)_iss[i].y()+1,
                (int64_t)_iss[i].z()+1,(int64_t)_iss[i].z()+1);
      else
        sprintf(buf,"f %ld %ld %ld\n",
                (int64_t)_iss[i].x()+1,
                (int64_t)_iss[i].y()+1,
                (int64_t)_iss[i].z()+1);
#endif
      groupOs[_issg[i]]+=buf;
    }

    std::map<int,std::string>::iterator
    begin=groupOs.begin(),end=groupOs.end();
    for(; begin!=end; begin++)
      os << begin->second;
  }

  return os.good();
}
bool OBJMESH::write(const std::string& path) const
{
  std::ofstream ofs(path);
  return write(ofs);
}
bool OBJMESH::writePov(const std::string& path,bool normal,const PovTexture* tex) const
{
  std::ofstream os(path);
  return writePov(os,normal,tex);
}
bool OBJMESH::writePov(std::ostream &os,bool normal,const PovTexture* tex) const
{
  if(_dim == 2)
    return false;
#define CACHE_SIZE 1024

  os << "mesh2 {\n";

  //vertices
  {
    os << "    vertex_vectors {" << _vss.size() << "\n";
    int nr=(int)_vss.size();
    std::ostringstream oss;

    for(int i=0; i<(int)_vss.size(); i++) {
      char buf[CACHE_SIZE];
      const PT3D v=_vss[i].cast<double>();

#ifdef _MSC_VER
      if(i == nr-1)
        sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>\n",v.x(),v.y(),v.z());
      else
        sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>,\n",v.x(),v.y(),v.z());
#else
      if(i == nr-1)
        sprintf(buf,"        <%f,%f,%f>\n",v.x(),v.y(),v.z());
      else
        sprintf(buf,"        <%f,%f,%f>,\n",v.x(),v.y(),v.z());
#endif
      oss << buf;
    }
    os << oss.str() << "    }\n";
  }

  //normals
  if(normal && _nss.size() == _vss.size()) {
    os << "    normal_vectors {" << _nss.size() << "\n";
    int nr=(int)_nss.size();
    std::ostringstream oss;

    for(int i=0; i<nr; i++) {
      char buf[CACHE_SIZE];
      const PT3D n=_nss[i].cast<double>();
#ifdef _MSC_VER
      if(i == nr-1)
        sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>\n",n.x(),n.y(),n.z());
      else
        sprintf_s<CACHE_SIZE>(buf,"        <%f,%f,%f>,\n",n.x(),n.y(),n.z());
#else
      if(i == nr-1)
        sprintf(buf,"        <%f,%f,%f>\n",n.x(),n.y(),n.z());
      else
        sprintf(buf,"        <%f,%f,%f>,\n",n.x(),n.y(),n.z());
#endif
      oss << buf;
    }

    os << oss.str() << "    }\n";
  }

  //uvs
  if(tex && tex->nrUV() > 0) {
    os << "    uv_vectors {" << tex->nrUV() << "\n";
    int nr=tex->nrUV();
    std::ostringstream oss;

    for(int i=0; i<nr; i++) {
      char buf[CACHE_SIZE];
      const PT2D n=tex->UV(i).cast<double>();
#ifdef _MSC_VER
      if(i == nr-1)
        sprintf_s<CACHE_SIZE>(buf,"        <%f,%f>\n",n.x(),n.y());
      else
        sprintf_s<CACHE_SIZE>(buf,"        <%f,%f>,\n",n.x(),n.y());
#else
      if(i == nr-1)
        sprintf(buf,"        <%f,%f>\n",n.x(),n.y());
      else
        sprintf(buf,"        <%f,%f>,\n",n.x(),n.y());
#endif
      oss << buf;
    }

    os << oss.str() << "    }\n";
  }

  //face indices
  if(tex && tex->nrMAT() > 0) {
    int nrT=tex->nrMAT();
    os << "    texture_list {" << nrT << "\n";
    for(int i=0; i<nrT; i++) {
      os << tex->MAT(i);
      if(i<nrT-1)
        os << ",";
      os << "\n";
    }
    os << "    }\n";
  }
  {
    os << "    face_indices {" << _iss.size() << "\n";
    int nr=(int)_iss.size();
    std::ostringstream oss;

    for(int i=0; i<nr; i++) {
      char buf[CACHE_SIZE];
      const Vec3i& f=_iss[i];
      if(tex && tex->nrMAT() > 0) {
#ifdef _MSC_VER
        if(i == nr-1)
          sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>%d\n" ,(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z(),tex->FMAT(i));
        else
          sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>%d,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z(),tex->FMAT(i));
#else
        if(i == nr-1)
          sprintf(buf,"        <%ld,%ld,%ld>%d\n", (int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z(),tex->FMAT(i));
        else
          sprintf(buf,"        <%ld,%ld,%ld>%d,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z(),tex->FMAT(i));
#endif
      } else {
#ifdef _MSC_VER
        if(i == nr-1)
          sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>\n" ,(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
        else
          sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
#else
        if(i == nr-1)
          sprintf(buf,"        <%ld,%ld,%ld>\n" ,(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
        else
          sprintf(buf,"        <%ld,%ld,%ld>,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
#endif
      }
      oss << buf;
    }
    os << oss.str() << "    }\n";
  }

  if(normal && _nss.size() == _vss.size()) {
    os << "    normal_indices {" << _iss.size() << "\n";
    int nr=(int)_iss.size();
    std::ostringstream oss;

    for(int i=0; i<nr; i++) {
      char buf[CACHE_SIZE];
      const Vec3i& f=_iss[i];
#ifdef _MSC_VER
      if(i == nr-1)
        sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>\n" ,(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
      else
        sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
#else
      if(i == nr-1)
        sprintf(buf,"        <%ld,%ld,%ld>\n" ,(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
      else
        sprintf(buf,"        <%ld,%ld,%ld>,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
#endif
      oss << buf;
    }

    os << oss.str() << "    }\n";
  }

  if(tex && tex->nrUV() > 0) {
    os << "    uv_indices {" << _iss.size() << "\n";
    int nr=(int)_iss.size();
    std::ostringstream oss;

    for(int i=0; i<nr; i++) {
      char buf[CACHE_SIZE];
      const Vec3i& f=tex->FUV(i);
      ASSERT(f[0] < tex->nrUV() && f[1] < tex->nrUV() && f[2] < tex->nrUV())
#ifdef _MSC_VER
      if(i == nr-1)
        sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>\n" ,(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
      else
        sprintf_s<CACHE_SIZE>(buf,"        <%I64d,%I64d,%I64d>,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
#else
      if(i == nr-1)
        sprintf(buf,"        <%ld,%ld,%ld>\n" ,(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
      else
        sprintf(buf,"        <%ld,%ld,%ld>,\n",(int64_t)f.x(),(int64_t)f.y(),(int64_t)f.z());
#endif
      oss << buf;
    }

    os << oss.str() << "    }\n";
  }
  os << "}\n";
  return os.good();

#undef CACHE_SIZE
}
bool OBJMESH::writePBRT(std::ostream &os) const
{
  if(!os.good())
    return false;
  if(_dim == 2)
    return false;

  os << "AttributeBegin" << std::endl;
  os << "Shape \"trianglemesh\"" << std::endl;
  if(!os.good())
    return false;

  os << "\"point P\" ["<< std::endl;
  for(int j=0; j<(int)_vss.size(); j++)
    os << _vss[j].x() << " " << _vss[j].y() << " " << _vss[j].z() << std::endl;
  os << "]" << std::endl;
  if(!os.good())
    return false;

  os << "\"integer indices\" [" << std::endl;
  for(int j=0; j<(int)_iss.size(); j++)
    os << _iss[j].x() << " " << _iss[j].y() << " " << _iss[j].z() << std::endl;
  os << "]" << std::endl;
  os << "AttributeEnd" << std::endl;
  if(!os.good())
    return false;

  return true;
}
bool OBJMESH::writeVTK(const std::string& path,bool binary,bool normal,bool vertexNormal,const PT4* color,const std::vector<T>* cellColor) const
{
  VTKWriter<T> os("WaveFront Obj Mesh",path,binary);
  if(_dim == 3)
    return writeVTK3D(os,normal,vertexNormal,color,cellColor);
  else return writeVTK2D(os,normal,vertexNormal,color,cellColor);
}
bool OBJMESH::writeVTK(VTKWriter<T>& os,bool normal,bool vertexNormal,const PT4* color,const std::vector<T>* cellColor) const
{
  if(_dim == 3)
    return writeVTK3D(os,normal,vertexNormal,color,cellColor);
  else return writeVTK2D(os,normal,vertexNormal,color,cellColor);
}
bool OBJMESH::writeVTK3D(VTKWriter<T>& os,bool normal,bool vertexNormal,const PT4* color,const std::vector<T>* cellColor) const
{
  {
    os.setRelativeIndex();
    os.appendPoints(getV().begin(),getV().end());
    os.appendCells(getI().begin(),getI().end(),VTKWriter<T>::TRIANGLE,true);
    if(color) {
      VTKWriter<T>::IteratorRepeat<PT4> beg(0,*color),end((sizeType)_vss.size(),*color);
      os.appendCustomPointColorData("color",beg,end);
    }
    if(cellColor) {
      ASSERT_MSG(cellColor->size() == getI().size(),"CellColor size error!")
      os.appendCustomData("cellColor",cellColor->begin(),cellColor->end());
    }
  }
  const T lNor=estimateNormalLength();
  if(_tnss.size() == _iss.size() && normal) {
    os.setRelativeIndex();
    std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
    for(int i=0; i<(int)_tnss.size(); i++) {
      PT3 ctr=getTC((int)i);
      vss.push_back(ctr);
      vss.push_back(ctr+_tnss[i]*lNor);
    }
    os.appendPoints(vss.begin(),vss.end());
    // original is
    /*
      VTKWriter<T>::IteratorIndex<Vec3i> beg(0,2,0),end(vss.size()/2,2,0);
     */
    VTKWriter<T>::IteratorIndex<Vec3i> beg(0,2,0);
    VTKWriter<T>::IteratorIndex<Vec3i> end(vss.size()/2,2,0);
    os.appendCells(beg,end,VTKWriter<T>::LINE,true);
  }
  if(_nss.size() == _vss.size() && vertexNormal) {
    os.setRelativeIndex();
    std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
    for(int i=0; i<(int)_vss.size(); i++) {
      vss.push_back(_vss[i]);
      vss.push_back(_vss[i]+_nss[i]*lNor);
    }
    os.appendPoints(vss.begin(),vss.end());
    VTKWriter<T>::IteratorIndex<Vec3i> beg(0,2,0);
    VTKWriter<T>::IteratorIndex<Vec3i> end(vss.size()/2,2,0);
    os.appendCells(beg,end,VTKWriter<T>::LINE,true);
  }
  return true;
}
bool OBJMESH::writeVTK2D(VTKWriter<T>& os,bool normal,bool vertexNormal,const PT4* color,const std::vector<T>* cellColor) const
{
  {
    os.setRelativeIndex();
    os.appendPoints(getV().begin(),getV().end());
    os.appendCells(getI().begin(),getI().end(),VTKWriter<T>::LINE,true);
    if(color) {
      VTKWriter<T>::IteratorRepeat<PT4> beg(0,*color),end((sizeType)_vss.size(),*color);
      os.appendCustomPointColorData("color",beg,end);
    }
    if(cellColor) {
      ASSERT_MSG(cellColor->size() == getI().size(),"CellColor size error!")
      os.appendCustomData("cellColor",cellColor->begin(),cellColor->end());
    }
  }
  const T lNor=estimateNormalLength();
  if(_tnss.size() == _iss.size() && normal) {
    os.setRelativeIndex();
    std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
    for(int i=0; i<(int)_iss.size(); i++) {
      Vec3i LI=_iss[i];
      PT3 ctr=(_vss[LI[0]]+_vss[LI[1]])/2.0f;
      vss.push_back(ctr);
      vss.push_back(ctr+_tnss[i]*lNor);
    }
    os.appendPoints(vss.begin(),vss.end());
    VTKWriter<T>::IteratorIndex<Vec3i> beg(0,2,0),end(vss.size()/2,2,0);
    os.appendCells(beg,end,VTKWriter<T>::LINE,true);
  }
  if(_nss.size() == _vss.size() && vertexNormal) {
    os.setRelativeIndex();
    std::vector<PT3,Eigen::aligned_allocator<PT3> > vss;
    for(int i=0; i<(int)_vss.size(); i++) {
      vss.push_back(_vss[i]);
      vss.push_back(_vss[i]+_nss[i]*lNor);
    }
    os.appendPoints(vss.begin(),vss.end());
    VTKWriter<T>::IteratorIndex<Vec3i> beg(0,2,0);
    VTKWriter<T>::IteratorIndex<Vec3i> end(vss.size()/2,2,0);
    os.appendCells(beg,end,VTKWriter<T>::LINE,true);
  }
  return true;
}
bool OBJMESH::writeBinary(std::ostream &os) const
{
  if(!os.good())
    return false;
  if(!writeBinaryData(_dim,os).good())
    return false;

  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& vss=getV();
  const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& iss=getI();
  const std::vector<int>& issg=getIG();
  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& nss=getN();
  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& fnss=getFN();
  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& tnss=getTN();
  const std::map<std::string,int>& gss=_gss;

  const int szV=(int)vss.size();
  const int szT=(int)iss.size();
  const int szTG=(int)issg.size();
  const int szG=(int)gss.size();

  if(!writeBinaryData(szV,os).good())
    return false;
  if(!writeBinaryData(szT,os).good())
    return false;
  if(!writeBinaryData(szTG,os).good())
    return false;
  if(!writeBinaryData(szG,os).good())
    return false;

  if(!writeBinaryData(vss,os).good())
    return false;
  if(!writeBinaryData(nss,os).good())
    return false;
  if(!writeBinaryData(fnss,os).good())
    return false;
  if(!writeBinaryData(tnss,os).good())
    return false;
  if(!writeBinaryData(iss,os).good())
    return false;
  if(!writeBinaryData(issg,os).good())
    return false;
  for(std::map<std::string,int>::const_iterator
      beg=gss.begin(),end=gss.end(); beg != end; beg++) {
    writeBinaryData(beg->first,os);
    writeBinaryData(beg->second,os);
    if(!os.good())
      return false;
  }

  if(!writeBinaryData(_ctrOff,os).good())
    return false;

  return true;
}
bool OBJMESH::write(std::ostream &vs,std::ostream &fs,int &index) const
{
  if(!vs.good() || !fs.good())
    return false;
  if(_dim == 2)
    return false;

  int start=index;
  WARNING("Group Info In OBJMESH Is Ignored!")

  char buf[4096];
  for(int i=0,sz=(int)_vss.size(); i<sz; i++) {
    const PT3D vd=_vss[i].cast<double>();
#ifdef _MSC_VER
    sprintf_s(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#else
    sprintf(buf,"v %f %f %f\n",vd.x(),vd.y(),vd.z());
#endif
    index++;
    vs << buf;
    if(!vs.good())
      return false;
  }

  for(int i=0,sz=(int)_iss.size(); i<sz; i++) {
#ifdef _MSC_VER
    sprintf_s(buf,"f %I64d %I64d %I64d\n", (int64_t)start+_iss[i].x(), (int64_t)start+_iss[i].y(), (int64_t)start+_iss[i].z());
#else
    sprintf(buf,  "f %ld %ld %ld\n"      , (int64_t)start+_iss[i].x(), (int64_t)start+_iss[i].y(), (int64_t)start+_iss[i].z());
#endif
    fs << buf;
    if(!fs.good())
      return false;
  }

  return true;
}
void OBJMESH::writeCsv(const std::string& path) const
{
#define WRITE_POINT(ID)	os << _vss[ID][0] << "," << _vss[ID][1] << "," << _vss[ID][2] << "," << std::endl;
  std::ofstream os(path);
  int nrT=(int)_iss.size();
  os << nrT << "," << std::endl;
  for(int i=0; i<nrT; i++) {
    WRITE_POINT(_iss[i][0])
    WRITE_POINT(_iss[i][1])
    if(getDim() == 3)
      WRITE_POINT(_iss[i][2])
    }
}
void OBJMESH::addMesh(const OBJMESH& mesh,const std::string& g)
{
  int vOff=(int)_vss.size();
  int nrg=(int)_gss.size();
  _vss.insert(_vss.end(),mesh._vss.begin(),mesh._vss.end());
  for(int i=0; i<(int)mesh._iss.size(); i++) {
    _iss.push_back(mesh._iss[i]+Vec3i::Constant(vOff));
    _issg.push_back((int)nrg);
  }
  ASSERT(_gss.find(g) == _gss.end())
  _gss[g]=nrg;
  _igss[nrg]=g;
}
void OBJMESH::addMesh(const OBJMESH& mesh)
{
  int vOff=(int)_vss.size();
  _vss.insert(_vss.end(),mesh._vss.begin(),mesh._vss.end());
  for(int i=0; i<(int)mesh._iss.size(); i++)
    _iss.push_back(mesh._iss[i]+Vec3i::Constant(vOff));
}
OBJMESH::PT3 OBJMESH::getTC(int i) const
{
  if(_dim == 3)
    return (_vss[_iss[i].x()]+_vss[_iss[i].y()]+_vss[_iss[i].z()])/3.0f;
  else return (_vss[_iss[i].x()]+_vss[_iss[i].y()])/2.0f;
}
const OBJMESH::PT3& OBJMESH::getV(int i) const
{
  return _vss[i];
}
const Vec3i& OBJMESH::getI(int i) const
{
  return _iss[i];
}
const int& OBJMESH::getIG(int i) const
{
  return _issg[i];
}
const OBJMESH::PT3& OBJMESH::getN(int i) const
{
  return _nss[i];
}
const OBJMESH::PT3& OBJMESH::getTN(int i) const
{
  return _tnss[i];
}
std::map<std::string,int>& OBJMESH::getGS()
{
  return _gss;
}
std::map<int,std::string>& OBJMESH::getIGS()
{
  return _igss;
}
OBJMESH::PT3& OBJMESH::getV(int i)
{
  return _vss[i];
}
Vec3i& OBJMESH::getI(int i)
{
  return _iss[i];
}
int& OBJMESH::getIG(int i)
{
  return _issg[i];
}
OBJMESH::PT3& OBJMESH::getN(int i)
{
  return _nss[i];
}
OBJMESH::PT3& OBJMESH::getTN(int i)
{
  return _tnss[i];
}
OBJMESH::T OBJMESH::getArea(int i) const
{
  if(_dim == 3)
    return TriangleTpl<T>(_vss[_iss[i][0]],_vss[_iss[i][1]],_vss[_iss[i][2]]).area();
  else return LineSegTpl<T>(_vss[_iss[i][0]],_vss[_iss[i][1]]).length();
}
const std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getV() const
{
  return _vss;
}
const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& OBJMESH::getI() const
{
  return _iss;
}
const std::vector<int>& OBJMESH::getIG() const
{
  return _issg;
}
const std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getN() const
{
  return _nss;
}
const std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getFN() const
{
  return _fnss;
}
const std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getTN() const
{
  return _tnss;
}
const std::map<std::string,int>& OBJMESH::getGS() const
{
  return _gss;
}
const std::map<int,std::string>& OBJMESH::getIGS() const
{
  return _igss;
}
std::string OBJMESH::getTG(int i) const
{
  if(_issg[i] == -1)
    return std::string();
  else return _igss.find((int)_issg[i])->second;
}
int OBJMESH::getGId(const std::string& name) const
{
  if(_gss.find(name) == _gss.end())
    return -1;
  else return _gss.find(name)->second;
}
std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getV()
{
  return _vss;
}
std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& OBJMESH::getI()
{
  return _iss;
}
std::vector<int>& OBJMESH::getIG()
{
  return _issg;
}
std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getN()
{
  return _nss;
}
std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getFN()
{
  return _fnss;
}
std::vector<OBJMESH::PT3,Eigen::aligned_allocator<OBJMESH::PT3> >& OBJMESH::getTN()
{
  return _tnss;
}
const OBJMESH::MAT3& OBJMESH::getT() const
{
  return _trans;
}
OBJMESH::MAT3& OBJMESH::getT()
{
  return _trans;
}
const OBJMESH::PT3& OBJMESH::getPos() const
{
  return _pos;
}
OBJMESH::PT3& OBJMESH::getPos()
{
  return _pos;
}
const OBJMESH::T& OBJMESH::getScale() const
{
  return _scale;
}
OBJMESH::T& OBJMESH::getScale()
{
  return _scale;
}
void OBJMESH::applyTrans()
{
  BBox<T> bb;
  for(int i=0; i<(int)_vss.size(); i++)
    bb.setUnion(_vss[i]);
  PT3 ctr=(bb._minC+bb._maxC)*0.5f;
  applyTrans(ctr);
}
void OBJMESH::applyTrans(const PT3& customCtr)
{
  for(int i=0,sz=(int)_vss.size(); i<sz; i++)
    _vss[i]=_trans*(_vss[i]-customCtr)*_scale+customCtr+_pos;
  _ctrOff=_trans*(_ctrOff)*_scale;

  _trans=MAT3::Identity();
  _pos=PT3(0.0f,0.0f,0.0f);
  _scale=1.0f;

  smooth();
}
BBox<OBJMESH::T> OBJMESH::getBB() const
{
  BBox<T> ret;
  for(int i=0; i<(int)_vss.size(); i++)
    ret.setUnion(_vss[i]);
  return ret;
}
const int& OBJMESH::getId() const
{
  return _id;
}
//physics properties
OBJMESH::T OBJMESH::getVolume() const
{
  T volume = 0.0;
  for(int it=0,nr=(int)_iss.size(); it<nr; it++) {
    if(_dim == 3) {
      const PT3 &p1 = _vss[_iss[it].x()];
      const PT3 &p2 = _vss[_iss[it].y()];
      const PT3 &p3 = _vss[_iss[it].z()];
      volume += TriangleTpl<T>(p1, p2, p3).signedVolume();
    } else {
      const PT3 &p1 = _vss[_iss[it].x()];
      const PT3 &p2 = _vss[_iss[it].y()];
      volume += LineSegTpl<T>(p1, p2).signedArea();
    }
  }
  return volume;
}
OBJMESH::PT3 OBJMESH::getCentroid() const
{
  //must be convex
  PT3 numerator(0.0f,0.0f,0.0f);
  T denominator(0.0f);
  for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
    T area=getArea(i);
    PT3 ctr=getTC(i);
    numerator+=ctr*area;
    denominator+=area;
  }
  return _ctrOff+numerator/denominator;
}
OBJMESH::PT3 OBJMESH::getVolumeCentroid() const
{
  //must be convex
  PT3 numerator(0.0f,0.0f,0.0f);
  T denominator(0.0f);
  for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
    Vec3i TI=_iss[i];
    T vol;
    if(_dim == 3) {
      TetrahedronTpl<T> tet(PT3::Zero(),_vss[TI[0]],_vss[TI[1]],_vss[TI[2]]);
      vol=tet.volume();
      if(tet._swap)
        vol*=-1.0f;
      numerator+=tet.masscenter()*vol;
    } else {
      TriangleTpl<T> tri(PT3::Zero(),_vss[TI[0]],_vss[TI[1]]);
      vol=tri.area();
      numerator+=tri.masscenter()*vol;
    }
    denominator+=vol;
  }
  return _ctrOff+numerator/denominator;
}
OBJMESH::PT3& OBJMESH::centroidOffset()
{
  return _ctrOff;
}
const OBJMESH::PT3& OBJMESH::centroidOffset() const
{
  return _ctrOff;
}
//simple utility
OBJMESH::T OBJMESH::getMass(const T& dens) const
{
  T vol=getVolume();
  ASSERT_MSG(vol>0,"Negative Volume!")
  return vol*dens;
}
void OBJMESH::subdivide(const int& nrIter,std::vector<std::pair<int,int> >* vssInfo)
{
  EdgeMap eMap;
  buildEdge(eMap);
  for(std::map<std::pair<int,int>,Edge,EdgeMap::LSS>::iterator
      beg=eMap._ess.begin(),end=eMap._ess.end(); beg != end; beg++)
    beg->second._subdId=beg->first;
  if(vssInfo)
    vssInfo->assign(_vss.size(),std::pair<int,int>(-1,-1));
  for(int it=0,mod=1; it<nrIter; it++,mod*=4) {
    if(_dim == 3)
      subdivideSingle3D(eMap,mod,vssInfo);
    else subdivideSingle2D();
  }
  smooth();
}
void OBJMESH::subdivideSingle3D(EdgeMap& eMap,int mod,std::vector<std::pair<int,int> >* vssInfo)
{
  EdgeMap eMapOut;

  //first for each edge insert an internode
  for(std::map<std::pair<int,int>,Edge,EdgeMap::LSS>::iterator
      beg=eMap._ess.begin(),end=eMap._ess.end(); beg != end; beg++) {
    //insert node
    beg->second._interNode=(int)_vss.size();
    _vss.push_back((_vss[beg->first.first]+_vss[beg->first.second])*0.5f);
    if(vssInfo)
      vssInfo->push_back(beg->second._subdId);
    //parent info
    {
      std::pair<int,int> PA(beg->second._interNode,beg->first.first);
      std::pair<int,int> PB(beg->second._interNode,beg->first.second);
      if(PA.first > PA.second)std::swap(PA.first,PA.second);
      if(PB.first > PB.second)std::swap(PB.first,PB.second);
      Edge EA;
      EA._subdId=beg->second._subdId;
      Edge EB;
      EB._subdId=beg->second._subdId;
      eMapOut._ess[PA]=EA;
      eMapOut._ess[PB]=EB;
    }
  }

  //second replace each triangle with four
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > issNew;
  std::vector<int> issgNew;
  PovTexture texNew;
  for(int i=0,nr=(int)_iss.size(); i<nr; i++) {
    const Vec3i& iT=_iss[i];

    const int v0=(int)iT.x();
    const int v1=(int)iT.y();
    const int v2=(int)iT.z();
    const int v3=getE(v0,v1,eMap)._interNode;
    const int v4=getE(v1,v2,eMap)._interNode;
    const int v5=getE(v2,v0,eMap)._interNode;

    if(_tex.nrUV() > 0) {
      PT2 t0=_tex._uv[_tex._fuv[i][0]];
      PT2 t1=_tex._uv[_tex._fuv[i][1]];
      PT2 t2=_tex._uv[_tex._fuv[i][2]];
      PT2 t3=(t0+t1)/2;
      PT2 t4=(t1+t2)/2;
      PT2 t5=(t2+t0)/2;

      int off=(int)texNew._uv.size();
      texNew._uv.push_back(t0);
      texNew._uv.push_back(t1);
      texNew._uv.push_back(t2);
      texNew._uv.push_back(t3);
      texNew._uv.push_back(t4);
      texNew._uv.push_back(t5);
      texNew._fuv.push_back(Vec3i(off+0,off+3,off+5));
      texNew._fuv.push_back(Vec3i(off+3,off+1,off+4));
      texNew._fuv.push_back(Vec3i(off+5,off+3,off+4));
      texNew._fuv.push_back(Vec3i(off+5,off+4,off+2));
    }
    issNew.push_back(Vec3i(v0,v3,v5));
    issNew.push_back(Vec3i(v3,v1,v4));
    issNew.push_back(Vec3i(v5,v3,v4));
    issNew.push_back(Vec3i(v5,v4,v2));

    //make same group
    if(_issg.size() == _iss.size()) {
      const int& igT=_issg[i];
      issgNew.push_back(igT);
      issgNew.push_back(igT);
      issgNew.push_back(igT);
      issgNew.push_back(igT);
    }

    //parent info
    {
      std::pair<int,int> PA(v3,v4);
      std::pair<int,int> PB(v3,v5);
      std::pair<int,int> PC(v4,v5);
      if(PA.first > PA.second)std::swap(PA.first,PA.second);
      if(PB.first > PB.second)std::swap(PB.first,PB.second);
      if(PC.first > PC.second)std::swap(PC.first,PC.second);
      Edge EA;
      EA._subdId=std::pair<int,int>(i/mod,-1);
      Edge EB;
      EB._subdId=std::pair<int,int>(i/mod,-1);
      Edge EC;
      EC._subdId=std::pair<int,int>(i/mod,-1);
      eMapOut._ess[PA]=EA;
      eMapOut._ess[PB]=EB;
      eMapOut._ess[PC]=EC;
    }
  }
  eMap=eMapOut;

  _iss=issNew;
  _tex=texNew;
  if(issgNew.size() == _iss.size())
    _issg=issgNew;
}
void OBJMESH::subdivideSingle2D()
{
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > issNew;
  for(int i=0; i<(int)_iss.size(); i++) {
    Vec3i LI=_iss[i];
    int back=(int)_vss.size();

    issNew.push_back(Vec3i(LI[0],back,0));
    issNew.push_back(Vec3i(back,LI[1],0));

    _vss.push_back((_vss[LI[0]]+_vss[LI[1]])/2.0f);
  }

  _iss=issNew;
}
void OBJMESH::marchingCube(const std::vector<T>& MCVal,OBJMESH& mesh)
{
  //marching cube on surface
  mesh=OBJMESH();
  //calculate normal dot view dir
  std::vector<int> vid(_vss.size(),-1);
  for(int i=0; i<(int)_vss.size(); i++) {
    if(MCVal[i] < 0.0f) {
      vid[i]=(int)mesh.getV().size();
      mesh.getV().push_back(_vss[i]);
    }
  }
  //generate vertex
  EdgeMap eMap;
  buildEdge(eMap);
  for(std::map<std::pair<int,int>,Edge,EdgeMap::LSS>::iterator beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++) {
    std::pair<int,int> k=beg->first;
    Edge& e=beg->second;
    e._interNode=-1;
    if((MCVal[k.first] <  0.0f && MCVal[k.second] >= 0.0f) ||
        (MCVal[k.first] >= 0.0f && MCVal[k.second] <  0.0f)) {
      e._interNode=(int)mesh.getV().size();
      mesh.getV().push_back((_vss[k.second]*MCVal[k.first]-_vss[k.first]*MCVal[k.second])/
                            (MCVal[k.first]-MCVal[k.second]));
    }
  }
  //generate index
  int nrT[8]= {1,2,2,1,2,1,1,0};
  int idT[8][6]= {
    {0,1,2, -1,-1,-1},
    {1,5,3, 1,2,5},
    {0,3,4, 0,4,2},
    {2,5,4, -1,-1,-1},

    {0,4,5, 0,1,4},
    {1,4,3, -1,-1,-1},
    {0,3,5, -1,-1,-1},
    {-1,-1,-1, -1,-1,-1},
  };
  for(int i=0; i<(int)_iss.size(); i++) {
    const Vec3i& I=_iss[i];
    unsigned char type=0;
    if(MCVal[I[0]] >= 0.0f)type+=1;
    if(MCVal[I[1]] >= 0.0f)type+=2;
    if(MCVal[I[2]] >= 0.0f)type+=4;
    for(int t=0; t<nrT[type]; t++) {
      Vec3i IBK;
      for(int V=0; V<3; V++) {
        int VBK=idT[type][t*3+V];
        if(VBK < 3)
          IBK[V]=vid[I[VBK]];
        else {
          std::pair<int,int> e((int)I[VBK-3],(int)I[(VBK-2)%3]);
          if(e.first > e.second)std::swap(e.first,e.second);
          IBK[V]=eMap._ess[e]._interNode;
        }
      }
      ASSERT(compGE(IBK,Vec3i::Zero()));
      mesh.getI().push_back(IBK);
    }
  }
  mesh.smooth();
}
//smooth
void OBJMESH::makeUnique()
{
  std::vector<bool> valid;
  valid.assign(_vss.size(),false);
  for(int i=0; i<int(_iss.size()); i++) {
    valid[_iss[i][0]]=true;
    valid[_iss[i][1]]=true;
    if(getDim() == 3)
      valid[_iss[i][2]]=true;
  }

  std::map<sizeType,sizeType> remap;
  sizeType index=0;
  for(int i=0; i<(int)_vss.size(); i++)
    if(valid[i]) {
      _vss[index]=_vss[i];
      remap[i]=index++;
    }
  _vss.resize(index);

  for(int i=0; i<int(_iss.size()); i++) {
    _iss[i][0]=remap[_iss[i][0]];
    _iss[i][1]=remap[_iss[i][1]];
    if(getDim() == 3)
      _iss[i][2]=remap[_iss[i][2]];
  }
}
void OBJMESH::makeUniform()
{
  //initialize edge
  EdgeMap eMap;
  buildEdge(eMap);
  //make uniform
  std::vector<bool> visited(_iss.size(),false);
  for(int i=0; i<(int)visited.size(); i++)
    if(!visited[i]) {
      std::stack<int> queue;
      queue.push(i);
      visited[i]=true;
      while(!queue.empty()) {
        int ti=queue.top();
        queue.pop();
        const Vec3i& I=_iss[ti];
        std::pair<int,int> e;
        for(int eid=0; eid<3; eid++) {
          e.first=(int)I[eid];
          e.second=(int)I[(eid+1)%3];
          if(e.first > e.second)std::swap(e.first,e.second);
          const Edge& edg=eMap._ess[e];
          for(int tid=0; tid<(int)edg._tris.size(); tid++) {
            int tj=edg._tris[tid];
            if(tj != ti && !visited[tj]) {
              makeUniform(ti,tj,e.first,e.second);
              queue.push(tj);
              visited[tj]=true;
            }
          }
        }
      }
    }
}
void OBJMESH::makeUniform(int i,int j,int v0,int v1)
{
  int v0i,v1i,v0j,v1j;
  for(int d=0; d<3; d++) {
    if(_iss[i][d] == v0)v0i=d;
    if(_iss[i][d] == v1)v1i=d;
    if(_iss[j][d] == v0)v0j=d;
    if(_iss[j][d] == v1)v1j=d;
  }
  bool isI=(v0i+1)%3 == v1i;
  bool isJ=(v0j+1)%3 == v1j;
  if(isI == isJ) {
    std::swap(_iss[j][1],_iss[j][2]);
    if(_tex.nrUV() > 0)
      std::swap(_tex._fuv[j][1],_tex._fuv[j][2]);
  }
}
void OBJMESH::smooth()
{
  T feps=ScalarUtil<T>::scalar_eps();
  //force copy and restrict
  std::vector<PT3,Eigen::aligned_allocator<PT3> > tmpNss;
  _nss.swap(tmpNss);

  std::vector<PT3,Eigen::aligned_allocator<PT3> > tmpTnss;
  _tnss.swap(tmpTnss);

  _nss.resize(_vss.size());
  _tnss.resize(_iss.size());
  std::vector<int> nr_face(_vss.size(),0);

  for(int k=0,sz=(int)_vss.size(); k<sz; k++)
    _nss[k].setConstant(0.0f);

  //iterate face
  if(_dim == 2) {
    for(int k=0,sz=(int)_iss.size(); k<sz; k++) {
      const Vec3i &f=_iss[k];

      const PT3 e0=_vss[f.y()]-_vss[f.x()];
      PT3 n(e0.y(),-e0.x(),0.0f);
      if(n.norm() > feps)
        n.normalize();
      else
        n=PT3::Zero();

      _tnss[k]=n;
      nr_face[f.x()]++;
      _nss[f.x()]+=n;
      nr_face[f.y()]++;
      _nss[f.y()]+=n;
    }
  } else {
    for(int k=0,sz=(int)_iss.size(); k<sz; k++) {
      const Vec3i &f=_iss[k];

      const PT3 e0=_vss[f.y()]-_vss[f.x()];
      const PT3 e1=_vss[f.z()]-_vss[f.x()];
      //const PT3 e2=_vss[f.z()]-_vss[f.y()];

      PT3 n=(e0).cross(e1);
      if(n.norm() > feps)
        n.normalize();
      else
        n=PT3::Zero();

      const T a0=getAngle3D<T>(_vss[f.y()]-_vss[f.x()],_vss[f.z()]-_vss[f.x()]);
      const T a1=getAngle3D<T>(_vss[f.z()]-_vss[f.y()],_vss[f.x()]-_vss[f.y()]);
      const T a2=getAngle3D<T>(_vss[f.x()]-_vss[f.z()],_vss[f.y()]-_vss[f.z()]);

      _tnss[k]=n;
      nr_face[f.x()]++;
      _nss[f.x()]+=n*a0;
      nr_face[f.y()]++;
      _nss[f.y()]+=n*a1;
      nr_face[f.z()]++;
      _nss[f.z()]+=n*a2;
    }
  }

  //compute normal
  for(int k=0,sz=(int)_nss.size(); k<sz; k++) {
    if(_nss[k].norm() > feps)
      _nss[k].normalize();
    else
      _nss[k]=PT3::Zero();
  }
}
void OBJMESH::insideOut()
{
  for(sizeType i=0; i<(sizeType)_iss.size(); i++) {
    std::swap(_iss[i][0],_iss[i][1]);
    if(_tex.nrUV() > 0)
      std::swap(_tex._fuv[i][0],_tex._fuv[i][1]);
  }
}
//topology
const OBJMESH::Edge& OBJMESH::getE(int a,int b,const EdgeMap& eMap) const
{
  ASSERT(!eMap._ess.empty())
  if(a<b) return eMap._ess.find(std::pair<int,int>(a,b))->second;
  else return eMap._ess.find(std::pair<int,int>(b,a))->second;
}
void OBJMESH::buildEdge(EdgeMap& eMap) const
{
  ASSERT_MSG(getDim() == 3,"Edge building is valid only for 3D mesh!")
  //edge
  eMap._ess.clear();
  for(int k=0,sz=(int)_iss.size(); k<sz; k++) {
    const Vec3i &f=_iss[k];
    addEdge((int)f.x(),(int)f.y(),k,eMap);
    addEdge((int)f.y(),(int)f.z(),k,eMap);
    addEdge((int)f.x(),(int)f.z(),k,eMap);
  }

  //edge normal
  for(std::map<std::pair<int,int>,Edge,EdgeMap::LSS>::iterator
      begin=eMap._ess.begin(),end=eMap._ess.end(); begin != end; begin++) {
    begin->second._nor.normalize();
  }
}
void OBJMESH::buildKRingV(std::vector<std::map<int,int> >& KRing,int r) const
{
  EdgeMap eMap;
  buildEdge(eMap);
  int nrV=(int)_vss.size();
  //initialize 1-ring triangle
  std::vector<std::set<int> > oneRingV;
  oneRingV.assign(nrV,std::set<int>());
  KRing.assign(nrV,std::map<int,int>());
  for(int i=0; i<nrV; i++)
    KRing[i][i]=-1;  //exclude self
  for(std::map<std::pair<int,int>,Edge,EdgeMap::LSS>::const_iterator
      beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++) {
    oneRingV[beg->first.first].insert(beg->first.second);
    oneRingV[beg->first.second].insert(beg->first.first);
    findInsertV(KRing[beg->first.first],beg->first.second,0);
    findInsertV(KRing[beg->first.second],beg->first.first,0);
  }
  //compute r ring triangle
  for(int i=1; i<r; i++) {
    std::vector<std::map<int,int> > KRingOld=KRing;
    for(int v=0; v<nrV; v++) {
      const std::map<int,int>& Ring=KRingOld[v];
      for(std::map<int,int>::const_iterator
          beg=Ring.begin(),end=Ring.end();
          beg!=end; beg++) {
        findInsert(KRing[v],oneRingV[beg->first],i);
      }
    }
  }
}
void OBJMESH::buildKRing(std::vector<std::map<int,int> >& KRing,int r) const
{
  EdgeMap eMap;
  buildEdge(eMap);
  int nrV=(int)_vss.size();
  //initialize 1-ring triangle
  std::vector<std::set<int> > oneRingT;
  KRing.assign(nrV,std::map<int,int>());
  oneRingT.assign(nrV,std::set<int>());
  for(std::map<std::pair<int,int>,Edge,EdgeMap::LSS>::const_iterator
      beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++) {
    oneRingT[beg->first.first].insert(beg->second._tris.begin(),beg->second._tris.end());
    oneRingT[beg->first.second].insert(beg->second._tris.begin(),beg->second._tris.end());
    findInsert(KRing[beg->first.first],beg->second._tris,0);
    findInsert(KRing[beg->first.second],beg->second._tris,0);
  }
  //compute r ring triangle
  for(int i=1; i<r; i++) {
    std::vector<std::map<int,int> > KRingOld=KRing;
    for(int v=0; v<nrV; v++) {
      const std::map<int,int>& Ring=KRingOld[v];
      for(std::map<int,int>::const_iterator
          beg=Ring.begin(),end=Ring.end();
          beg!=end; beg++) {
        const Vec3i& t=_iss[beg->first];
        findInsert(KRing[v],oneRingT[t[0]],i);
        findInsert(KRing[v],oneRingT[t[1]],i);
        findInsert(KRing[v],oneRingT[t[2]],i);
      }
    }
  }
}
void OBJMESH::findInsertV(std::map<int,int>& Ring,const int& v,int currR) const
{
  if(Ring.find(v) == Ring.end())
    Ring.insert(std::pair<int,int>(v,currR));
}
OBJMESH OBJMESH::cutOpen(T dihedral)
{
  EdgeMap eMap;
  buildEdge(eMap);
  typedef std::map<std::pair<int,int>,Edge,EdgeMap::LSS> ELIST;

  //cluster
  smooth();
  DisjointSet<char> cluster(_iss.size());
  for(ELIST::const_iterator beg=eMap._ess.begin(),end=eMap._ess.end(); beg!=end; beg++) {
    const std::vector<int>& tris=beg->second._tris;
    for(int i=1; i<(int)tris.size(); i++)
      if(_tnss[tris[i]].dot(_tnss[tris[0]]) > dihedral)
        cluster.joinSafe(tris[i],tris[0]);
  }

  //find patch
  INFOV("Cut into %ld patches!",cluster.numSets())
  std::map<int,std::vector<sizeType> > patches;
  for(int i=0; i<(int)_iss.size(); i++)
    patches[cluster.find(i)].push_back(i);
  ASSERT((sizeType)patches.size() == cluster.numSets())

  //output
  int pid=0;
  OBJMESH output,tmp;
  for(std::map<int,std::vector<sizeType> >::const_iterator
      beg=patches.begin(),end=patches.end(); beg!=end; beg++) {
    tmp=*this;
    tmp._iss.clear();
    const std::vector<sizeType>& patch=beg->second;
    for(int f=0; f<(int)patch.size(); f++)
      tmp._iss.push_back(_iss[patch[f]]);
    tmp.makeUnique();

    std::ostringstream oss;
    oss << (pid++);
    output.addMesh(tmp,oss.str());
  }
  output.smooth();
  return output;
}
//dimension
const int& OBJMESH::getDim() const
{
  return _dim;
}
int& OBJMESH::getDim()
{
  return _dim;
}
void OBJMESH::addEdge(int a,int b,int tri,EdgeMap& eMap) const
{
  if(a>b)std::swap(a,b);
  Edge& edg=eMap._ess[std::pair<int,int>(a,b)];
  edg._nor+=_tnss[tri];
  edg._tris.push_back(tri);
}
OBJMESH::T OBJMESH::estimateNormalLength() const
{
  T total=0;
  for(int i=0; i<(int)_iss.size(); i++)
    total+=getArea(i);
  total/=(T)_iss.size();
  if(_dim == 2)
    return total;
  else return std::sqrt(total);
}

#endif
