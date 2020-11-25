#ifdef OBJMESH

class OBJMESH
{
public:
  typedef SCALAR_NAME T;
  typedef Eigen::Matrix<T,2,1> PT2;
  typedef Eigen::Matrix<T,3,1> PT3;
  typedef Eigen::Matrix<T,4,1> PT4;
  typedef Eigen::Matrix<T,3,3> MAT3;
  typedef Eigen::Matrix<double,2,1> PT2D;
  typedef Eigen::Matrix<double,3,1> PT3D;
  struct Edge {
    //subdivision data
    //location of this edge in the
    //original mesh:
    //case 1: _subdId.second == -1
    //	if this edge is inside a triangle,
    //	in this case _subdId.first is the
    //	offset in the original _iss list
    //case 2: _subdId.second != -1
    //	if this edge is a subedge of the
    //	original mesh in this case _subdId
    //	store the edge-end vertex offset in original vss
    std::pair<int,int> _subdId;
    //offset in the subdivided _vss of
    //the vertex created for that edge
    int _interNode;
    //info data
    PT3 _nor;
    std::vector<int> _tris;
    Edge();
  };
  struct EdgeMap {
    friend class OBJMESH;
    struct LSS {
      bool operator()(const std::pair<int,int>& a,const std::pair<int,int>& b) const;
    };
    std::map<std::pair<int,int>,Edge,LSS> _ess;	//edge
  };
  struct PovTexture {
    int nrMAT() const;
    const int& FMAT(int i) const;
    const std::string& MAT(int i) const;
    int nrUV() const;
    const Vec3i& FUV(int i) const;
    const PT2& UV(int i) const;
    //material
    std::vector<int> _fmat;
    std::vector<std::string> _mat;
    //uv
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _fuv;
    std::vector<PT2,Eigen::aligned_allocator<PT2> > _uv;
  };
  OBJMESH();
  OBJMESH(const int& id);
  void setDim(int dim);
  virtual ~OBJMESH();
  template <typename T2>
  void cast(typename ObjMeshTraits<T2>::Type& other) const {
    other.getI()=getI();
    other.getV().resize(_vss.size());
    for(sizeType i=0; i<(sizeType)_vss.size(); i++)
      other.getV()[i]=_vss[i].cast<T2>();
    other.getDim()=getDim();
    other.smooth();
  }
  bool read(std::istream &is,bool move=true,bool scale=true,int assertFaceType=-1,int assertTriangle=-1);
  bool readBinary(std::istream& is);
  bool write(std::ostream &os) const;
  bool write(const std::string& path) const;
  bool writePov(const std::string& path,bool normal=false,const PovTexture* tex=NULL) const;
  bool writePov(std::ostream &os,bool normal=false,const PovTexture* tex=NULL) const;
  bool writePBRT(std::ostream &os) const;
  bool writeVTK(const std::string& path,bool binary,bool normal=false,bool vertexNormal=false,const PT4* color=NULL,const std::vector<T>* cellColor=NULL) const;
  bool writeVTK(VTKWriter<T>& os,bool normal=false,bool vertexNormal=false,const PT4* color=NULL,const std::vector<T>* cellColor=NULL) const;
  bool writeVTK3D(VTKWriter<T>& os,bool normal=false,bool vertexNormal=false,const PT4* color=NULL,const std::vector<T>* cellColor=NULL) const;
  bool writeVTK2D(VTKWriter<T>& os,bool normal=false,bool vertexNormal=false,const PT4* color=NULL,const std::vector<T>* cellColor=NULL) const;
  bool writeBinary(std::ostream &os) const;
  bool write(std::ostream &vs,std::ostream &fs,int &index) const;
  void writeCsv(const std::string& path) const;
  void addMesh(const OBJMESH& mesh,const std::string& g);
  void addMesh(const OBJMESH& mesh);
  PT3 getTC(int i) const;
  const PT3& getV(int i) const;
  const Vec3i& getI(int i) const;
  const int& getIG(int i) const;
  const PT3& getN(int i) const;
  const PT3& getTN(int i) const;
  std::map<std::string,int>& getGS();
  std::map<int,std::string>& getIGS();
  PT3& getV(int i);
  Vec3i& getI(int i);
  int& getIG(int i);
  PT3& getN(int i);
  PT3& getTN(int i);
  T getArea(int i) const;
  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& getV() const;
  const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& getI() const;
  const std::vector<int>& getIG() const;
  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& getN() const;
  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& getFN() const;
  const std::vector<PT3,Eigen::aligned_allocator<PT3> >& getTN() const;
  const std::map<std::string,int>& getGS() const;
  const std::map<int,std::string>& getIGS() const;
  std::string getTG(int i) const;
  int getGId(const std::string& name) const;
  std::vector<PT3,Eigen::aligned_allocator<PT3> >& getV();
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& getI();
  std::vector<int>& getIG();
  std::vector<PT3,Eigen::aligned_allocator<PT3> >& getN();
  std::vector<PT3,Eigen::aligned_allocator<PT3> >& getFN();
  std::vector<PT3,Eigen::aligned_allocator<PT3> >& getTN();
  const MAT3& getT() const;
  MAT3& getT();
  const PT3& getPos() const;
  PT3& getPos();
  const T& getScale() const;
  T& getScale();
  void applyTrans();
  void applyTrans(const PT3& customCtr);
  BBox<T> getBB() const;
  const int& getId() const;
  //physics properties
  T getVolume() const;
  PT3 getCentroid() const;
  PT3 getVolumeCentroid() const;
  PT3& centroidOffset();
  const PT3& centroidOffset() const;
  //simple utility
  T getMass(const T& dens) const;
  void subdivide(const int& nrIter,std::vector<std::pair<int,int> >* vssInfo=NULL);
  void subdivideSingle3D(EdgeMap& eMap,int mod,std::vector<std::pair<int,int> >* vssInfo=NULL) ;
  void subdivideSingle2D();
  void marchingCube(const std::vector<T>& MCVal,OBJMESH& mesh);
  //smooth
  void makeUnique();
  void makeUniform();
  void makeUniform(int i,int j,int v0,int v1);
  void smooth();
  void insideOut();
  //topology
  const Edge& getE(int a,int b,const EdgeMap& eMap) const;
  void buildEdge(EdgeMap& eMap) const;
  void buildKRingV(std::vector<std::map<int,int> >& KRing,int r) const;
  void buildKRing(std::vector<std::map<int,int> >& KRing,int r) const;
  void findInsertV(std::map<int,int>& Ring,const int& v,int currR) const;
  OBJMESH cutOpen(T dihedral);
  //dimension
  const int& getDim() const;
  int& getDim();
  //texture
  PovTexture _tex;
protected:
  void addEdge(int a,int b,int tri,EdgeMap& eMap) const;
  T estimateNormalLength() const;
  std::vector<PT3,Eigen::aligned_allocator<PT3> > _vss;	//vert
  std::vector<PT3,Eigen::aligned_allocator<PT3> > _nss;	//normal
  std::vector<PT3,Eigen::aligned_allocator<PT3> > _fnss;//tree normals per face
  std::vector<PT3,Eigen::aligned_allocator<PT3> > _tnss;//tri_normal
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > _iss;//index
  std::vector<int> _issg;
  std::map<std::string,int> _gss;//groups
  std::map<int,std::string> _igss;//inverse groups
  MAT3 _trans;
  int _id;
  PT3 _pos;
  T _scale;
  PT3 _ctrOff;
  int _dim;
};
template <>
struct ObjMeshTraits<SCALAR_NAME>
{
  typedef OBJMESH Type;
};

#endif
