#ifndef GRID_BASIC_H
#define GRID_BASIC_H

#include "MathBasic.h"
#include "CollisionDetection.h"
#include "IOFwd.h"

PRJ_BEGIN

template <typename T,typename TI,typename TG=std::vector<T,Eigen::aligned_allocator<T> > >
struct Grid : public Serializable {
public:
  using Serializable::read;
  using Serializable::write;
  typedef T value_type;
  typedef typename Eigen::Matrix<T,3,1> ValueType;
  typedef typename Eigen::Matrix<TI,3,1> IndexType;
  typedef typename Eigen::Matrix<T,3,3> MatrixType;
  EIGEN_DEVICE_FUNC Grid();
  EIGEN_DEVICE_FUNC virtual ~Grid();
  virtual std::shared_ptr<SerializableBase> copy() const;
  virtual bool write(std::ostream &os,IOData* dat) const ;
  virtual bool read(std::istream &is,IOData* dat);
  virtual bool write(std::ostream &os) const;
  virtual bool read(std::istream& is);
  bool writeASCII1D(std::ostream& os) const;
  bool writeASCII2D(std::ostream& os) const;
  void reset(const Vec3i& nrCell,const BBox<TI>& bb,const T& val,bool center=true,sizeType align=4,bool shadow=false,bool ZFirst=true);
  template <typename T2,typename TI2>
  void makeSameGeometry(const Grid<T2,TI2>& other,bool shadow=false,bool ZFirst=true,sizeType align=-1);
  void expand(const Vec3i& nr,const T& val);
  void decimate(const Grid<T,TI>& child,const sizeType& decimate,const T& val,bool uniformScale,sizeType align=4);
public:
  sizeType getAlign() const;
  sizeType getIndexNoAlign(const Vec3i& index) const;
  sizeType getIndex(Vec3i index,const Vec3i* warpMode=NULL) const;
  sizeType getIndexSafe(Vec3i index,const Vec3i* warpMode=NULL) const;
  Vec3i getIndex(const sizeType& index) const;
  const T& operator[](const Vec3i& index) const;
  const T& operator[](const sizeType& index) const;
  T& operator[](const Vec3i& index);
  T& operator[](const sizeType& index);
  const T& getOptional(Vec3i index,const T& defVal,const Vec3i* warpMode=NULL) const;
  const T& get(const Vec3i& index,const Vec3i* warpMode=NULL) const;
  const T& getSafe(const Vec3i& index,const Vec3i* warpMode=NULL) const;
  T& get(const Vec3i& index,const Vec3i* warpMode=NULL);
  T& getSafe(const Vec3i& index,const Vec3i* warpMode=NULL);
  bool isSafeIndex(Vec3i index,const Vec3i* warpMode=NULL) const;
  const Vec3i& getNrPoint() const;
  const Vec3i getMaxInterp() const;
  Vec3i getNrCell() const;
  const BBox<TI>& getBB() const;
  const IndexType& getInvCellSize() const;
  const IndexType& getCellSize() const;
  IndexType getIndexFrac(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  IndexType getIndexFracSafe(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  IndexType getPt(const Vec3i& cell) const;
  IndexType getCellCtr(const Vec3i& cell) const;
  IndexType getCellCorner(const Vec3i& cell) const;
  sizeType getSzLinear() const;
  sizeType getSzLinearNoAlign() const;
  //sample value
  T sampleSafe(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  T sampleSafe3D(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  T sampleSafe2D(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  T sampleSafe1D(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  //sample stencil
  sizeType getSampleStencilSafe(const IndexType& pos,TI* coefs,Vec3i* pts=NULL,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  sizeType getSampleStencilSafe3D(const IndexType& pos,TI* coefs,Vec3i* pts=NULL,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  sizeType getSampleStencilSafe2D(const IndexType& pos,TI* coefs,Vec3i* pts=NULL,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  sizeType getSampleStencilSafe1D(const IndexType& pos,TI* coefs,Vec3i* pts=NULL,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  //sample value with minmax
  T sampleSafe(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode=NULL) const;
  T sampleSafe3D(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode=NULL) const;
  T sampleSafe2D(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode=NULL) const;
  T sampleSafe1D(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode=NULL) const;
  //sample value with default
  T sampleSafe(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode=NULL) const;
  T sampleSafe3D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode=NULL) const;
  T sampleSafe2D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode=NULL) const;
  T sampleSafe1D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode=NULL) const;
  //sample gradient
  ValueType sampleSafeGrad(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  ValueType sampleSafe3DGrad(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  ValueType sampleSafe2DGrad(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  ValueType sampleSafe1DGrad(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  //sample gradient stencil
  void getSampleStencilSafeGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  void getSampleStencil3DGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  void getSampleStencil2DGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  void getSampleStencil1DGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS=NULL,const Vec3i* warpMode=NULL) const;
  //sample laplace
  void sampleLaplaceSafe3D(const IndexType& pos,MatrixType& lap,const Vec3i* warpMode=NULL) const;
  //sample value cubic
  T sampleSafeCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  T sampleSafe3DCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  T sampleSafe2DCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  T sampleSafe1DCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  //sample stencil cubic
  void getSampleStencilSafeCubic(const IndexType& pos,TI coefs[4][4][4],sizeType offS[4][4][4],const Vec3i* warpMode=NULL) const;
  void getSampleStencilSafe3DCubic(const IndexType& pos,TI coefs[4][4][4],sizeType index[4][4][4],const Vec3i* warpMode=NULL) const;
  void getSampleStencilSafe2DCubic(const IndexType& pos,TI coefs[4][4],sizeType index[4][4],const Vec3i* warpMode=NULL) const;
  void getSampleStencilSafe1DCubic(const IndexType& pos,TI coefs[4],sizeType index[4],const Vec3i* warpMode=NULL) const;
  //sample gradient cubic
  ValueType sampleSafeGradCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  ValueType sampleSafe3DGradCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  ValueType sampleSafe2DGradCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  ValueType sampleSafe1DGradCubic(const IndexType& pos,const Vec3i* warpMode=NULL) const;
  //simple operation
  T dot(const Grid& other) const;
  T getAbsMax() const;
  T sum() const;
  T squaredDistTo(const Grid& other) const;
  T squaredNorm() const;
  void minMax(T& minV,T& maxV) const;
  Grid& add(const Grid& other);
  Grid& add(const T& coef);
  Grid& sub(const Grid& other);
  Grid& min(const Grid& other);
  Grid& sub(const T& coef);
  Grid& addScaled(const Grid& other,const T& coef);
  Grid& mul(const T& coef);
  Grid& clamp(const T& maxVal);
  //dimension reduction
  void getSlice(Grid<T,TI>& slice,const sizeType& dim0,const sizeType& dim1,const sizeType& dim2,const T& x) const;
  void getSliceYZ(Grid<T,TI>& slice,const T& x) const;
  void getSliceXZ(Grid<T,TI>& slice,const T& y) const;
  void getSliceXY(Grid<T,TI>& slice,const T& z) const;
  Grid<T,TI> getSliceYZ(const T& x) const;
  Grid<T,TI> getSliceXZ(const T& y) const;
  Grid<T,TI> getSliceXY(const T& z) const;
  //dangerous methods
  const TG& data() const;
  TG& dataRef();
  void setData(const Grid<T,TI>& other);
  const T& get(const sizeType& index) const;
  T& get(const sizeType& index);
  const Vec3i& getStride() const;
  bool isCenter() const;
  void init(const T& val);
  const sizeType& getDim() const;
  bool getZFirst() const;
  void swap(Grid& other);
  void setBB(const BBox<TI>& bb);
protected:
  template <typename TT,typename ALLOC>
  static void assign(std::vector<TT,ALLOC>& data,sizeType nr,T val);
  template <typename TT>
  static void assign(Eigen::Matrix<TT,-1,1>& data,sizeType nr,T val);
  //param
  TI _off;
  IndexType _szCell;
  IndexType _invSzCell;
  Vec3i _szPoint;
  BBox<TI> _bb;
  //memory
  Vec3i _szPointAligned;
  Vec3i _stride;
  sizeType _align;
  TG _grid;
  sizeType _dim;
  bool _ZFirst;
};

template <typename T,typename TI,typename TG=std::vector<T,Eigen::aligned_allocator<T> > >
struct MACGrid : public Serializable {
  friend class dataStructureCL;
public:
  typedef typename Eigen::Matrix<T,3,1> ValueType;
  typedef typename Eigen::Matrix<TI,3,1> IndexType;
  EIGEN_DEVICE_FUNC MACGrid();
  EIGEN_DEVICE_FUNC ~MACGrid();
  virtual std::shared_ptr<SerializableBase> copy() const;
  virtual bool write(std::ostream &os,IOData* dat) const;
  virtual bool read(std::istream &is,IOData* dat);
  virtual bool write(std::ostream &os) const;
  virtual bool read(std::istream &is);
  template <typename T2,typename TI2,typename TG2>
  void reset(const Grid<T2,TI2,TG2>& ref,bool shadow=false,bool edge=false);
  template <typename T2,typename TI2,typename TG2>
  void makeSameGeometry(const MACGrid<T2,TI2,TG2>& other);
public:
  sizeType getAlign() const;
  const Vec3i& getNrCell() const;
  const BBox<TI>& getBB() const;
  const IndexType& getInvCellSize() const;
  const IndexType& getCellSize() const;
  ValueType get(const Vec3i& index) const;
  ValueType getSafe(const Vec3i& index) const;
  ValueType get1D(const Vec3i& index) const;
  ValueType getSafe1D(const Vec3i& index) const;
  ValueType get2D(const Vec3i& index) const;
  ValueType getSafe2D(const Vec3i& index) const;
  ValueType get3D(const Vec3i& index) const;
  ValueType getSafe3D(const Vec3i& index) const;
  //safe version
  ValueType sample(const IndexType& pos,const Vec3c* warpMode=NULL) const;
  ValueType sampleSafe(const IndexType& pos,const Vec3c* warpMode=NULL) const;
#define SAMPLE(D,DIM) \
ValueType sample##DIM(const IndexType& pos,const Vec3c* warpMode=NULL) const; \
ValueType sampleSafe##DIM(const IndexType& pos,const Vec3c* warpMode=NULL) const;
  SAMPLE(1,1D)
  SAMPLE(2,2D)
  SAMPLE(3,3D)
#undef SAMPLE
  //safe version with default value
#define SAMPLE(D,DIM) \
ValueType sample##DIM(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def,const Vec3c* warpMode=NULL) const; \
ValueType sampleSafe##DIM(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def,const Vec3c* warpMode=NULL) const;
  SAMPLE(1,1D)
  SAMPLE(2,2D)
  SAMPLE(3,3D)
#undef SAMPLE
  //simple operation
  T dot(const MACGrid& other) const;
  T squaredDistTo(const MACGrid& other) const;
  T squaredNorm() const;
  ValueType getAbsMax() const;
  void minMax(ValueType& minV,ValueType& maxV) const;
  MACGrid& add(const MACGrid& other);
  MACGrid& sub(const MACGrid& other);
  MACGrid& min(const MACGrid& other);
  MACGrid& mul(const T& other);
  MACGrid& addScaled(const MACGrid& other,const T& coef);
  MACGrid& clamp(const T& maxVal);
  //dangerous methods
  void setData(const MACGrid<T,TI>& other);
  void init(const ValueType& val);
  const sizeType& getDim() const;
  Grid<T,TI,TG>& getGu();
  Grid<T,TI,TG>& getGv();
  Grid<T,TI,TG>& getGw();
  Grid<T,TI,TG>& getComp(const sizeType& i);
  const Grid<T,TI,TG>& getGu() const;
  const Grid<T,TI,TG>& getGv() const;
  const Grid<T,TI,TG>& getGw() const;
  const Grid<T,TI,TG>& getComp(const sizeType& i) const;
  void swap(MACGrid& other);
protected:
  template <typename T2,typename TI2,typename TG2>
  void resetGu(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge);
  template <typename T2,typename TI2,typename TG2>
  void resetGv(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge);
  template <typename T2,typename TI2,typename TG2>
  void resetGw(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge);
  //params
  BBox<TI> _bb;
  IndexType _cellSz;
  Vec3i _nrCell;
  sizeType _dim;
  Grid<T,TI,TG> _g[3];
};

typedef Grid<unsigned char,scalarF> TagFieldF;
typedef Grid<scalarF,scalarF> ScalarFieldF;
typedef Grid<Vec3f,scalarF> VectorFieldF;
typedef MACGrid<unsigned char,scalarF> MACTagFieldF;
typedef MACGrid<scalarF,scalarF> MACVelocityFieldF;

typedef Grid<unsigned char,scalarD> TagFieldD;
typedef Grid<scalarD,scalarD> ScalarFieldD;
typedef Grid<Vec3d,scalarD> VectorFieldD;
typedef MACGrid<unsigned char,scalarD> MACTagFieldD;
typedef MACGrid<scalarD,scalarD> MACVelocityFieldD;

typedef Grid<unsigned char,scalar> TagField;
typedef Grid<scalar,scalar> ScalarField;
typedef Grid<Vec3,scalar> VectorField;
typedef MACGrid<unsigned char,scalar> MACTagField;
typedef MACGrid<scalar,scalar> MACVelocityField;

PRJ_END

#endif
