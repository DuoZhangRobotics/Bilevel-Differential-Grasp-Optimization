#include "GridBasic.h"
#include "Interp.h"
#include "IO.h"

PRJ_BEGIN

//Grid
template <typename T,typename TI,typename TG>
EIGEN_DEVICE_FUNC Grid<T,TI,TG>::Grid():Serializable(typeid(Grid).name()),_szPoint(0,0,0) {}
template <typename T,typename TI,typename TG>
EIGEN_DEVICE_FUNC Grid<T,TI,TG>::~Grid() {}
template <typename T,typename TI,typename TG>
std::shared_ptr<SerializableBase> Grid<T,TI,TG>::copy() const
{
  return std::shared_ptr<SerializableBase>(new Grid());
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::write(std::ostream &os,IOData* dat) const
{
  if(!writeBinaryData(_off,os).good())
    return false;
  if(!writeBinaryData(_szCell,os).good())
    return false;
  if(!writeBinaryData(_invSzCell,os).good())
    return false;
  if(!writeBinaryData(_szPoint,os).good())
    return false;
  if(!writeBinaryData(_bb,os).good())
    return false;
  if(!writeBinaryData(_szPointAligned,os).good())
    return false;
  if(!writeBinaryData(_stride,os).good())
    return false;
  if(!writeBinaryData(_align,os).good())
    return false;
  if(!writeBinaryData(_grid,os).good())
    return false;
  if(!writeBinaryData(_dim,os).good())
    return false;

  return true;
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::read(std::istream &is,IOData* dat)
{
  if(!readBinaryData(_off,is).good())
    return false;
  if(!readBinaryData(_szCell,is).good())
    return false;
  if(!readBinaryData(_invSzCell,is).good())
    return false;
  if(!readBinaryData(_szPoint,is).good())
    return false;
  if(!readBinaryData(_bb,is).good())
    return false;
  if(!readBinaryData(_szPointAligned,is).good())
    return false;
  if(!readBinaryData(_stride,is).good())
    return false;
  if(!readBinaryData(_align,is).good())
    return false;
  if(!readBinaryData(_grid,is).good())
    return false;
  if(!readBinaryData(_dim,is).good())
    return false;

  return true;
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::write(std::ostream &os) const
{
  return write(os,NULL);
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::read(std::istream& is)
{
  return read(is,NULL);
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::writeASCII1D(std::ostream& os) const
{
  os << getNrPoint().x() << std::endl;
  for(sizeType x=0; x<getNrPoint().x(); x++)
    os << (T)get(Vec3i(x,0,0)) << " "  << std::endl;
  return os.good();
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::writeASCII2D(std::ostream& os) const
{
  //y
  //|
  //|
  //|
  //0--------x
  os << getNrPoint().x() << " " << getNrPoint().y() << std::endl;
  for(sizeType y=getNrPoint().y()-1; y>=0; y--) {
    for(sizeType x=0; x<getNrPoint().x(); x++)
      os << (T)get(Vec3i(x,y,0)) << " ";
    os << std::endl;
  }
  return os.good();
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::reset(const Vec3i& nrCell,const BBox<TI>& bb,const T& val,bool center,sizeType align,bool shadow,bool ZFirst)
{
  const IndexType extent=bb.getExtent();
  _dim=3;
  if(extent.z() == 0.0f)_dim--;
  if(extent.y() == 0.0f)_dim--;

  _bb=bb;
  _szPoint=center ? nrCell : (nrCell+Vec3i::Constant(1));
  //fill alignment
  _szPointAligned=_szPoint;
  if(_dim == 3)
    _szPointAligned[2]=((_szPoint[2]+align-1)/align)*align;
  else if(_dim == 2)
    _szPointAligned[1]=((_szPoint[1]+align-1)/align)*align;
  else
    _szPointAligned[0]=((_szPoint[0]+align-1)/align)*align;
  //set all pending dimension to 1
  for(sizeType j=_dim; j<3; j++)
    _szPoint(j)=_szPointAligned(j)=1;
  if((_ZFirst=ZFirst))
    _stride=Vec3i(_szPointAligned.y()*_szPointAligned.z(),_szPointAligned.z(),1);
  else
    _stride=Vec3i(1,_szPointAligned.x(),_szPointAligned.x()*_szPointAligned.y());
  _align=align;

  _off=center ? 0.5f : 0.0f;
  _szCell=IndexType(extent.x()/TI(nrCell.x()),
                    (_dim < 2) ? 0.0f : extent.y()/TI(nrCell.y()),
                    (_dim < 3) ? 0.0f : extent.z()/TI(nrCell.z()));
  _invSzCell=IndexType(1.0f/_szCell.x(),
                       (_dim < 2) ? ScalarUtil<TI>::scalar_max() : 1.0f/_szCell.y(),
                       (_dim < 3) ? ScalarUtil<TI>::scalar_max() : 1.0f/_szCell.z());
  if(!shadow)
    init(val);
  else {
    TG tmp;
    _grid.swap(tmp);
  }
}
template <typename T,typename TI,typename TG>
template <typename T2,typename TI2>
void Grid<T,TI,TG>::makeSameGeometry(const Grid<T2,TI2>& other,bool shadow,bool ZFirst,sizeType align)
{
  BBox<TI> bb;
  bb.copy(other.getBB());
  if(align == -1)
    align=other.getAlign();
  reset(other.getNrCell(),bb,T(),other.isCenter(),align,shadow,ZFirst);
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::expand(const Vec3i& nr,const T& val)
{
  IndexType nrD(_szCell.x()*(scalar)nr.x(),
                _szCell.y()*(scalar)nr.y(),
                _szCell.y()*(scalar)nr.z());

  //expand cell size
  Vec3i nrCell=getNrCell()+nr*2;

  //expand bounding box size
  BBox<TI> bb=getBB();
  bb._minC-=nrD.template cast<TI>();
  bb._maxC+=nrD.template cast<TI>();
  for(sizeType j=2; j>=_dim; j--)
    bb._maxC[j]=bb._minC[j];

  reset(nrCell,bb,val,isCenter(),getAlign(),false);//,_ZFirst,getMove());
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::decimate(const Grid<T,TI>& child,const sizeType& decimate,const T& val,bool uniformScale,sizeType align)
{
  _align=align;

  //we only support node grid decimate
  //ASSERT(child._off == 0.0f);
  ASSERT(decimate > 1);

  if(uniformScale) {
    //cell size
    IndexType cellSz=child.getCellSize()*(TI)decimate;
    Vec3i nrCell=ceilV((IndexType)((child._bb.getExtent().array()/cellSz.array()).matrix()));
    nrCell=compMax(nrCell,Vec3i(2,2,2));	//ensure correct interpolation
    ASSERT(nrCell.x()*nrCell.y()*nrCell.z() > 1)

    BBox<TI> bb;
    bb._minC=child._bb._minC;
    bb._maxC=(cellSz.array()*IndexType((TI)nrCell.x(),(TI)nrCell.y(),(TI)nrCell.z()).array()).matrix()+bb._minC;
    for(sizeType j=child.getDim(); j<3; j++)
      bb._maxC(j)=bb._minC(j);
    reset(nrCell,bb,val,child.isCenter(),align,false);//child.getMove());
  } else {
    Vec3i nrCell=child.getNrCell()/decimate;
    nrCell=compMax(nrCell,Vec3i(2,2,2));
    ASSERT(nrCell.x()*nrCell.y()*nrCell.z() > 1)

    reset(nrCell,child._bb,val,child.isCenter(),align,false);//child.getMove());
  }
}
#define APPLY_WARP  \
if(warpMode)  \
for(sizeType d=0; d<_dim; d++)  \
if((*warpMode)[d] > 0) {  \
while(index[d] < 0)  \
index[d]+=(*warpMode)[d];  \
while(index[d] >= (*warpMode)[d])  \
index[d]-=(*warpMode)[d];}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getAlign() const
{
  return _align;
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getIndexNoAlign(const Vec3i& index) const
{
  return Vec3i(_szPoint.y()*_szPoint.z(),_szPoint.z(),1).dot(index);
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getIndex(Vec3i index,const Vec3i* warpMode) const
{
  APPLY_WARP
  return _stride.dot(index);
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getIndexSafe(Vec3i index,const Vec3i* warpMode) const
{
  APPLY_WARP
  return _stride.dot(compMax(compMin(index,getNrPoint()-Vec3i::Ones()),Vec3i(0,0,0)));
}
template <typename T,typename TI,typename TG>
Vec3i Grid<T,TI,TG>::getIndex(const sizeType& index) const
{
  return Vec3i(index/_stride.x(),(index%_stride.x())/_stride.y(),index%_stride.y());
}
template <typename T,typename TI,typename TG>
const T& Grid<T,TI,TG>::operator[](const Vec3i& index) const
{
  return get(index);
}
template <typename T,typename TI,typename TG>
const T& Grid<T,TI,TG>::operator[](const sizeType& index) const
{
  return _grid[index];
}
template <typename T,typename TI,typename TG>
T& Grid<T,TI,TG>::operator[](const Vec3i& index)
{
  return get(index);
}
template <typename T,typename TI,typename TG>
T& Grid<T,TI,TG>::operator[](const sizeType& index)
{
  return _grid[index];
}
template <typename T,typename TI,typename TG>
const T& Grid<T,TI,TG>::getOptional(Vec3i index,const T& defVal,const Vec3i* warpMode) const
{
  APPLY_WARP
  if(compGE(index,Vec3i(0,0,0)) && compLE(index,getNrPoint()-Vec3i::Ones()))
    return get(index);
  else return defVal;
}
template <typename T,typename TI,typename TG>
const T& Grid<T,TI,TG>::get(const Vec3i& index,const Vec3i* warpMode) const
{
  return _grid[getIndex(index,warpMode)];
}
template <typename T,typename TI,typename TG>
const T& Grid<T,TI,TG>::getSafe(const Vec3i& index,const Vec3i* warpMode) const
{
  return _grid[getIndexSafe(index,warpMode)];
}
template <typename T,typename TI,typename TG>
T& Grid<T,TI,TG>::get(const Vec3i& index,const Vec3i* warpMode)
{
  sizeType id=getIndex(index,warpMode);
  ASSERT(id >= 0 && id < (sizeType)_grid.size())
  return _grid[id];
}
template <typename T,typename TI,typename TG>
T& Grid<T,TI,TG>::getSafe(const Vec3i& index,const Vec3i* warpMode)
{
  return _grid[getIndexSafe(index,warpMode)];
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::isSafeIndex(Vec3i index,const Vec3i* warpMode) const
{
  APPLY_WARP
  return compGE(index,Vec3i(0,0,0)) && compLE(index,getNrPoint()-Vec3i::Ones());
}
#undef APPLY_WARP
template <typename T,typename TI,typename TG>
const Vec3i& Grid<T,TI,TG>::getNrPoint() const
{
  return _szPoint;
}
template <typename T,typename TI,typename TG>
const Vec3i Grid<T,TI,TG>::getMaxInterp() const
{
  Vec3i maxIndex=_szPoint-Vec3i::Constant(1);
  Vec3i maxInterp=maxIndex;
  if(_dim >= 3)
    maxInterp.z()--;
  if(_dim >= 2)
    maxInterp.y()--;
  if(_dim >= 1)
    maxInterp.x()--;
  return maxInterp;
}
template <typename T,typename TI,typename TG>
Vec3i Grid<T,TI,TG>::getNrCell() const
{
  if(_dim == 1)
    return (_off == 0.5f) ? _szPoint : _szPoint-Vec3i(1,0,0);
  else if(_dim == 2)
    return (_off == 0.5f) ? _szPoint : _szPoint-Vec3i(1,1,0);
  else
    return (_off == 0.5f) ? _szPoint : _szPoint-Vec3i::Constant(1);
}
template <typename T,typename TI,typename TG>
const BBox<TI>& Grid<T,TI,TG>::getBB() const
{
  return _bb;
}
template <typename T,typename TI,typename TG>
const typename Grid<T,TI,TG>::IndexType& Grid<T,TI,TG>::getInvCellSize() const
{
  return _invSzCell;
}
template <typename T,typename TI,typename TG>
const typename Grid<T,TI,TG>::IndexType& Grid<T,TI,TG>::getCellSize() const
{
  return _szCell;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::IndexType Grid<T,TI,TG>::getIndexFrac(const IndexType& pos,const Vec3i* warpMode) const
{
  IndexType ret=((pos-_bb._minC).array()*getInvCellSize().array()).matrix()-IndexType::Constant(_off);
  if(_dim < 2)
    ret.y()=0.0f;
  if(_dim < 3)
    ret.z()=0.0f;
  return ret;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::IndexType Grid<T,TI,TG>::getIndexFracSafe(const IndexType& pos,const Vec3i* warpMode) const
{
  IndexType ret=getIndexFrac(pos);
  if(warpMode) {
    for(sizeType d=0; d<_dim; d++)
      if((*warpMode)[d] <= 0)
        ret[d]=std::max<TI>(std::min<TI>(ret[d],(TI)(_szPoint[d]-1-1E-4f)),0);
  } else {
    for(sizeType d=0; d<_dim; d++)
      ret[d]=std::max<TI>(std::min<TI>(ret[d],(TI)(_szPoint[d]-1-1E-4f)),0);
  }
  return ret;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::IndexType Grid<T,TI,TG>::getPt(const Vec3i& cell) const
{
  return _bb._minC+(getCellSize().array()*IndexType(TI(cell.x())+_off,TI(cell.y())+_off,TI(cell.z())+_off).array()).matrix();
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::IndexType Grid<T,TI,TG>::getCellCtr(const Vec3i& cell) const
{
  return _bb._minC+(getCellSize().array()*IndexType(TI(cell.x())+0.5f,TI(cell.y())+0.5f,TI(cell.z())+0.5f).array()).matrix();
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::IndexType Grid<T,TI,TG>::getCellCorner(const Vec3i& cell) const
{
  return _bb._minC+(getCellSize().array()*IndexType(TI(cell.x())     ,TI(cell.y())     ,TI(cell.z())).array()).matrix();
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getSzLinear() const
{
  return _szPointAligned.prod();
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getSzLinearNoAlign() const
{
  return _szPoint.x()*_szPoint.y()*_szPoint.z();
}
//determine 3D interpolation value
#define DETERMINE_VALUE3D(FRAC_FUNC)  \
ASSERT(_dim == 3) \
ASSERT(_szPoint.x() > 1)  \
ASSERT(_szPoint.y() > 1)  \
ASSERT(_szPoint.z() > 1); \
IndexType frac=FRAC_FUNC(pos,warpMode); \
const Vec3i base=floorV(frac);  \
frac-=IndexType((TI)base.x(),(TI)base.y(),(TI)base.z());  \
const sizeType id000=getIndex(base+Vec3i(0,0,0),warpMode);  \
const sizeType id100=getIndex(base+Vec3i(1,0,0),warpMode);  \
const sizeType id010=getIndex(base+Vec3i(0,1,0),warpMode);  \
const sizeType id110=getIndex(base+Vec3i(1,1,0),warpMode);  \
const sizeType id001=getIndex(base+Vec3i(0,0,1),warpMode);  \
const sizeType id101=getIndex(base+Vec3i(1,0,1),warpMode);  \
const sizeType id011=getIndex(base+Vec3i(0,1,1),warpMode);  \
const sizeType id111=getIndex(base+Vec3i(1,1,1),warpMode);
//determine 2D interpolation value
#define DETERMINE_VALUE2D(FRAC_FUNC)  \
ASSERT(_dim == 2) \
ASSERT(_szPoint.x() > 1)  \
ASSERT(_szPoint.y() > 1); \
IndexType frac=FRAC_FUNC(pos,warpMode); \
const Vec3i base=floorV(frac);  \
frac-=IndexType((TI)base.x(),(TI)base.y(),0.0f);  \
const sizeType id000=getIndex(base+Vec3i(0,0,0),warpMode);  \
const sizeType id100=getIndex(base+Vec3i(1,0,0),warpMode);  \
const sizeType id010=getIndex(base+Vec3i(0,1,0),warpMode);  \
const sizeType id110=getIndex(base+Vec3i(1,1,0),warpMode);
//determine 1D interpolation value
#define DETERMINE_VALUE1D(FRAC_FUNC)  \
ASSERT(_dim == 1) \
ASSERT(_szPoint.x() > 1); \
IndexType frac=FRAC_FUNC(pos,warpMode); \
const Vec3i base=floorV(frac);  \
frac-=IndexType((TI)base.x(),0.0f,0.0f);  \
const sizeType id000=getIndex(base+Vec3i(0,0,0),warpMode);  \
const sizeType id100=getIndex(base+Vec3i(1,0,0),warpMode);
#define DETERMINE_VALUE3D_CUBIC(FRAC_FUNC)  \
ASSERT(_dim == 3) \
ASSERT(_szPoint.x() > 1)  \
ASSERT(_szPoint.y() > 1)  \
IndexType frac=FRAC_FUNC(pos,warpMode); \
const Vec3i base=floorV(frac);  \
frac-=IndexType((TI)base.x(),(TI)base.y(),(TI)base.z());  \
for(int i=-1;i<=2;i++)  \
for(int j=-1;j<=2;j++)  \
  for(int k=-1;k<=2;k++)  \
    index[i+1][j+1][k+1]=getIndexSafe(base+Vec3i(k,j,i),warpMode);
#define DETERMINE_VALUE2D_CUBIC(FRAC_FUNC)  \
ASSERT(_dim == 2) \
ASSERT(_szPoint.x() > 1)  \
ASSERT(_szPoint.y() > 1)  \
IndexType frac=FRAC_FUNC(pos,warpMode); \
const Vec3i base=floorV(frac);  \
frac-=IndexType((TI)base.x(),(TI)base.y(),0.0f);  \
for(int j=-1;j<=2;j++)  \
for(int k=-1;k<=2;k++)  \
  index[j+1][k+1]=getIndexSafe(base+Vec3i(k,j,0),warpMode);
#define DETERMINE_VALUE1D_CUBIC(FRAC_FUNC)  \
ASSERT(_dim == 1) \
ASSERT(_szPoint.x() > 1)  \
IndexType frac=FRAC_FUNC(pos,warpMode); \
const Vec3i base=floorV(frac);  \
frac-=IndexType((TI)base.x(),0.0f,0.0f);  \
for(int k=-1;k<=2;k++)  \
index[k+1]=getIndexSafe(base+Vec3i(k,0,0),warpMode);
//sample value
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe(const IndexType& pos,const Vec3i* warpMode) const
{
  return _dim == 1 ? sampleSafe1D(pos,warpMode) :
         _dim == 2 ? sampleSafe2D(pos,warpMode) :
         sampleSafe3D(pos,warpMode);
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe3D(const IndexType& pos,const Vec3i* warpMode) const
{
  DETERMINE_VALUE3D(getIndexFracSafe)
  return interp3D(_grid[id000],_grid[id100],
                  _grid[id010],_grid[id110],	 //lower z plane
                  _grid[id001],_grid[id101],
                  _grid[id011],_grid[id111],
                  frac.x(),frac.y(),frac.z());//higher z plane
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe2D(const IndexType& pos,const Vec3i* warpMode) const
{
  DETERMINE_VALUE2D(getIndexFracSafe)
  return interp2D(_grid[id000],_grid[id100],
                  _grid[id010],_grid[id110],
                  frac.x(),frac.y());
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe1D(const IndexType& pos,const Vec3i* warpMode) const
{
  DETERMINE_VALUE1D(getIndexFracSafe)
  return interp1D(_grid[id000],_grid[id100],frac.x());
}
//sample stencil
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getSampleStencilSafe(const IndexType& pos,TI* coefs,Vec3i* pts,sizeType* offS,const Vec3i* warpMode) const
{
  return _dim == 1 ? getSampleStencilSafe1D(pos,coefs,pts,offS,warpMode) :
         _dim == 2 ? getSampleStencilSafe2D(pos,coefs,pts,offS,warpMode) :
         getSampleStencilSafe3D(pos,coefs,pts,offS,warpMode);
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getSampleStencilSafe3D(const IndexType& pos,TI* coefs,Vec3i* pts,sizeType* offS,const Vec3i* warpMode) const
{
  DETERMINE_VALUE3D(getIndexFracSafe)
  if(pts) {
    pts[0]=base+Vec3i(0,0,0);
    pts[1]=base+Vec3i(1,0,0);
    pts[2]=base+Vec3i(0,1,0);
    pts[3]=base+Vec3i(1,1,0);
    pts[4]=base+Vec3i(0,0,1);
    pts[5]=base+Vec3i(1,0,1);
    pts[6]=base+Vec3i(0,1,1);
    pts[7]=base+Vec3i(1,1,1);
  }
  if(offS) {
    offS[0]=id000;
    offS[1]=id100;
    offS[2]=id010;
    offS[3]=id110;
    offS[4]=id001;
    offS[5]=id101;
    offS[6]=id011;
    offS[7]=id111;
  }
  return stencil3D(coefs,frac.x(),frac.y(),frac.z());//higher z plane
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getSampleStencilSafe2D(const IndexType& pos,TI* coefs,Vec3i* pts,sizeType* offS,const Vec3i* warpMode) const
{
  DETERMINE_VALUE2D(getIndexFracSafe)
  if(pts) {
    pts[0]=base+Vec3i(0,0,0);
    pts[1]=base+Vec3i(1,0,0);
    pts[2]=base+Vec3i(0,1,0);
    pts[3]=base+Vec3i(1,1,0);
  }
  if(offS) {
    offS[0]=id000;
    offS[1]=id100;
    offS[2]=id010;
    offS[3]=id110;
  }
  return stencil2D(coefs,frac.x(),frac.y());//higher z plane
}
template <typename T,typename TI,typename TG>
sizeType Grid<T,TI,TG>::getSampleStencilSafe1D(const IndexType& pos,TI* coefs,Vec3i* pts,sizeType* offS,const Vec3i* warpMode) const
{
  DETERMINE_VALUE1D(getIndexFracSafe)
  if(pts) {
    pts[0]=base+Vec3i(0,0,0);
    pts[1]=base+Vec3i(1,0,0);
  }
  if(offS) {
    offS[0]=id000;
    offS[1]=id100;
  }
  return stencil1D(coefs,frac.x());//higher z plane
}
//sample value with minmax
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode) const
{
  return _dim == 1 ? sampleSafe1D(pos,minV,maxV,warpMode) :
         _dim == 2 ? sampleSafe2D(pos,minV,maxV,warpMode) :
         sampleSafe3D(pos,minV,maxV,warpMode);
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe3D(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode) const
{
  DETERMINE_VALUE3D(getIndexFracSafe)
  return interp3D(_grid[id000],_grid[id100],
                  _grid[id010],_grid[id110],	 //lower z plane
                  _grid[id001],_grid[id101],
                  _grid[id011],_grid[id111],
                  frac.x(),frac.y(),frac.z(),minV,maxV);//higher z plane
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe2D(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode) const
{
  DETERMINE_VALUE2D(getIndexFracSafe)
  return interp2D(_grid[id000],_grid[id100],
                  _grid[id010],_grid[id110],
                  frac.x(),frac.y(),minV,maxV);
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe1D(const IndexType& pos,T& minV,T& maxV,const Vec3i* warpMode) const
{
  DETERMINE_VALUE1D(getIndexFracSafe)
  return interp1D(_grid[id000],_grid[id100],frac.x(),minV,maxV);
}
//sample value with default
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode) const
{
  return _dim == 1 ? sampleSafe1D(pos,valid,def,warpMode) :
         _dim == 2 ? sampleSafe2D(pos,valid,def,warpMode) :
         sampleSafe3D(pos,valid,def,warpMode);
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe3D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode) const
{
  DETERMINE_VALUE3D(getIndexFracSafe)
  if(valid[id000] == 0 || valid[id100] == 0 ||
      valid[id010] == 0 || valid[id110] == 0 ||
      valid[id001] == 0 || valid[id101] == 0 ||
      valid[id011] == 0 || valid[id111] == 0)
    return def;

  return interp3D(_grid[id000],_grid[id100],
                  _grid[id010],_grid[id110],	 //lower z plane
                  _grid[id001],_grid[id101],
                  _grid[id011],_grid[id111],
                  frac.x(),frac.y(),frac.z());//higher z plane
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe2D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode) const
{
  DETERMINE_VALUE2D(getIndexFracSafe)
  if(valid[id000] == 0 || valid[id100] == 0 ||
      valid[id010] == 0 || valid[id110] == 0)
    return def;

  return interp2D(_grid[id000],_grid[id100],
                  _grid[id010],_grid[id110],
                  frac.x(),frac.y());
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe1D(const IndexType& pos,const Grid<unsigned char,TI>& valid,const T& def,const Vec3i* warpMode) const
{
  DETERMINE_VALUE1D(getIndexFracSafe)
  if(valid[id000] == 0 || valid[id100] == 0)
    return def;

  return interp1D(_grid[id000],_grid[id100],frac.x());
}
//sample gradient
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafeGrad(const IndexType& pos,const Vec3i* warpMode) const
{
  return _dim == 1 ? sampleSafe1DGrad(pos,warpMode) :
         _dim == 2 ? sampleSafe2DGrad(pos,warpMode) :
         sampleSafe3DGrad(pos,warpMode);
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafe3DGrad(const IndexType& pos,const Vec3i* warpMode) const
{
  DETERMINE_VALUE3D(getIndexFracSafe)
  return
    ValueType(interp2D<T,typename IndexType::Scalar>
              (_grid[id100]-_grid[id000],
               _grid[id110]-_grid[id010],
               _grid[id101]-_grid[id001],
               _grid[id111]-_grid[id011],
               frac.y(),frac.z())*(typename EigenTraits<T>::ScalarType)_invSzCell.x(),
              interp2D<T,typename IndexType::Scalar>
              (_grid[id010]-_grid[id000],
               _grid[id110]-_grid[id100],
               _grid[id011]-_grid[id001],
               _grid[id111]-_grid[id101],
               frac.x(),frac.z())*(typename EigenTraits<T>::ScalarType)_invSzCell.y(),
              interp2D<T,typename IndexType::Scalar>
              (_grid[id001]-_grid[id000],
               _grid[id101]-_grid[id100],
               _grid[id011]-_grid[id010],
               _grid[id111]-_grid[id110],
               frac.x(),frac.y())*(typename EigenTraits<T>::ScalarType)_invSzCell.z());
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafe2DGrad(const IndexType& pos,const Vec3i* warpMode) const
{
  DETERMINE_VALUE2D(getIndexFracSafe)
  return ValueType(interp1D<T,typename IndexType::Scalar>(_grid[id100]-_grid[id000],_grid[id110]-_grid[id010],frac.y())*(typename EigenTraits<T>::ScalarType)_invSzCell.x(),
                   interp1D<T,typename IndexType::Scalar>(_grid[id010]-_grid[id000],_grid[id110]-_grid[id100],frac.x())*(typename EigenTraits<T>::ScalarType)_invSzCell.y(),
                   EigenTraits<T>::value());
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafe1DGrad(const IndexType& pos,const Vec3i* warpMode) const
{
  DETERMINE_VALUE1D(getIndexFracSafe)
  return ValueType((_grid[id100]-_grid[id000])*(typename EigenTraits<T>::ScalarType)_invSzCell.x(),
                   EigenTraits<T>::value(),EigenTraits<T>::value());
}
//sample gradient stencil
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencilSafeGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS,const Vec3i* warpMode) const
{
  _dim == 1 ? getSampleStencil1DGrad(pos,stencil,offS,warpMode) :
  _dim == 2 ? getSampleStencil2DGrad(pos,stencil,offS,warpMode) :
  getSampleStencil3DGrad(pos,stencil,offS,warpMode);
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencil3DGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS,const Vec3i* warpMode) const
{
  DETERMINE_VALUE3D(getIndexFracSafe)
  stencil.resize(8,3);
  Eigen::Matrix<TI,4,1> coefs;
  stencil2D(coefs.data(),frac.y(),frac.z());
  coefs*=_invSzCell.x();
  stencil(0,0)=-coefs[0];
  stencil(1,0)= coefs[0];
  stencil(2,0)=-coefs[1];
  stencil(3,0)= coefs[1];
  stencil(4,0)=-coefs[2];
  stencil(5,0)= coefs[2];
  stencil(6,0)=-coefs[3];
  stencil(7,0)= coefs[3];
  stencil2D(coefs.data(),frac.x(),frac.z());
  coefs*=_invSzCell.y();
  stencil(0,1)=-coefs[0];
  stencil(2,1)= coefs[0];
  stencil(1,1)=-coefs[1];
  stencil(3,1)= coefs[1];
  stencil(4,1)=-coefs[2];
  stencil(6,1)= coefs[2];
  stencil(5,1)=-coefs[3];
  stencil(7,1)= coefs[3];
  stencil2D(coefs.data(),frac.x(),frac.y());
  coefs*=_invSzCell.z();
  stencil(0,2)=-coefs[0];
  stencil(4,2)= coefs[0];
  stencil(1,2)=-coefs[1];
  stencil(5,2)= coefs[1];
  stencil(2,2)=-coefs[2];
  stencil(6,2)= coefs[2];
  stencil(3,2)=-coefs[3];
  stencil(7,2)= coefs[3];
  if(offS) {
    offS[0]=id000;
    offS[1]=id100;
    offS[2]=id010;
    offS[3]=id110;
    offS[4]=id001;
    offS[5]=id101;
    offS[6]=id011;
    offS[7]=id111;
  }
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencil2DGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS,const Vec3i* warpMode) const
{
  DETERMINE_VALUE2D(getIndexFracSafe)
  stencil.resize(4,2);
  Eigen::Matrix<TI,2,1> coefs;
  stencil1D(coefs.data(),frac.y());
  coefs*=_invSzCell.x();
  stencil(0,0)=-coefs[0];
  stencil(1,0)= coefs[0];
  stencil(2,0)=-coefs[1];
  stencil(3,0)= coefs[1];
  stencil1D(coefs.data(),frac.x());
  coefs*=_invSzCell.y();
  stencil(0,1)=-coefs[0];
  stencil(2,1)= coefs[0];
  stencil(1,1)=-coefs[1];
  stencil(3,1)= coefs[1];
  if(offS) {
    offS[0]=id000;
    offS[1]=id100;
    offS[2]=id010;
    offS[3]=id110;
  }
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencil1DGrad(const IndexType& pos,Eigen::Matrix<TI,-1,-1>& stencil,sizeType* offS,const Vec3i* warpMode) const
{
  DETERMINE_VALUE1D(getIndexFracSafe)
  stencil.resize(2,1);
  stencil(0,0)=-_invSzCell.x();
  stencil(1,0)= _invSzCell.x();
  if(offS) {
    offS[0]=id000;
    offS[1]=id100;
  }
}
//sample laplace
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::sampleLaplaceSafe3D(const IndexType& pos,MatrixType& lap,const Vec3i* warpMode) const
{
  DETERMINE_VALUE3D(getIndexFracSafe)
  lap.row(0)=ValueType(EigenTraits<T>::value(),
                       interp1D<T,typename IndexType::Scalar>
                       ((_grid[id110]-_grid[id010])-(_grid[id100]-_grid[id000]),
                        (_grid[id111]-_grid[id011])-(_grid[id101]-_grid[id001]),frac.z())*(typename EigenTraits<T>::ScalarType)(_invSzCell.y()*_invSzCell.z()),
                       interp1D<T,typename IndexType::Scalar>
                       ((_grid[id101]-_grid[id001])-(_grid[id100]-_grid[id000]),
                        (_grid[id111]-_grid[id011])-(_grid[id110]-_grid[id010]),frac.y())*(typename EigenTraits<T>::ScalarType)(_invSzCell.y()*_invSzCell.z()));

  lap.row(1)=ValueType(interp1D<T,typename IndexType::Scalar>
                       ((_grid[id110]-_grid[id100])-(_grid[id010]-_grid[id000]),
                        (_grid[id111]-_grid[id101])-(_grid[id011]-_grid[id001]),frac.z())*(typename EigenTraits<T>::ScalarType)(_invSzCell.x()*_invSzCell.z()),
                       EigenTraits<T>::value(),
                       interp1D<T,typename IndexType::Scalar>
                       ((_grid[id011]-_grid[id001])-(_grid[id010]-_grid[id000]),
                        (_grid[id111]-_grid[id101])-(_grid[id110]-_grid[id100]),frac.x())*(typename EigenTraits<T>::ScalarType)(_invSzCell.x()*_invSzCell.z()));

  lap.row(2)=ValueType(interp1D<T,typename IndexType::Scalar>
                       ((_grid[id101]-_grid[id100])-(_grid[id001]-_grid[id000]),
                        (_grid[id111]-_grid[id110])-(_grid[id011]-_grid[id010]),frac.y())*(typename EigenTraits<T>::ScalarType)(_invSzCell.x()*_invSzCell.y()),
                       interp1D<T,typename IndexType::Scalar>
                       ((_grid[id011]-_grid[id010])-(_grid[id001]-_grid[id000]),
                        (_grid[id111]-_grid[id110])-(_grid[id101]-_grid[id100]),frac.x())*(typename EigenTraits<T>::ScalarType)(_invSzCell.x()*_invSzCell.y()),
                       EigenTraits<T>::value());
}
//sample value cubic
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafeCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  return _dim == 1 ? sampleSafe1DCubic(pos,warpMode) :
         _dim == 2 ? sampleSafe2DCubic(pos,warpMode) :
         sampleSafe3DCubic(pos,warpMode);
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe3DCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  sizeType index[4][4][4];
  DETERMINE_VALUE3D_CUBIC(getIndexFracSafe)
  T val=interp3DCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),frac.z(),&valLinear,&minV,&maxV);
  if(compGE(val,minV) && compLE(val,maxV))
    return val;
  else return valLinear;
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe2DCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  sizeType index[4][4];
  DETERMINE_VALUE2D_CUBIC(getIndexFracSafe)
  T val=interp2DCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),&valLinear,&minV,&maxV);
  if(compGE(val,minV) && compLE(val,maxV))
    return val;
  else return valLinear;
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sampleSafe1DCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  sizeType index[4];
  DETERMINE_VALUE1D_CUBIC(getIndexFracSafe)
  T val=interp1DCubic<T,TI,TG>(index,_grid,frac.x(),&valLinear,&minV,&maxV);
  if(compGE(val,minV) && compLE(val,maxV))
    return val;
  else return valLinear;
}
//sample stencil cubic
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencilSafeCubic(const IndexType& pos,TI coefs[4][4][4],sizeType offS[4][4][4],const Vec3i* warpMode) const
{
  _dim == 1 ? getSampleStencilSafe1DCubic(pos,coefs[0][0],offS[0][0],warpMode) :
  _dim == 2 ? getSampleStencilSafe2DCubic(pos,coefs[0],offS[0],warpMode) :
  getSampleStencilSafe3DCubic(pos,coefs,offS,warpMode);
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencilSafe3DCubic(const IndexType& pos,TI coefs[4][4][4],sizeType index[4][4][4],const Vec3i* warpMode) const
{
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  DETERMINE_VALUE3D_CUBIC(getIndexFracSafe)
  T val=interp3DCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),frac.z(),&valLinear,&minV,&maxV);
  stencil3DCubic<TI>(coefs,frac.x(),frac.y(),frac.z(),!(compGE(val,minV) && compLE(val,maxV)));
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencilSafe2DCubic(const IndexType& pos,TI coefs[4][4],sizeType index[4][4],const Vec3i* warpMode) const
{
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  DETERMINE_VALUE2D_CUBIC(getIndexFracSafe)
  T val=interp2DCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),&valLinear,&minV,&maxV);
  stencil2DCubic<TI>(coefs,frac.x(),frac.y(),!(compGE(val,minV) && compLE(val,maxV)));
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSampleStencilSafe1DCubic(const IndexType& pos,TI coefs[4],sizeType index[4],const Vec3i* warpMode) const
{
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  DETERMINE_VALUE1D_CUBIC(getIndexFracSafe)
  T val=interp1DCubic<T,TI,TG>(index,_grid,frac.x(),&valLinear,&minV,&maxV);
  stencil1DCubic<TI>(coefs,frac.x(),!(compGE(val,minV) && compLE(val,maxV)));
}
//sample gradient cubic
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafeGradCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  return _dim == 1 ? sampleSafe1DGradCubic(pos,warpMode) :
         _dim == 2 ? sampleSafe2DGradCubic(pos,warpMode) :
         sampleSafe3DGradCubic(pos,warpMode);
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafe3DGradCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  sizeType index[4][4][4];
  ValueType ret;
  for(sizeType i=0; i<ret.size(); i++)
    ret[i]=EigenTraits<T>::value();
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  DETERMINE_VALUE3D_CUBIC(getIndexFracSafe)
  T val=interp3DCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),frac.z(),&valLinear,&minV,&maxV);
  ret=interp3DGradCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),frac.z(),!(compGE(val,minV) && compLE(val,maxV)));
  for(sizeType i=0; i<ret.size(); i++)
    ret[i]*=_invSzCell[i];
  return ret;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafe2DGradCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  sizeType index[4][4];
  ValueType ret;
  for(sizeType i=0; i<ret.size(); i++)
    ret[i]=EigenTraits<T>::value();
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  DETERMINE_VALUE2D_CUBIC(getIndexFracSafe)
  T val=interp2DCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),&valLinear,&minV,&maxV);
  ret.template block<2,1>(0,0)=interp2DGradCubic<T,TI,TG>(index,_grid,frac.x(),frac.y(),!(compGE(val,minV) && compLE(val,maxV)));
  for(sizeType i=0; i<ret.size(); i++)
    ret[i]*=_invSzCell[i];
  return ret;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::ValueType Grid<T,TI,TG>::sampleSafe1DGradCubic(const IndexType& pos,const Vec3i* warpMode) const
{
  sizeType index[4];
  ValueType ret;
  for(sizeType i=0; i<ret.size(); i++)
    ret[i]=EigenTraits<T>::value();
  T valLinear,minV=EigenTraits<T>::max(),maxV=-minV;
  DETERMINE_VALUE1D_CUBIC(getIndexFracSafe)
  T val=interp1DCubic<T,TI,TG>(index,_grid,frac.x(),&valLinear,&minV,&maxV);
  ret.template block<1,1>(0,0)=interp1DGradCubic<T,TI,TG>(index,_grid,frac.x(),!(compGE(val,minV) && compLE(val,maxV)));
  for(sizeType i=0; i<ret.size(); i++)
    ret[i]*=_invSzCell[i];
  return ret;
}
#undef DETERMINE_VALUE3D
#undef DETERMINE_VALUE2D
#undef DETERMINE_VALUE1D
#undef DETERMINE_VALUE3D_CUBIC
#undef DETERMINE_VALUE2D_CUBIC
#undef DETERMINE_VALUE1D_CUBIC
//simple operation
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::dot(const Grid& other) const
{
  const sizeType n=(sizeType)_grid.size();
  T ret=EigenTraits<T>::value();
  for(sizeType i=0; i<n; i++)
    ret+=EigenTraits<T>::mul(_grid[i],other._grid[i]);
  return ret;
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::getAbsMax() const
{
  const sizeType n=(sizeType)_grid.size();
  T maxV=-EigenTraits<T>::max();
  for(sizeType i=0; i<n; i++) {
    const T& val=_grid[i];
    if(compG(EigenTraits<T>::abs(val),maxV))
      maxV=EigenTraits<T>::abs(val);
  }
  return maxV;
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::sum() const
{
  T ret=EigenTraits<T>::value();
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    ret+=_grid[i];
  return ret;
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::squaredDistTo(const Grid& other) const
{
  T ret=EigenTraits<T>::value();
  const TG& otherG=other._grid;
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    ret+=EigenTraits<T>::mul(_grid[i]-otherG[i],_grid[i]-otherG[i]);
  return ret;
}
template <typename T,typename TI,typename TG>
T Grid<T,TI,TG>::squaredNorm() const
{
  T ret=EigenTraits<T>::value();
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    ret+=EigenTraits<T>::mul(_grid[i],_grid[i]);
  return ret;
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::minMax(T& minV,T& maxV) const
{
  const sizeType n=(sizeType)_grid.size();
  minV=EigenTraits<T>::max();
  maxV=-minV;
  for(sizeType i=0; i<n; i++) {
    const T& val=_grid[i];
    minV=compMin(val,minV);
    maxV=compMax(val,maxV);
  }
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::add(const Grid& other)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    _grid[i]+=other._grid[i];
  return *this;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::add(const T& coef)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    _grid[i]+=coef;
  return *this;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::sub(const Grid& other)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    _grid[i]-=other._grid[i];
  return *this;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::min(const Grid& other)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    _grid[i]=compMin(_grid[i],other._grid[i]);
  return *this;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::sub(const T& coef)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    _grid[i]-=coef;
  return *this;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::addScaled(const Grid& other,const T& coef)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    _grid[i]+=EigenTraits<T>::mul(other._grid[i],coef);
  return *this;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::mul(const T& coef)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++)
    _grid[i]=EigenTraits<T>::mul(_grid[i],coef);
  return *this;
}
template <typename T,typename TI,typename TG>
typename Grid<T,TI,TG>::Grid& Grid<T,TI,TG>::clamp(const T& maxVal)
{
  const sizeType n=(sizeType)_grid.size();
  for(sizeType i=0; i<n; i++) {
    T val=EigenTraits<T>::abs(_grid[i]);
    if(compG(val,maxVal))
      _grid[i]=EigenTraits<T>::mul(_grid[i],EigenTraits<T>::div(maxVal,val));
  }
  return *this;
}
//dimension reduction
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSlice(Grid<T,TI>& slice,const sizeType& dim0,const sizeType& dim1,const sizeType& dim2,const T& x) const
{
  Vec3i nrCell(getNrCell()[dim0],getNrCell()[dim1],0);
  BBox<TI> bb;
  bb._minC.x()=_bb._minC[dim0];
  bb._maxC.x()=_bb._maxC[dim0];
  bb._minC.y()=_bb._minC[dim1];
  bb._maxC.y()=_bb._maxC[dim1];
  bb._minC.z()=0.0f;
  bb._maxC.z()=0.0f;
  slice.reset(nrCell,bb,get(Vec3i(0,0,0)),isCenter(),getAlign(),false);

  T val(EigenTraits<T>::sub(x,_bb._minC[dim2])/(typename EigenTraits<T>::ScalarType)getCellSize()[dim2]);
  sizeType sliceId=std::convert<sizeType>()(floorV(val));

  Vec3i nrPoint=slice.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++) {
      Vec3i id;
      id[dim0]=x;
      id[dim1]=y;
      id[dim2]=sliceId;
      slice.get(Vec3i(x,y,0))=get(id);
    }
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSliceYZ(Grid<T,TI>& slice,const T& x) const
{
  getSlice(slice,1,2,0,x);
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSliceXZ(Grid<T,TI>& slice,const T& y) const
{
  getSlice(slice,0,2,1,y);
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::getSliceXY(Grid<T,TI>& slice,const T& z) const
{
  getSlice(slice,0,1,2,z);
}
template <typename T,typename TI,typename TG>
Grid<T,TI> Grid<T,TI,TG>::getSliceYZ(const T& x) const
{
  Grid<T,TI> ret;
  getSlice(ret,1,2,0,x);
  return ret;
}
template <typename T,typename TI,typename TG>
Grid<T,TI> Grid<T,TI,TG>::getSliceXZ(const T& y) const
{
  Grid<T,TI> ret;
  getSlice(ret,0,2,1,y);
  return ret;
}
template <typename T,typename TI,typename TG>
Grid<T,TI> Grid<T,TI,TG>::getSliceXY(const T& z) const
{
  Grid<T,TI> ret;
  getSlice(ret,0,1,2,z);
  return ret;
}
//dangerous methods
template <typename T,typename TI,typename TG>
const TG& Grid<T,TI,TG>::data() const
{
  return _grid;
}
template <typename T,typename TI,typename TG>
TG& Grid<T,TI,TG>::dataRef()
{
  return _grid;
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::setData(const Grid<T,TI>& other)
{
  for(sizeType x=0; x<_szPoint.x(); x++)
    for(sizeType y=0; y<_szPoint.y(); y++)
      for(sizeType z=0; z<_szPoint.z(); z++)
        get(Vec3i(x,y,z))=other.get(Vec3i(x,y,z));
}
template <typename T,typename TI,typename TG>
const T& Grid<T,TI,TG>::get(const sizeType& index) const
{
  return _grid[index];
}
template <typename T,typename TI,typename TG>
T& Grid<T,TI,TG>::get(const sizeType& index)
{
  return _grid[index];
}
template <typename T,typename TI,typename TG>
const Vec3i& Grid<T,TI,TG>::getStride() const
{
  return _stride;
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::isCenter() const
{
  return _off == 0.5f;
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::init(const T& val)
{
  {
    TG tmp;
    _grid.swap(tmp);
  }
  assign(_grid,_szPointAligned.x()*_szPointAligned.y()*_szPointAligned.z(),val);
}
template <typename T,typename TI,typename TG>
const sizeType& Grid<T,TI,TG>::getDim() const
{
  return _dim;
}
template <typename T,typename TI,typename TG>
bool Grid<T,TI,TG>::getZFirst() const
{
  return _ZFirst;
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::swap(Grid& other)
{
  ASSERT(_grid.size() == other._grid.size())
  std::swap(_off,other._off);
  std::swap(_szCell,other._szCell);
  std::swap(_invSzCell,other._invSzCell);
  std::swap(_szPoint,other._szPoint);
  std::swap(_bb,other._bb);
  std::swap(_szPointAligned,other._szPointAligned);
  std::swap(_stride,other._stride);
  std::swap(_align,other._align);
  _grid.swap(other._grid);
  std::swap(_dim,other._dim);
}
template <typename T,typename TI,typename TG>
void Grid<T,TI,TG>::setBB(const BBox<TI>& bb)
{
  ASSERT((bb.getExtent()-_bb.getExtent()).norm() < ScalarUtil<TI>::scalar_eps());
  _bb._minC=bb._minC;
  _bb._maxC=bb._maxC;
}
template <typename T,typename TI,typename TG>
template <typename TT,typename ALLOC>
void Grid<T,TI,TG>::assign(std::vector<TT,ALLOC>& data,sizeType nr,T val)
{
  data.assign(nr,val);
}
template <typename T,typename TI,typename TG>
template <typename TT>
void Grid<T,TI,TG>::assign(Eigen::Matrix<TT,-1,1>& data,sizeType nr,T val)
{
  data.setConstant(nr,val);
}

//MACGrid
template <typename T,typename TI,typename TG>
EIGEN_DEVICE_FUNC MACGrid<T,TI,TG>::MACGrid():Serializable(typeid(MACGrid).name()) {}
template <typename T,typename TI,typename TG>
EIGEN_DEVICE_FUNC MACGrid<T,TI,TG>::~MACGrid() {}
template <typename T,typename TI,typename TG>
std::shared_ptr<SerializableBase> MACGrid<T,TI,TG>::copy() const
{
  return std::shared_ptr<SerializableBase>(new MACGrid());
}
template <typename T,typename TI,typename TG>
bool MACGrid<T,TI,TG>::write(std::ostream &os,IOData* dat) const
{
  if(!writeBinaryData(_bb,os).good())
    return false;
  if(!writeBinaryData(_cellSz,os).good())
    return false;
  if(!writeBinaryData(_nrCell,os).good())
    return false;
  if(!writeBinaryData(_dim,os).good())
    return false;

  bool ret=true;
  if(_dim >= 1 && ret)
    ret=ret && _g[0].write(os);
  if(_dim >= 2 && ret)
    ret=ret && _g[1].write(os);
  if(_dim >= 3 && ret)
    ret=ret && _g[2].write(os);
  return ret;
}
template <typename T,typename TI,typename TG>
bool MACGrid<T,TI,TG>::read(std::istream &is,IOData* dat)
{
  if(!readBinaryData(_bb,is).good())
    return false;
  if(!readBinaryData(_cellSz,is).good())
    return false;
  if(!readBinaryData(_nrCell,is).good())
    return false;
  if(!readBinaryData(_dim,is).good())
    return false;

  bool ret=true;
  if(_dim >= 1 && ret)
    ret=ret && _g[0].read(is);
  if(_dim >= 2 && ret)
    ret=ret && _g[1].read(is);
  if(_dim >= 3 && ret)
    ret=ret && _g[2].read(is);
  return ret;
}
template <typename T,typename TI,typename TG>
bool MACGrid<T,TI,TG>::write(std::ostream &os) const
{
  return write(os,NULL);
}
template <typename T,typename TI,typename TG>
bool MACGrid<T,TI,TG>::read(std::istream &is)
{
  return read(is,NULL);
}
template <typename T,typename TI,typename TG>
template <typename T2,typename TI2,typename TG2>
void MACGrid<T,TI,TG>::reset(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge)
{
  //ASSERT(ref.isCenter())
  _cellSz=IndexType((TI)ref.getCellSize().x(),(TI)ref.getCellSize().y(),(TI)ref.getCellSize().z());
  _nrCell=ref.getNrCell();

  _dim=ref.getDim();
  _bb.copy(ref.getBB());

  if(ref.getDim() >= 1)
    resetGu(ref,shadow,edge);
  if(ref.getDim() >= 2)
    resetGv(ref,shadow,edge);
  if(ref.getDim() >= 3)
    resetGw(ref,shadow,edge);
}
template <typename T,typename TI,typename TG>
template <typename T2,typename TI2,typename TG2>
void MACGrid<T,TI,TG>::makeSameGeometry(const MACGrid<T2,TI2,TG2>& other)
{
  _bb.copy(other.getBB());

  typename MACGrid<T2,TI2,TG2>::IndexType cellSz=other.getCellSize();
  _cellSz=IndexType(cellSz.x(),cellSz.y(),cellSz.z());
  _nrCell=other.getNrCell();
  _dim=other.getDim();

  if(_dim >= 1)
    getGu().makeSameGeometry(other.getGu());
  if(_dim >= 2)
    getGv().makeSameGeometry(other.getGv());
  if(_dim >= 3)
    getGw().makeSameGeometry(other.getGw());
}
template <typename T,typename TI,typename TG>
sizeType MACGrid<T,TI,TG>::getAlign() const
{
  return _g[0].getAlign();
}
template <typename T,typename TI,typename TG>
const Vec3i& MACGrid<T,TI,TG>::getNrCell() const
{
  return _nrCell;
}
template <typename T,typename TI,typename TG>
const BBox<TI>& MACGrid<T,TI,TG>::getBB() const
{
  return _bb;
}
template <typename T,typename TI,typename TG>
const typename MACGrid<T,TI,TG>::IndexType& MACGrid<T,TI,TG>::getInvCellSize() const
{
  return _g[0].getInvCellSize();
}
template <typename T,typename TI,typename TG>
const typename MACGrid<T,TI,TG>::IndexType& MACGrid<T,TI,TG>::getCellSize() const
{
  return _g[0].getCellSize();
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::get(const Vec3i& index) const
{
  if(getDim() == 1)
    return get1D(index);
  else if(getDim() == 2)
    return get2D(index);
  else
    return get3D(index);
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::getSafe(const Vec3i& index) const
{
  if(getDim() == 1)
    return getSafe1D(index);
  else if(getDim() == 2)
    return getSafe2D(index);
  else
    return getSafe3D(index);
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::get1D(const Vec3i& index) const
{
  ASSERT(getDim() == 1)
  return ValueType((_g[0].get(index)+_g[0].get(index+Vec3i(1,0,0)))*0.5f,
                   0.0f,0.0f);
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::getSafe1D(const Vec3i& index) const
{
  return get1D(compMin(compMax(index,Vec3i::Zero()),getNrCell()-Vec3i::Ones()));
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::get2D(const Vec3i& index) const
{
  ASSERT(getDim() == 2)
  return ValueType((_g[0].get(index)+_g[0].get(index+Vec3i(1,0,0)))*0.5f,
                   (_g[1].get(index)+_g[1].get(index+Vec3i(0,1,0)))*0.5f,
                   0.0f);
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::getSafe2D(const Vec3i& index) const
{
  return get2D(compMin(compMax(index,Vec3i::Zero()),getNrCell()-Vec3i::Ones()));
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::get3D(const Vec3i& index) const
{
  ASSERT(getDim() == 3)
  return ValueType((_g[0].get(index)+_g[0].get(index+Vec3i(1,0,0)))*0.5f,
                   (_g[1].get(index)+_g[1].get(index+Vec3i(0,1,0)))*0.5f,
                   (_g[2].get(index)+_g[2].get(index+Vec3i(0,0,1)))*0.5f);
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::getSafe3D(const Vec3i& index) const
{
  return get3D(compMin(compMax(index,Vec3i::Zero()),getNrCell()-Vec3i::Ones()));
}
//safe version
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::sample(const IndexType& pos,const Vec3c* warpMode) const
{
  if(getDim() == 1)
    return sample1D(pos,warpMode);
  else if(getDim() == 2)
    return sample2D(pos,warpMode);
  else return sample3D(pos,warpMode);
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::sampleSafe(const IndexType& pos,const Vec3c* warpMode) const
{
  if(getDim() == 1)
    return sampleSafe1D(pos,warpMode);
  else if(getDim() == 2)
    return sampleSafe2D(pos,warpMode);
  else return sampleSafe3D(pos,warpMode);
}
#define SAMPLE(D,DIM) \
template <typename T,typename TI,typename TG> \
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::sample##DIM(const IndexType& pos,const Vec3c* warpMode) const {  \
ASSERT(getDim() == D);  \
ValueType ret=ValueType::Zero();  \
if(warpMode) {  \
  Vec3i nrCell=Vec3i::Zero(); \
  for(sizeType d=0; d<D; d++) \
    if((*warpMode)[d])  \
      nrCell[d]=_nrCell[d]; \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos,&nrCell);  \
} else {  \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos);  \
} \
return ret; \
} \
template <typename T,typename TI,typename TG> \
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::sampleSafe##DIM(const IndexType& pos,const Vec3c* warpMode) const  \
{ \
ASSERT(getDim() == D);  \
ValueType ret=ValueType::Zero();  \
if(warpMode) {  \
  Vec3i nrCell=Vec3i::Zero(); \
  for(sizeType d=0; d<D; d++) \
    if((*warpMode)[d])  \
      nrCell[d]=_nrCell[d]; \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos,&nrCell);  \
} else {  \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos);  \
} \
return ret; \
}
SAMPLE(1,1D)
SAMPLE(2,2D)
SAMPLE(3,3D)
#undef SAMPLE
//safe version with default value
#define SAMPLE(D,DIM) \
template <typename T,typename TI,typename TG> \
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::sample##DIM(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def,const Vec3c* warpMode) const {  \
ASSERT(getDim() == D);  \
ValueType ret=ValueType::Zero();  \
if(warpMode) {  \
  Vec3i nrCell=Vec3i::Zero(); \
  for(sizeType d=0; d<D; d++) \
    if((*warpMode)[d])  \
      nrCell[d]=_nrCell[d]; \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos,valid.getComp(d),def[d],&nrCell);  \
} else {  \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos,valid.getComp(d),def[d]);  \
} \
return ret; \
} \
template <typename T,typename TI,typename TG> \
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::sampleSafe##DIM(const IndexType& pos,const MACGrid<unsigned char,TI>& valid,const ValueType& def,const Vec3c* warpMode) const  \
{ \
ASSERT(getDim() == D);  \
ValueType ret=ValueType::Zero();  \
if(warpMode) {  \
  Vec3i nrCell=Vec3i::Zero(); \
  for(sizeType d=0; d<D; d++) \
    if((*warpMode)[d])  \
      nrCell[d]=_nrCell[d]; \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos,valid.getComp(d),def[d],&nrCell);  \
} else {  \
  for(sizeType d=0; d<D; d++) \
    ret[d]=_g[d].sampleSafe##DIM(pos,valid.getComp(d),def[d]);  \
} \
return ret; \
}
SAMPLE(1,1D)
SAMPLE(2,2D)
SAMPLE(3,3D)
#undef SAMPLE
//simple operation
template <typename T,typename TI,typename TG>
T MACGrid<T,TI,TG>::dot(const MACGrid& other) const
{
  T ret=0.0f;
  if(getDim() >= 1)
    ret+=_g[0].dot(other._g[0]);
  if(getDim() >= 2)
    ret+=_g[1].dot(other._g[1]);
  if(getDim() >= 3)
    ret+=_g[2].dot(other._g[2]);
  return ret;
}
template <typename T,typename TI,typename TG>
T MACGrid<T,TI,TG>::squaredDistTo(const MACGrid& other) const
{
  T ret=0;
  for(sizeType d=0; d<_dim; d++)
    ret+=_g[d].squaredDistTo(other._g[d]);
  return ret;
}
template <typename T,typename TI,typename TG>
T MACGrid<T,TI,TG>::squaredNorm() const
{
  T ret=0;
  for(sizeType d=0; d<_dim; d++)
    ret+=_g[d].squaredNorm();
  return ret;
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::ValueType MACGrid<T,TI,TG>::getAbsMax() const
{
  ValueType ret=ValueType::Zero();
  if(getDim() >= 1)
    ret.x()=_g[0].getAbsMax();
  if(getDim() >= 2)
    ret.y()=_g[1].getAbsMax();
  if(getDim() >= 3)
    ret.z()=_g[2].getAbsMax();
  return ret;
}
template <typename T,typename TI,typename TG>
void MACGrid<T,TI,TG>::minMax(ValueType& minV,ValueType& maxV) const
{
  minV=maxV=ValueType::Zero();
  if(getDim() >= 1)
    _g[0].minMax(minV.x(),maxV.x());
  if(getDim() >= 2)
    _g[1].minMax(minV.y(),maxV.y());
  if(getDim() >= 3)
    _g[2].minMax(minV.z(),maxV.z());
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::MACGrid& MACGrid<T,TI,TG>::add(const MACGrid& other)
{
  if(getDim() >= 1)
    _g[0].add(other._g[0]);
  if(getDim() >= 2)
    _g[1].add(other._g[1]);
  if(getDim() >= 3)
    _g[2].add(other._g[2]);
  return *this;
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::MACGrid& MACGrid<T,TI,TG>::sub(const MACGrid& other)
{
  if(getDim() >= 1)
    _g[0].sub(other._g[0]);
  if(getDim() >= 2)
    _g[1].sub(other._g[1]);
  if(getDim() >= 3)
    _g[2].sub(other._g[2]);
  return *this;
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::MACGrid& MACGrid<T,TI,TG>::min(const MACGrid& other)
{
  if(getDim() >= 1)
    _g[0].min(other._g[0]);
  if(getDim() >= 2)
    _g[1].min(other._g[1]);
  if(getDim() >= 3)
    _g[2].min(other._g[2]);
  return *this;
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::MACGrid& MACGrid<T,TI,TG>::mul(const T& other)
{
  if(getDim() >= 1)
    _g[0].mul(other);
  if(getDim() >= 2)
    _g[1].mul(other);
  if(getDim() >= 3)
    _g[2].mul(other);
  return *this;
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::MACGrid& MACGrid<T,TI,TG>::addScaled(const MACGrid& other,const T& coef)
{
  if(getDim() >= 1)
    _g[0].addScaled(other._g[0],coef);
  if(getDim() >= 2)
    _g[1].addScaled(other._g[1],coef);
  if(getDim() >= 3)
    _g[2].addScaled(other._g[2],coef);
  return *this;
}
template <typename T,typename TI,typename TG>
typename MACGrid<T,TI,TG>::MACGrid& MACGrid<T,TI,TG>::clamp(const T& maxVal)
{
  if(getDim() >= 1)
    _g[0].clamp(maxVal);
  if(getDim() >= 2)
    _g[1].clamp(maxVal);
  if(getDim() >= 3)
    _g[2].clamp(maxVal);
  return *this;
}
//dangerous methods
template <typename T,typename TI,typename TG>
void MACGrid<T,TI,TG>::setData(const MACGrid<T,TI>& other)
{
  if(_dim>=1)getGu().setData(other.getGu());
  if(_dim>=2)getGv().setData(other.getGv());
  if(_dim>=2)getGw().setData(other.getGw());
}
template <typename T,typename TI,typename TG>
void MACGrid<T,TI,TG>::init(const ValueType& val)
{
  if(getDim() >= 1)
    _g[0].init(val.x());
  if(getDim() >= 2)
    _g[1].init(val.y());
  if(getDim() >= 3)
    _g[2].init(val.z());
}
template <typename T,typename TI,typename TG>
const sizeType& MACGrid<T,TI,TG>::getDim() const
{
  return _dim;
}
template <typename T,typename TI,typename TG>
Grid<T,TI,TG>& MACGrid<T,TI,TG>::getGu()
{
  return _g[0];
}
template <typename T,typename TI,typename TG>
Grid<T,TI,TG>& MACGrid<T,TI,TG>::getGv()
{
  return _g[1];
}
template <typename T,typename TI,typename TG>
Grid<T,TI,TG>& MACGrid<T,TI,TG>::getGw()
{
  return _g[2];
}
template <typename T,typename TI,typename TG>
Grid<T,TI,TG>& MACGrid<T,TI,TG>::getComp(const sizeType& i)
{
  return _g[i];
}
template <typename T,typename TI,typename TG>
const Grid<T,TI,TG>& MACGrid<T,TI,TG>::getGu() const
{
  return _g[0];
}
template <typename T,typename TI,typename TG>
const Grid<T,TI,TG>& MACGrid<T,TI,TG>::getGv() const
{
  return _g[1];
}
template <typename T,typename TI,typename TG>
const Grid<T,TI,TG>& MACGrid<T,TI,TG>::getGw() const
{
  return _g[2];
}
template <typename T,typename TI,typename TG>
const Grid<T,TI,TG>& MACGrid<T,TI,TG>::getComp(const sizeType& i) const
{
  return _g[i];
}
template <typename T,typename TI,typename TG>
void MACGrid<T,TI,TG>::swap(MACGrid& other)
{
  std::swap(_bb,other._bb);
  std::swap(_cellSz,other._cellSz);
  std::swap(_nrCell,other._nrCell);
  std::swap(_dim,other._dim);
  _g[0].swap(other._g[0]);
  _g[1].swap(other._g[1]);
  _g[2].swap(other._g[2]);
}
template <typename T,typename TI,typename TG>
template <typename T2,typename TI2,typename TG2>
void MACGrid<T,TI,TG>::resetGu(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge)
{
  Vec3i nrCell=ref.getNrCell();
  BBox<TI> bb=ref.getBB();
  if(edge) {
    ASSERT(ref.getDim() > 1)
    if(ref.getDim() == 3) {
      bb._minC.x()+=_cellSz.x()*0.5f;
      bb._maxC.x()-=_cellSz.x()*0.5f;
      nrCell.x()--;
    }
  } else {
    if(ref.getDim() >= 2) {
      bb._minC.y()+=_cellSz.y()*0.5f;
      bb._maxC.y()-=_cellSz.y()*0.5f;
      nrCell.y()--;
    }
    if(ref.getDim() >= 3) {
      bb._minC.z()+=_cellSz.z()*0.5f;
      bb._maxC.z()-=_cellSz.z()*0.5f;
      nrCell.z()--;
    }
  }

  _g[0].reset(nrCell,bb,(T)0.0f,false,ref.getAlign(),shadow,true);
}
template <typename T,typename TI,typename TG>
template <typename T2,typename TI2,typename TG2>
void MACGrid<T,TI,TG>::resetGv(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge)
{
  Vec3i nrCell=ref.getNrCell();
  BBox<TI> bb=ref.getBB();
  if(edge) {
    ASSERT(ref.getDim() > 1)
    if(ref.getDim() == 3) {
      bb._minC.y()+=_cellSz.y()*0.5f;
      bb._maxC.y()-=_cellSz.y()*0.5f;
      nrCell.y()--;
    }
  } else {
    if(ref.getDim() >= 1) {
      bb._minC.x()+=_cellSz.x()*0.5f;
      bb._maxC.x()-=_cellSz.x()*0.5f;
      nrCell.x()--;
    }
    if(ref.getDim() >= 3) {
      bb._minC.z()+=_cellSz.z()*0.5f;
      bb._maxC.z()-=_cellSz.z()*0.5f;
      nrCell.z()--;
    }
  }

  _g[1].reset(nrCell,bb,(T)0.0f,false,ref.getAlign(),shadow,true);
}
template <typename T,typename TI,typename TG>
template <typename T2,typename TI2,typename TG2>
void MACGrid<T,TI,TG>::resetGw(const Grid<T2,TI2,TG2>& ref,bool shadow,bool edge)
{
  Vec3i nrCell=ref.getNrCell();
  BBox<TI> bb=ref.getBB();
  if(edge) {
    ASSERT(ref.getDim() > 1)
    bb._minC.z()+=_cellSz.z()*0.5f;
    bb._maxC.z()-=_cellSz.z()*0.5f;
    nrCell.z()--;
  } else {
    if(ref.getDim() >= 1) {
      bb._minC.x()+=_cellSz.x()*0.5f;
      bb._maxC.x()-=_cellSz.x()*0.5f;
      nrCell.x()--;
    }
    if(ref.getDim() >= 2) {
      bb._minC.y()+=_cellSz.y()*0.5f;
      bb._maxC.y()-=_cellSz.y()*0.5f;
      nrCell.y()--;
    }
  }

  _g[2].reset(nrCell,bb,(T)0.0f,false,ref.getAlign(),shadow,true);
}

#define DEF_GRID(T,TI)  \
template struct Grid<T,TI>;  \
template void Grid<T,TI>::makeSameGeometry<sizeType,TI>(const Grid<sizeType,TI>& other,bool shadow,bool ZFirst,sizeType align); \
template void Grid<T,TI>::makeSameGeometry<unsigned char,TI>(const Grid<unsigned char,TI>& other,bool shadow,bool ZFirst,sizeType align); \
template void Grid<T,TI>::makeSameGeometry<scalarF,TI>(const Grid<scalarF,TI>& other,bool shadow,bool ZFirst,sizeType align); \
template void Grid<T,TI>::makeSameGeometry<scalarD,TI>(const Grid<scalarD,TI>& other,bool shadow,bool ZFirst,sizeType align);	\
template void Grid<T,TI>::makeSameGeometry<Vec2i,TI>(const Grid<Vec2i,TI>& other,bool shadow,bool ZFirst,sizeType align);	\
template void Grid<T,TI>::makeSameGeometry<Vec3i,TI>(const Grid<Vec3i,TI>& other,bool shadow,bool ZFirst,sizeType align);

//instance
DEF_GRID(sizeType,sizeType)
DEF_GRID(unsigned char,sizeType)
DEF_GRID(scalarF,sizeType)
DEF_GRID(scalarD,sizeType)
DEF_GRID(Vec2i,sizeType)
DEF_GRID(Vec3i,sizeType)
DEF_GRID(Vec3uc,sizeType)
DEF_GRID(Vec3f,sizeType)
DEF_GRID(Vec3d,sizeType)
template class MACGrid<sizeType,sizeType>;
template class MACGrid<unsigned char,sizeType>;
template class MACGrid<scalarF,sizeType>;
template class MACGrid<scalarD,sizeType>;

DEF_GRID(sizeType,unsigned char)
DEF_GRID(unsigned char,unsigned char)
DEF_GRID(scalarF,unsigned char)
DEF_GRID(scalarD,unsigned char)
DEF_GRID(Vec2i,unsigned char)
DEF_GRID(Vec3i,unsigned char)
DEF_GRID(Vec3uc,unsigned char)
DEF_GRID(Vec3f,unsigned char)
DEF_GRID(Vec3d,unsigned char)
template class MACGrid<sizeType,unsigned char>;
template class MACGrid<unsigned char,unsigned char>;
template class MACGrid<scalarF,unsigned char>;
template class MACGrid<scalarD,unsigned char>;

DEF_GRID(sizeType,scalarF)
DEF_GRID(unsigned char,scalarF)
DEF_GRID(scalarF,scalarF)
DEF_GRID(scalarD,scalarF)
DEF_GRID(Vec2i,scalarF)
DEF_GRID(Vec3i,scalarF)
DEF_GRID(Vec3uc,scalarF)
DEF_GRID(Vec3f,scalarF)
DEF_GRID(Vec3d,scalarF)
template class MACGrid<sizeType,scalarF>;
template class MACGrid<unsigned char,scalarF>;
template class MACGrid<scalarF,scalarF>;
template class MACGrid<scalarD,scalarF>;

DEF_GRID(sizeType,scalarD)
DEF_GRID(unsigned char,scalarD)
DEF_GRID(scalarF,scalarD)
DEF_GRID(scalarD,scalarD)
DEF_GRID(Vec2i,scalarD)
DEF_GRID(Vec3i,scalarD)
DEF_GRID(Vec3uc,scalarD)
DEF_GRID(Vec3f,scalarD)
DEF_GRID(Vec3d,scalarD)
template class MACGrid<sizeType,scalarD>;
template class MACGrid<unsigned char,scalarD>;
template class MACGrid<scalarF,scalarD>;
template class MACGrid<scalarD,scalarD>;

PRJ_END
