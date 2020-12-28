#include "GridOp.h"
#include "ObjMesh.h"
#include "ImplicitFuncInterface.h"

#include "IO.h"
#include "Heap.h"
#include <algorithm>
#include <set>

PRJ_BEGIN

//pde
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE(GridType& from,const T& dist,bool preserveSign)
{
  GridType tmp1=from;
  GridType tmp2=from;
  if(from.getDim() == 3)
    solveEikonalPDE3D(from,tmp1,tmp2,dist,preserveSign);
  else if(from.getDim() == 2)
    solveEikonalPDE2D(from,tmp1,tmp2,dist,preserveSign);
  else ASSERT(false);
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smoothEikonalPDE3DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr)
{
  for(sizeType y=0; y<nrPoint.y(); y++)
    for(sizeType z=0; z<nrPoint.z(); z++) {
      T ctr=from.get(Vec3i(x,y,z));
      smooth.get(Vec3i(x,y,z))=ctr/std::sqrt(ctr*ctr+maxCellSzSqr);
    }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smoothEikonalPDE2DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr)
{
  for(sizeType y=0; y<nrPoint.y(); y++) {
    T ctr=from.get(Vec3i(x,y,0));
    smooth.get(Vec3i(x,y,0))=ctr/std::sqrt(ctr*ctr+maxCellSzSqr);
  }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE3DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt)
{
  ValueType invCellSzV=invCellSz.template cast<typename ValueType::Scalar>();
  for(sizeType y=0; y<nrPoint.y(); y++)
    for(sizeType z=0; z<nrPoint.z(); z++) {
      T ctr=from.get(Vec3i(x,y,z));
      T smoothVal=smooth.get(Vec3i(x,y,z));
      ValueType upNorm=(ValueType(upwindingGrad(from.getSafe(Vec3i(x+1,y,z))-ctr,ctr-from.getSafe(Vec3i(x-1,y,z)),smoothVal),
                                  upwindingGrad(from.getSafe(Vec3i(x,y+1,z))-ctr,ctr-from.getSafe(Vec3i(x,y-1,z)),smoothVal),
                                  upwindingGrad(from.getSafe(Vec3i(x,y,z+1))-ctr,ctr-from.getSafe(Vec3i(x,y,z-1)),smoothVal)).array()*invCellSzV.array()).matrix();
      to.get(Vec3i(x,y,z))=ctr-dt*smoothVal*(upNorm.norm()-1.0f);
    }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE2DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt)
{
  ValueType invCellSzV=invCellSz.template cast<typename ValueType::Scalar>();
  for(sizeType y=0; y<nrPoint.y(); y++) {
    T ctr=from.get(Vec3i(x,y,0));
    T smoothVal=smooth.get(Vec3i(x,y,0));
    ValueType upNorm=(ValueType(upwindingGrad(from.getSafe(Vec3i(x+1,y,0))-ctr,ctr-from.getSafe(Vec3i(x-1,y,0)),smoothVal),
                                upwindingGrad(from.getSafe(Vec3i(x,y+1,0))-ctr,ctr-from.getSafe(Vec3i(x,y-1,0)),smoothVal),0.0f).array()*invCellSzV.array()).matrix();
    upNorm.z()=0.0f;
    to.get(Vec3i(x,y,0))=ctr-dt*smoothVal*(upNorm.norm()-1.0f);
  }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE3D(GridType& from,GridType& to,GridType& smooth,const T& dist,bool preserveSign)
{
  const T dt=EigenTraits<T>::constant(from.getCellSize().minCoeff()*0.3f);
  const IndexType invCellSz=from.getInvCellSize();
  const T maxCellSzSqr=EigenTraits<T>::constant(pow(from.getCellSize().maxCoeff(),(TI)2.0f));

  const sizeType nrIter=std::convert<sizeType>()(ceilV(T(dist/dt)));
  const Vec3i nrPoint=from.getNrPoint();

  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<nrPoint.x(); x++)
    smoothEikonalPDE3DLayer(x,from,smooth,nrPoint,maxCellSzSqr);

  for(sizeType iter=0; iter<nrIter; iter++) {
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++)
      solveEikonalPDE3DLayer(x,from,smooth,to,nrPoint,invCellSz,dt);

    if(preserveSign) {
      OMP_PARALLEL_FOR_
      for(sizeType x=0; x<from.getNrPoint().x(); x++)
        for(sizeType y=0; y<from.getNrPoint().y(); y++)
          for(sizeType z=0; z<from.getNrPoint().z(); z++) {
            sizeType id=from.getIndex(Vec3i(x,y,z));
            T& vf=from.get(id);
            T& vt=to.get(id);
            if(vf*vt > 0)
              vf=vt;
            else vf=sgn(vf)*from.getCellSize().maxCoeff()*0.5f;
          }
    } else {
      from.swap(to);
    }
  }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::solveEikonalPDE2D(GridType& from,GridType& to,GridType& smooth,const T& dist,bool preserveSign)
{
  const T dt=std::min<TI>(from.getCellSize().x(),from.getCellSize().y())*0.3f;
  const IndexType invCellSz=from.getInvCellSize();
  const T maxCellSzSqr=std::pow(std::max<TI>(from.getCellSize().x(),from.getCellSize().y()),(TI)2.0f);

  const sizeType nrIter=std::convert<sizeType>()(ceil(dist/dt));
  const Vec3i nrPoint=from.getNrPoint();

  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<nrPoint.x(); x++)
    smoothEikonalPDE2DLayer(x,from,smooth,nrPoint,maxCellSzSqr);

  for(sizeType iter=0; iter<nrIter; iter++) {
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++)
      solveEikonalPDE2DLayer(x,from,smooth,to,nrPoint,invCellSz,dt);

    if(preserveSign) {
      OMP_PARALLEL_FOR_
      for(sizeType x=0; x<from.getNrPoint().x(); x++)
        for(sizeType y=0; y<from.getNrPoint().y(); y++)
          for(sizeType z=0; z<from.getNrPoint().z(); z++) {
            sizeType id=from.getIndex(Vec3i(x,y,z));
            T& vf=from.get(id);
            T& vt=to.get(id);
            if(vf*vt > 0)
              vf=vt;
            else vf=sgn(vf)*from.getCellSize().maxCoeff()*0.5f;
          }
    } else {
      from.swap(to);
    }
  }
}
template <typename T,typename TI,typename TV>
T GridOp<T,TI,TV>::upwindingGrad(T partialPlus,T partialMinus,const T& phi0)
{
  T determinePlus=partialPlus*phi0;
  T determineMinus=partialMinus*phi0;

  if(determineMinus <= 0 && determinePlus <= 0) {} //no modify
  else if(determineMinus >= 0 && determinePlus >= 0)
    partialPlus=partialMinus;
  else if(determineMinus <= 0 && determinePlus >= 0)
    partialPlus=0.0;
  else if(determineMinus >= 0 && determinePlus <= 0) {
    if(std::abs(determinePlus) >= std::abs(determineMinus)) {}//no modify
    else partialPlus=partialMinus;
  }
  return partialPlus;
}
//fast march
template <typename T,typename TI,typename TV>
template<typename VT2>
bool GridOp<T,TI,TV>::solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& normals,Grid<IndexType,TI>& normalExtra)
{
#define CHECK_NEIGH_EXTRA																		\
{																								\
IndexType dir=from.getPt(nPos)-IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z());				\
TI dist=std::abs(dir.dot(IndexType(normals[i].x(),normals[i].y(),normals[i].z())));				\
if(from.get(nPos) == -1.0f){																	\
	knowns.push_back(nPos);																		\
	from.get(nPos)=dist;																		\
    normalExtra.get(nPos)=IndexType(normals[i].x(),normals[i].y(),normals[i].z());				\
}else{																							\
        from.get(nPos)=std::min<TI>(dist,from.get(nPos));										\
    normalExtra.get(nPos)=IndexType(normals[i].x(),normals[i].y(),normals[i].z());				\
}																								\
}

  from.init(-1.0f);

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > knowns;
  for(sizeType i=0; i<(sizeType)nodes.size(); i++) {
    Vec3i base=floor(from.getIndexFracSafe(IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z())));
    Vec3i nPos;
    if(from.getDim() == 2) {
      nPos=base;
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(1,0,0);
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(1,1,0);
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(0,1,0);
      CHECK_NEIGH_EXTRA
    } else if(from.getDim() == 3) {
      nPos=base;
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(1,0,0);
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(1,1,0);
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(0,1,0);
      CHECK_NEIGH_EXTRA

      nPos=base+Vec3i(0,0,1);
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(1,0,1);
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(1,1,1);
      CHECK_NEIGH_EXTRA
      nPos=base+Vec3i(0,1,1);
      CHECK_NEIGH_EXTRA
    }
  }

  return solveEikonalFM<IndexType>(from,knowns,NULL);
#undef CHECK_NEIGH_EXTRA
}
template <typename T,typename TI,typename TV>
template<typename VT2>
bool GridOp<T,TI,TV>::solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes)
{
#define CHECK_NEIGH																				\
{																								\
TI dist=(from.getPt(nPos)-IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z())).norm();			\
if(from.get(nPos) == -1.0f){																	\
	knowns.push_back(nPos);																		\
	from.get(nPos)=dist;																		\
}else{																							\
	from.get(nPos)=std::min<T>(dist,from.get(nPos));											\
}																								\
}

  from.init(-1.0f);

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > knowns;
  for(sizeType i=0; i<(sizeType)nodes.size(); i++) {
    Vec3i base=floor(from.getIndexFracSafe(IndexType(nodes[i].x(),nodes[i].y(),nodes[i].z())));
    Vec3i nPos;
    if(from.getDim() == 2) {
      nPos=base;
      CHECK_NEIGH
      nPos=base+Vec3i(1,0,0);
      CHECK_NEIGH
      nPos=base+Vec3i(1,1,0);
      CHECK_NEIGH
      nPos=base+Vec3i(0,1,0);
      CHECK_NEIGH
    } else if(from.getDim() == 3) {
      nPos=base;
      CHECK_NEIGH
      nPos=base+Vec3i(1,0,0);
      CHECK_NEIGH
      nPos=base+Vec3i(1,1,0);
      CHECK_NEIGH
      nPos=base+Vec3i(0,1,0);
      CHECK_NEIGH

      nPos=base+Vec3i(0,0,1);
      CHECK_NEIGH
      nPos=base+Vec3i(1,0,1);
      CHECK_NEIGH
      nPos=base+Vec3i(1,1,1);
      CHECK_NEIGH
      nPos=base+Vec3i(0,1,1);
      CHECK_NEIGH
    }
  }

  return solveEikonalFM<IndexType>(from,knowns,NULL);
#undef CHECK_NEIGH
}
template <typename T,typename TI,typename TV>
template <typename T2>
bool GridOp<T,TI,TV>::solveEikonalFM(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra,T thres,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary)
{
  if(from.getDim() == 3)
    return solveEikonalFM3D<T2>(from,known,extra,thres,boundary);
  else if(from.getDim() == 2)
    return solveEikonalFM2D<T2>(from,known,extra,thres,boundary);
  else ASSERT(false);
  return false;
}
template <typename T,typename TI,typename TV>
template<typename T2>
bool GridOp<T,TI,TV>::solveEikonalFM3D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra,T thres,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary)
{
  //tags and marks
  TagGridType tags;
  tags.makeSameGeometry(from);
  tags.init(UNKNOWN);

  MarkGridType heapOffsets;
  heapOffsets.makeSameGeometry(from);
  heapOffsets.init(-1);

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > closes;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > heap;

  //initialize heap
  const sizeType nrKnown=(sizeType)known.size();
  for(sizeType i=0; i<nrKnown; i++)
    tags.get(known[i])=KNOWN;
  if(boundary)
    for(sizeType i=0; i<(sizeType)boundary->size(); i++)
      tags.get(boundary->at(i))=BOUNDARY;
  for(sizeType i=0; i<nrKnown; i++) {
    const Vec3i& id=known[i];
    if(tags.get(id)!=KNOWN)
      continue;
    if(id.x() > 0) {
      Vec3i nPos=Vec3i(id.x()-1,id.y(),id.z());
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.x() < tags.getNrPoint().x()-1) {
      Vec3i nPos=Vec3i(id.x()+1,id.y(),id.z());
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.y() > 0) {
      Vec3i nPos=Vec3i(id.x(),id.y()-1,id.z());
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.y() < tags.getNrPoint().y()-1) {
      Vec3i nPos=Vec3i(id.x(),id.y()+1,id.z());
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.z() > 0) {
      Vec3i nPos=Vec3i(id.x(),id.y(),id.z()-1);
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.z() < tags.getNrPoint().z()-1) {
      Vec3i nPos=Vec3i(id.x(),id.y(),id.z()+1);
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
  }

  //add all close to heap
  const sizeType nrClose=closes.size();
  for(sizeType i=0; i<nrClose; i++)
    if(!updateClose3D(from,tags,heapOffsets,heap,closes[i]))
      return false;

  //march until heap empty
  while(!heap.empty()) {
    Vec3i pos=popHeapAbs(from,heapOffsets,heap,Vec3i::Constant(-1));
    tags.get(pos)=KNOWN;

    if(thres > 0.0f && std::abs(from.get(pos)) > thres)
      break;

    if(extra)
      extrapolate3D(from,tags,*extra,pos);

    if(!updateKnown3D(from,tags,heapOffsets,heap,pos))
      return false;
  }

  if(!heap.empty())
    floodFill3D(from,heap,tags);
  return true;
}
template <typename T,typename TI,typename TV>
template<typename T2>
bool GridOp<T,TI,TV>::solveEikonalFM2D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra,T thres,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary)
{
  //tags and marks
  TagGridType tags;
  tags.makeSameGeometry(from);
  tags.init(UNKNOWN);

  MarkGridType heapOffsets;
  heapOffsets.makeSameGeometry(from);
  heapOffsets.init(-1);

  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > closes;
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > heap;

  //initialize heap
  const sizeType nrKnown=(sizeType)known.size();
  for(sizeType i=0; i<nrKnown; i++)
    tags.get(known[i])=KNOWN;
  if(boundary)
    for(sizeType i=0; i<(sizeType)boundary->size(); i++)
      tags.get(boundary->at(i))=BOUNDARY;
  for(sizeType i=0; i<nrKnown; i++) {
    const Vec3i& id=known[i];
    if(tags.get(id)!=KNOWN)
      continue;
    if(id.x() > 0) {
      Vec3i nPos=Vec3i(id.x()-1,id.y(),0);
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.x() < tags.getNrPoint().x()-1) {
      Vec3i nPos=Vec3i(id.x()+1,id.y(),0);
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.y() > 0) {
      Vec3i nPos=Vec3i(id.x(),id.y()-1,0);
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
    if(id.y() < tags.getNrPoint().y()-1) {
      Vec3i nPos=Vec3i(id.x(),id.y()+1,0);
      unsigned char& tag=tags.get(nPos);
      if(tag == UNKNOWN) {
        tag=CLOSE;
        closes.push_back(nPos);
      }
    }
  }

  //add all close to heap
  const sizeType nrClose=closes.size();
  for(sizeType i=0; i<nrClose; i++)
    if(!updateClose2D(from,tags,heapOffsets,heap,closes[i]))
      return false;

  //march until heap empty
  while(!heap.empty()) {
    Vec3i pos=popHeapAbs(from,heapOffsets,heap,Vec3i::Constant(-1));
    tags.get(pos)=KNOWN;

    if(thres > 0.0f && std::abs(from.get(pos)) > thres)
      break;

    if(extra)
      extrapolate2D(from,tags,*extra,pos);

    if(!updateKnown2D(from,tags,heapOffsets,heap,pos))
      return false;
  }

  if(!heap.empty())
    floodFill2D(from,heap,tags);
  return true;
}
template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::extrapolate3D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos)
{
#define CHECK_DIR_EXTRA(DIR,N,P)									\
has=false;															\
if(pos.DIR() > 0 && tags.get(N) == KNOWN)							\
{																	\
    extraDenomTmp=std::abs(from.get(pos)-from.get(N))*				\
				  invCellSz.DIR();									\
	extraNumTmp=extraRef.get(N)*extraDenomTmp;						\
	has=true;														\
}																	\
if(pos.DIR() < nrPoint.DIR()-1 && tags.get(P) == KNOWN)				\
{																	\
    if(!has || std::abs(from.get(P)) < std::abs(from.get(N)))		\
	{																\
        extraDenomTmp=std::abs(from.get(pos)-from.get(P))*			\
					  invCellSz.DIR();								\
		extraNumTmp=extraRef.get(P)*extraDenomTmp;					\
	}																\
	has=true;														\
}																	\
if(has)																\
{																	\
	extraNum+=extraNumTmp;											\
	extraDenom+=extraDenomTmp;										\
}

  Vec3i nrPoint=from.getNrPoint();
  IndexType invCellSz=from.getInvCellSize();
  Grid<T2,TI>& extraRef=extra;

  T2 extraNum=EigenTraits<T2>::value();
  T extraDenom=0.0f;
  bool has;

  Vec3i NX=pos-Vec3i::Unit(0);
  Vec3i PX=pos+Vec3i::Unit(0);
  Vec3i NY=pos-Vec3i::Unit(1);
  Vec3i PY=pos+Vec3i::Unit(1);
  Vec3i NZ=pos-Vec3i::Unit(2);
  Vec3i PZ=pos+Vec3i::Unit(2);

  T2 extraNumTmp=EigenTraits<T2>::value();
  T extraDenomTmp=0.0f;
  CHECK_DIR_EXTRA(x,NX,PX)
  CHECK_DIR_EXTRA(y,NY,PY)
  CHECK_DIR_EXTRA(z,NZ,PZ)
  extraRef.get(pos)=extraNum/std::max<T>(extraDenom,ScalarUtil<T>::scalar_eps());
}
template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::extrapolate2D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos)
{
  Vec3i nrPoint=from.getNrPoint();
  IndexType invCellSz=from.getInvCellSize();
  Grid<T2,TI>& extraRef=extra;

  T2 extraNum=EigenTraits<T2>::value();
  T extraDenom=0.0f;
  bool has;

  Vec3i NX=pos-Vec3i::Unit(0);
  Vec3i PX=pos+Vec3i::Unit(0);
  Vec3i NY=pos-Vec3i::Unit(1);
  Vec3i PY=pos+Vec3i::Unit(1);

  T2 extraNumTmp=EigenTraits<T2>::value();
  T extraDenomTmp=0.0f;
  CHECK_DIR_EXTRA(x,NX,PX)
  CHECK_DIR_EXTRA(y,NY,PY)
  extraRef.get(pos)=extraNum/std::max<T>(extraDenom,ScalarUtil<T>::scalar_eps());
#undef CHECK_DIR_EXTRA
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::floodFill2D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags)
{
#define CHECK_DIR(dir)				\
if(tags.isSafeIndex(dir))			\
{									\
	if(tags.get(dir) == KNOWN)		\
	{tags.get(pos)=KNOWN;			\
	 from.get(pos)=from.get(dir);}	\
	else{id.push_back(dir);}		\
}

  while(!id.empty()) {
    const Vec3i pos=id.back();
    id.pop_back();
    CHECK_DIR(pos+Vec3i(-1,0,0))
    CHECK_DIR(pos+Vec3i(1,0,0))
    CHECK_DIR(pos+Vec3i(0,-1,0))
    CHECK_DIR(pos+Vec3i(0,1,0))
  }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::floodFill3D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags)
{
  while(!id.empty()) {
    const Vec3i pos=id.back();
    id.pop_back();
    CHECK_DIR(pos+Vec3i(-1,0,0))
    CHECK_DIR(pos+Vec3i(1,0,0))
    CHECK_DIR(pos+Vec3i(0,-1,0))
    CHECK_DIR(pos+Vec3i(0,1,0))
    CHECK_DIR(pos+Vec3i(0,0,-1))
    CHECK_DIR(pos+Vec3i(0,0,1))
  }

#undef CHECK_DIR
}
template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateClose2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close)
{
#define CHECK_DIR(DIR,N,P)											\
exist=false;														\
if(close.DIR() > 0 && tags.get(N) == KNOWN){						\
	B[off]=result.get(N);											\
	exist=true;														\
}																	\
if(close.DIR() < tags.getNrPoint().DIR()-1 && tags.get(P) == KNOWN)	\
{																	\
	tmp=result.get(P);												\
	if(exist){														\
        if(std::abs(tmp) < std::abs(B[off]))						\
			B[off]=tmp;												\
	}																\
	else{															\
		B[off]=tmp;													\
	}																\
	exist=true;														\
}																	\
if(exist)															\
        off++;

  //locals
  T B[3];
  sizeType off=0;
  T a,b,c,tmp,det;
  bool exist;
  ValueType invCellSz=tags.getInvCellSize().template cast<typename ValueType::Scalar>();
  ValueType invCellSzSqr=(invCellSz.array()*invCellSz.array()).matrix();
  Vec3i NX=close-Vec3i::Unit(0);
  Vec3i PX=close+Vec3i::Unit(0);
  Vec3i NY=close-Vec3i::Unit(1);
  Vec3i PY=close+Vec3i::Unit(1);

  //check every direction
  CHECK_DIR(x,NX,PX)
  CHECK_DIR(y,NY,PY)

  //check same sign
  bool positive=B[0] > 0.0f;
  for(sizeType i=1; i<off; i++) {
    if(B[0]*B[i]<0.0f) {
      WARNING("In-Out Test Fail For Fast Marching,Mesh Not Watertight!")
      return false;
    }
  }

  //sort the minimal direction
  if(off == 2) {
    if(std::abs(B[1]) < std::abs(B[0]))
      std::swap(B[0],B[1]);
  }

  //try solving
  while(off>0) {
    //build a,b,c
    a=b=0.0f;
    c=(T)-1.0f;
    for(sizeType i=0; i<off; i++) {
      a+=invCellSzSqr[i];
      b-=2.0f*B[i]*invCellSzSqr[i];
      c+=B[i]*B[i]*invCellSzSqr[i];
    }

    det=b*b-4.0f*a*c;
    if(det < 0.0f) {
      off--;
    } else {
      //solution to the quadratic
      a=std::max<T>(2.0f*a,ScalarUtil<T>::scalar_eps());
      result.get(close)=positive ? ((-b+std::sqrt(det))/a) : ((-b-std::sqrt(det))/a);
      result.get(close)=std::abs(result.get(close));
      if(!positive)
        result.get(close)*=-1.0f;

      //update heap
      if(heapOffsets.get(close) < 0)
        pushHeapAbs(result,heapOffsets,heap,close);
      else
        updateHeapAbs(result,heapOffsets,heap,close);

      break;
    }
  }

  if(off == 0) {
    WARNING("Numerically Impossible!")
    return false;
  }
  return true;
}
template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateClose3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close)
{
  //locals
  T B[3];
  sizeType off=0;
  T a,b,c,tmp,det;
  bool exist;
  ValueType invCellSz=tags.getInvCellSize().template cast<typename ValueType::Scalar>();
  ValueType invCellSzSqr=(invCellSz.array()*invCellSz.array()).matrix();
  Vec3i NX=close-Vec3i::Unit(0);
  Vec3i PX=close+Vec3i::Unit(0);
  Vec3i NY=close-Vec3i::Unit(1);
  Vec3i PY=close+Vec3i::Unit(1);
  Vec3i NZ=close-Vec3i::Unit(2);
  Vec3i PZ=close+Vec3i::Unit(2);

  //check every direction
  CHECK_DIR(x,NX,PX)
  CHECK_DIR(y,NY,PY)
  CHECK_DIR(z,NZ,PZ)

  //check same sign
  bool positive=B[0] > 0.0f;
  for(sizeType i=1; i<off; i++) {
    if(B[0]*B[i]<0.0f) {
      WARNING("In-Out Test Fail For Fast Marching,Mesh Not Watertight!")
      return false;
    }
  }

  //sort the minimal direction
  if(off == 3) {
    if(std::abs(B[1]) < std::abs(B[0]))
      std::swap(B[0],B[1]);

    if(std::abs(B[2]) < std::abs(B[1]))
      std::swap(B[1],B[2]);

    if(std::abs(B[1]) < std::abs(B[0]))
      std::swap(B[0],B[1]);
  } else if(off == 2) {
    if(std::abs(B[1]) < std::abs(B[0]))
      std::swap(B[0],B[1]);
  }

  //try solving
  while(off>0) {
    //build a,b,c
    a=b=0.0f;
    c=(T)-1.0f;
    for(sizeType i=0; i<off; i++) {
      a+=invCellSzSqr[i];
      b-=2.0f*B[i]*invCellSzSqr[i];
      c+=B[i]*B[i]*invCellSzSqr[i];
    }

    det=b*b-4.0f*a*c;
    if(det < 0.0f) {
      off--;
    } else {
      //solution to the quadratic
      a=std::max<T>(2.0f*a,ScalarUtil<T>::scalar_eps());
      result.get(close)=positive ? ((-b+std::sqrt(det))/a) : ((-b-std::sqrt(det))/a);
      result.get(close)=std::abs(result.get(close));
      if(!positive)
        result.get(close)*=-1.0f;

      //update heap
      if(heapOffsets.get(close) < 0)
        pushHeapAbs(result,heapOffsets,heap,close);
      else
        updateHeapAbs(result,heapOffsets,heap,close);

      break;
    }
  }

  if(off == 0) {
    WARNING("Numerically Impossible!")
    return false;
  }
  return true;
#undef CHECK_DIR
}
template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateKnown2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id)
{
  bool ret=true;
  if(id.x() > 0) {
    Vec3i nPos=Vec3i(id.x()-1,id.y(),0);
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.x() < tags.getNrPoint().x()-1) {
    Vec3i nPos=Vec3i(id.x()+1,id.y(),0);
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.y() > 0) {
    Vec3i nPos=Vec3i(id.x(),id.y()-1,0);
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.y() < tags.getNrPoint().y()-1) {
    Vec3i nPos=Vec3i(id.x(),id.y()+1,0);
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose2D(result,tags,heapOffsets,heap,nPos);
    }
  }
  return ret;
}
template <typename T,typename TI,typename TV>
bool GridOp<T,TI,TV>::updateKnown3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id)
{
  bool ret=true;
  if(id.x() > 0) {
    Vec3i nPos=Vec3i(id.x()-1,id.y(),id.z());
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.x() < tags.getNrPoint().x()-1) {
    Vec3i nPos=Vec3i(id.x()+1,id.y(),id.z());
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.y() > 0) {
    Vec3i nPos=Vec3i(id.x(),id.y()-1,id.z());
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.y() < tags.getNrPoint().y()-1) {
    Vec3i nPos=Vec3i(id.x(),id.y()+1,id.z());
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.z() > 0) {
    Vec3i nPos=Vec3i(id.x(),id.y(),id.z()-1);
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
    }
  }
  if(id.z() < tags.getNrPoint().z()-1) {
    Vec3i nPos=Vec3i(id.x(),id.y(),id.z()+1);
    unsigned char& tag=tags.get(nPos);
    if(tag != KNOWN && tag!= BOUNDARY) {
      tag=CLOSE;
      ret=ret&&updateClose3D(result,tags,heapOffsets,heap,nPos);
    }
  }
  return ret;
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::reinitialize(GridType& from)
{
  const T thres=from.getCellSize().maxCoeff()*2.0f;
  const Vec3i nrP=from.getNrPoint();
  solveEikonalPDE(from,thres*3.0f,true);
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > known;
  if(from.getDim() == 3) {
    for(sizeType x=0; x<nrP.x(); x++)
      for(sizeType y=0; y<nrP.y(); y++)
        for(sizeType z=0; z<nrP.z(); z++) {
          const T val=from.get(Vec3i(x,y,z));
          if(std::abs(val)<thres)
            known.push_back(Vec3i(x,y,z));
        }
  } else {
    for(sizeType x=0; x<nrP.x(); x++)
      for(sizeType y=0; y<nrP.y(); y++) {
        const T val=from.get(Vec3i(x,y,0));
        if(std::abs(val)<thres)
          known.push_back(Vec3i(x,y,0));
      }
  }
  solveEikonalFM<IndexType>(from,known);
}
//smooth
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smooth(GridType& from,const TagGridType* tagField,unsigned char tag)
{
  GridType tmp=from;
  if(from.getDim() == 3)
    smooth3D(from,tmp,tagField,tag);
  else if(from.getDim() == 2)
    smooth2D(from,tmp,tagField,tag);
  else ASSERT(false);
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smooth3D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag)
{
  const scalar weight=1.0f;//0.6f;
  const Vec3i nrPoint=from.getNrPoint();

  OMP_PARALLEL_FOR_
  for(sizeType x=1; x<nrPoint.x()-1; x++)
    for(sizeType y=1; y<nrPoint.y()-1; y++)
      for(sizeType z=1; z<nrPoint.z()-1; z++) {
        if(!tagField || tagField->get(Vec3i(x,y,z)) == tag)
          tmp.get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z))*(1.0f-weight)+
                                (from.get(Vec3i(x-1,y,z))+from.get(Vec3i(x+1,y,z))+
                                 from.get(Vec3i(x,y-1,z))+from.get(Vec3i(x,y+1,z))+
                                 from.get(Vec3i(x,y,z-1))+from.get(Vec3i(x,y,z+1)))*weight/6.0f;
      }
  from.swap(tmp);
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::smooth2D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag)
{
  const scalar weight=1.0f;//0.6f;
  const Vec3i nrPoint=from.getNrPoint();

  OMP_PARALLEL_FOR_
  for(sizeType x=1; x<nrPoint.x()-1; x++)
    for(sizeType y=1; y<nrPoint.y()-1; y++) {
      if(!tagField || tagField->get(Vec3i(x,y,0)) == tag)
        tmp.get(Vec3i(x,y,0))=from.get(Vec3i(x,y,0))*(1.0f-weight)+
                              (from.get(Vec3i(x-1,y,0))+from.get(Vec3i(x+1,y,0))+
                               from.get(Vec3i(x,y-1,0))+from.get(Vec3i(x,y+1,0)))*weight/4.0f;
    }
  from.swap(tmp);
}
//io
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarGridObj(ObjMesh& mesh,const GridType& grd,bool moveVertex)
{
  ASSERT(grd.getDim() == 2)
  std::vector<Vec3,Eigen::aligned_allocator<Vec3> >& pos=mesh.getV();
  std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& index=mesh.getI();

  const Vec3i nrPoint=grd.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++) {
      IndexType p=grd.getPt(Vec3i(x,y,0));
      if(moveVertex)
        p.z()=(TI)grd.get(Vec3i(x,y,0));
      pos.push_back(p.template cast<scalar>());
    }

  for(sizeType x=0; x<nrPoint.x()-1; x++)
    for(sizeType y=0; y<nrPoint.y()-1; y++) {
#define GI(a,b) ((a)*nrPoint.y()+(b))
      index.push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
      index.push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
#undef GI
    }

  mesh.getDim()=3;
  mesh.smooth();
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarGridVTK(const std::string& path,const GridType& grd,bool moveVertex)
{
  ASSERT(grd.getDim() == 2)
  std::vector<IndexType> pos;
  std::vector<Vec3i> index;

  std::vector<T> pointData;
  const Vec3i nrPoint=grd.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++) {
      IndexType p=grd.getPt(Vec3i(x,y,0));
      if(moveVertex)
        p.z()=(TI)grd.get(Vec3i(x,y,0));
      pos.push_back(p);
      pointData.push_back(grd.get(Vec3i(x,y,0)));
    }

  for(sizeType x=0; x<nrPoint.x()-1; x++)
    for(sizeType y=0; y<nrPoint.y()-1; y++) {
#define GI(a,b) ((a)*nrPoint.y()+(b))
      index.push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
      index.push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
#undef GI
    }

  VTKWriter<TI> writer("2D Grid",path,true);
  writer.appendPoints(pos.begin(),pos.end());
  writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::TRIANGLE);
  writer.appendCustomPointData("value",pointData.begin(),pointData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarBarChartVTK(const std::string& path,const GridType& grd,bool moveVertex,const TI& scale)
{
  //typedef Eigen::Matrix<sizeType,8,1> IDS;

  ASSERT(grd.getDim() == 2)
  IndexType cellSize=grd.getCellSize()*0.5f*scale;
  cellSize.z()=0.0f;

  std::vector<IndexType> pos;
  std::vector<TI> cellData;

  const Vec3i nrPoint=grd.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++) {
      T val=grd.get(Vec3i(x,y,0));
      IndexType zOff(0.0f,0.0f,(TI)val);
      if(!moveVertex)
        zOff.z()=1E-3f;
      IndexType p=grd.getPt(Vec3i(x,y,0));

      pos.push_back(p-cellSize);
      pos.push_back(p+cellSize+zOff);

      cellData.push_back((TI)val);
    }

  VTKWriter<TI> writer("2D Grid BarChart",path,true);
  writer.appendVoxels(pos.begin(),pos.end(),true);
  writer.appendCustomData("value",cellData.begin(),cellData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DMarkGridVTK(const std::string& path,const MarkGridType& grd,bool moveVertex)
{
  ASSERT(grd.getDim() == 2)
  std::vector<IndexType> pos;
  std::vector<Vec3i> index;

  std::vector<T> pointData;
  const Vec3i nrPoint=grd.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++) {
      IndexType p=grd.getPt(Vec3i(x,y,0));
      if(moveVertex)
        p.z()=(T)grd.get(Vec3i(x,y,0));
      pos.push_back(p);
      pointData.push_back((T)grd.get(Vec3i(x,y,0)));
    }

  for(sizeType x=0; x<nrPoint.x()-1; x++)
    for(sizeType y=0; y<nrPoint.y()-1; y++) {
#define GI(a,b) ((a)*nrPoint.y()+(b))
      index.push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
      index.push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
#undef GI
    }

  VTKWriter<TI> writer("2D Grid",path,true);
  writer.appendPoints(pos.begin(),pos.end());
  writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::TRIANGLE);
  writer.appendCustomPointData("value",pointData.begin(),pointData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarGridGradVTK(const std::string& path,const GridType& grd,const sizeType& sampleInterval,const TI& len)
{
  ASSERT(grd.getDim() == 2)
  std::vector<IndexType> pos;
  std::vector<Vec3i> index;

  std::vector<TI> cellData;
  const Vec3i nrPoint=grd.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x+=sampleInterval)
    for(sizeType y=0; y<nrPoint.y(); y+=sampleInterval) {
      IndexType p=grd.getPt(Vec3i(x,y,0));
      p.z()=0.0f;

      ValueType delta=grd.sampleSafeGrad(p);
      if(delta.norm() < ScalarUtil<T>::scalar_eps())
        continue;

      delta.normalize();
      delta*=len;
      IndexType end(p.x()+delta.x(),p.y()+delta.y(),p.z()+delta.z());
      end.z()=0.0f;

      index.push_back(Vec3i(pos.size(),pos.size()+1,0));
      cellData.push_back(grd.get(Vec3i(x,y,0)));
      pos.push_back(p);
      pos.push_back(end);
    }

  VTKWriter<TI> writer("2D Grid Grad",path,false);
  writer.appendPoints(pos.begin(),pos.end());
  writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::LINE);
  writer.appendCustomData("value",cellData.begin(),cellData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DScalarGridGradMagVTK(const std::string& path,const GridType& grd)
{
  ASSERT(grd.getDim() == 2)
  std::vector<IndexType> pos;
  std::vector<Vec3i> index;

  std::vector<T> pointData;
  const Vec3i nrPoint=grd.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++) {
      IndexType p=grd.getPt(Vec3i(x,y,0));
      pos.push_back(p);
      pointData.push_back(grd.sampleSafeGrad(grd.getPt(Vec3i(x,y,0))).norm());
    }

  for(sizeType x=0; x<nrPoint.x()-1; x++)
    for(sizeType y=0; y<nrPoint.y()-1; y++) {
#define GI(a,b) ((a)*nrPoint.y()+(b))
      index.push_back(Vec3i(GI(x,y),GI(x+1,y),GI(x+1,y+1)));
      index.push_back(Vec3i(GI(x,y),GI(x+1,y+1),GI(x,y+1)));
#undef GI
    }

  VTKWriter<TI> writer("2D Grid",path,true);
  writer.appendPoints(pos.begin(),pos.end());
  writer.appendCells(index.begin(),index.end(),VTKWriter<TI>::TRIANGLE);
  writer.appendCustomPointData("value",pointData.begin(),pointData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write3DScalarGridVTK(const std::string& path,const GridType& grd)
{
  ASSERT(grd.getDim() == 3)

  std::vector<TI> pointData;
  const Vec3i nrPoint=grd.getNrPoint();
  for(sizeType z=0; z<nrPoint.z(); z++)
    for(sizeType y=0; y<nrPoint.y(); y++)
      for(sizeType x=0; x<nrPoint.x(); x++)
        pointData.push_back(grd.get(Vec3i(x,y,z)));

  VTKWriter<TI> writer("3D Grid",path,true,grd.getBB(),grd.getNrCell(),grd.isCenter());
  writer.appendDatas("value",pointData.begin(),pointData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DVectorGridVTK(const std::string& path,const Grid<IndexType,TI>& vel,const TI& len)
{
  ASSERT(vel.getDim() == 2)

  std::vector<IndexType,Eigen::aligned_allocator<IndexType> > vertices;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;

  {
    const Vec3i nrPoints=vel.getNrPoint();
    for(sizeType xx=0; xx<nrPoints.x(); xx++)
      for(sizeType yy=0; yy<nrPoints.y(); yy++) {
        IndexType pt=vel.getPt(Vec3i(xx,yy,0));
        IndexType velDir=vel.sampleSafe2D(pt);
        if(velDir.norm() < ScalarUtil<T>::scalar_eps())
          continue;
        if(len > 0.0f) {
          velDir.normalize();
          velDir*=len;
        }

        indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
        vertices.push_back(pt);
        vertices.push_back(pt+velDir);
      }
  }

  VTKWriter<TI> writer("Vec",path,true);
  writer.appendPoints(vertices.begin(),vertices.end());
  writer.appendCells(indices.begin(),indices.end(),VTKWriter<TI>::LINE);
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DMACGridVTK(const std::string& path,const MACGridType& mac)
{
  ASSERT(mac.getDim() == 2)

  const IndexType cellSz=mac.getCellSize();
  const IndexType offX(cellSz.x(),0.0f,0.0f);
  const IndexType offY(0.0f,cellSz.y(),0.0f);

  std::vector<IndexType,Eigen::aligned_allocator<IndexType> > vertices;
  std::vector<TI> colorData;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;

  {
    const GridType& uu=mac.getGu();
    const Vec3i nrPoints=mac.getGu().getNrPoint();
    for(sizeType xx=0; xx<nrPoints.x(); xx++)
      for(sizeType yy=0; yy<nrPoints.y(); yy++) {
        IndexType pt=uu.getPt(Vec3i(xx,yy,0))-offY*0.5f;
        indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
        vertices.push_back(pt);
        vertices.push_back(pt+offY);
        colorData.push_back(uu.get(Vec3i(xx,yy,0)));
        colorData.push_back(uu.get(Vec3i(xx,yy,0)));
      }
  }

  {
    const GridType& vv=mac.getGv();
    const Vec3i nrPoints=mac.getGv().getNrPoint();
    for(sizeType xx=0; xx<nrPoints.x(); xx++)
      for(sizeType yy=0; yy<nrPoints.y(); yy++) {
        IndexType pt=vv.getPt(Vec3i(xx,yy,0))-offX*0.5f;
        indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
        vertices.push_back(pt);
        vertices.push_back(pt+offX);
        colorData.push_back(vv.get(Vec3i(xx,yy,0)));
        colorData.push_back(vv.get(Vec3i(xx,yy,0)));
      }
  }

  VTKWriter<TI> writer("Vel",path,true);
  writer.appendPoints(vertices.begin(),vertices.end());
  writer.appendCells(indices.begin(),indices.end(),VTKWriter<typename IndexType::Scalar>::LINE);
  writer.appendCustomPointData("gridData",colorData.begin(),colorData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DVelocityGridVTK(const std::string& path,const MACGridType& vel,const TI& len)
{
  ASSERT(vel.getDim() == 2)

  std::vector<IndexType,Eigen::aligned_allocator<IndexType> > vertices;
  std::vector<TI> colorData;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indicesPt;

  {
    const GridType& uu=vel.getGu();
    const Vec3i nrPoints=vel.getGu().getNrPoint();
    for(sizeType xx=0; xx<nrPoints.x(); xx++)
      for(sizeType yy=0; yy<nrPoints.y(); yy++) {
        IndexType pt=uu.getPt(Vec3i(xx,yy,0));
        IndexType velDir=vel.sampleSafe2D(pt).template cast<typename IndexType::Scalar>();
        if(velDir.norm() < ScalarUtil<T>::scalar_eps())
          continue;
        if(len > 0.0f) {
          velDir.normalize();
          velDir*=len;
        }

        indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
        indicesPt.push_back(Vec2i(vertices.size(),0));
        vertices.push_back(pt);
        vertices.push_back(pt+velDir);
        colorData.push_back(0.0f);
        colorData.push_back(1.0f);
      }
  }

  {
    const GridType& vv=vel.getGv();
    const Vec3i nrPoints=vel.getGv().getNrPoint();
    for(sizeType xx=0; xx<nrPoints.x(); xx++)
      for(sizeType yy=0; yy<nrPoints.y(); yy++) {
        IndexType pt=vv.getPt(Vec3i(xx,yy,0));
        IndexType velDir=vel.sampleSafe2D(pt).template cast<typename IndexType::Scalar>();
        if(velDir.norm() < ScalarUtil<T>::scalar_eps())
          continue;
        if(len > 0.0f) {
          velDir.normalize();
          velDir*=len;
        }

        indices.push_back(Vec2i(vertices.size(),vertices.size()+1));
        indicesPt.push_back(Vec2i(vertices.size(),0));
        vertices.push_back(pt);
        vertices.push_back(pt+velDir);
        colorData.push_back(0.0f);
        colorData.push_back(1.0f);
      }
  }

  VTKWriter<TI> writer("Vel",path,true);
  writer.appendPoints(vertices.begin(),vertices.end());
  writer.appendCells(indices.begin(),indices.end(),VTKWriter<TI>::LINE);
  writer.appendCells(indicesPt.begin(),indicesPt.end(),VTKWriter<TI>::POINT);
  writer.appendCustomPointData("Color",colorData.begin(),colorData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write2DVelocityGridVTK(const std::string& path,const MACGridType& vel,const TI& time,const IndexType& sample,bool jitter)
{
  std::vector<IndexType,Eigen::aligned_allocator<IndexType> > vertices;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;
  std::vector<TI> colorData;

  for(T x=vel.getBB()._minC.x(); x<vel.getBB()._maxC.x(); x+=sample.x())
    for(T y=vel.getBB()._minC.y(); y<vel.getBB()._maxC.y(); y+=sample.y()) {
      IndexType pt(x,y,0);
      if(jitter) {
        pt.x()+=(rand()/(TI)RAND_MAX-0.5f)*sample.x()*0.5f;
        pt.y()+=(rand()/(TI)RAND_MAX-0.5f)*sample.y()*0.5f;
      }
      advectRK2(pt,time,vel,vertices,indices,colorData);
    }

  VTKWriter<TI> writer("Vel",path,true);
  writer.appendPoints(vertices.begin(),vertices.end());
  writer.appendCells(indices.begin(),indices.end(),VTKWriter<TI>::LINE);
  writer.appendCustomPointData("Color",colorData.begin(),colorData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::write3DVelocityGridVTK(const std::string& path,const MACGridType& vel,const TI& time,const IndexType& sample,bool jitter)
{
  std::vector<IndexType,Eigen::aligned_allocator<IndexType> > vertices;
  std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> > indices;
  std::vector<TI> colorData;

  for(T x=vel.getBB()._minC.x(); x<vel.getBB()._maxC.x(); x+=sample.x())
    for(T y=vel.getBB()._minC.y(); y<vel.getBB()._maxC.y(); y+=sample.y())
      for(T z=vel.getBB()._minC.z(); z<vel.getBB()._maxC.z(); z+=sample.z()) {
        IndexType pt(x,y,z);
        if(jitter) {
          pt.x()+=(rand()/(TI)RAND_MAX-0.5f)*sample.x()*0.5f;
          pt.y()+=(rand()/(TI)RAND_MAX-0.5f)*sample.y()*0.5f;
          pt.z()+=(rand()/(TI)RAND_MAX-0.5f)*sample.z()*0.5f;
        }
        advectRK2(pt,time,vel,vertices,indices,colorData);
      }

  VTKWriter<TI> writer("Vel",path,true);
  writer.appendPoints(vertices.begin(),vertices.end());
  writer.appendCells(indices.begin(),indices.end(),VTKWriter<TI>::LINE);
  writer.appendCustomPointData("Color",colorData.begin(),colorData.end());
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::advectRK2(const IndexType& from,const TI& time,const MACGridType& vel,
                                std::vector<IndexType,Eigen::aligned_allocator<IndexType> >& vertices,
                                std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> >& indices,
                                std::vector<TI>& colorData)
{
#define SUB_D 2
  IndexType curr=from;
  sizeType nrV=(sizeType)vertices.size();
  vertices.push_back(curr);
  //add vertices
  for(TI t=0; t<time;) {
    //determine step size
    ValueType v=vel.sampleSafe(curr);
    TI stepSz=(time-t),vNorm=v.norm();
    TI dist=vel.getCellSize().maxCoeff()/(scalar)SUB_D;
    if(stepSz*vNorm > dist)
      stepSz=dist/std::max<TI>(vNorm,ScalarUtil<TI>::scalar_eps());
    //RK2
    IndexType mid=curr+IndexType(v.x(),v.y(),v.z())*stepSz*0.5f;
    v=vel.sampleSafe(mid);
    curr+=IndexType(v.x(),v.y(),v.z())*stepSz;
    //add vertex
    indices.push_back(Vec2i(vertices.size()-1,vertices.size()));
    vertices.push_back(curr);
    //update time
    t+=stepSz;
  }
  //add colors
  sizeType nrVN=(sizeType)vertices.size();
  for(sizeType v=nrV; v<nrVN; v++)
    colorData.push_back((scalar)(v-nrV)/(scalar)(nrVN-1-nrV));
#undef SUB_D
}
//miscellaneous
template <typename T,typename TI,typename T2,typename TV>
struct ImplicitFuncTraits {
  typedef typename GridOp<T,TI,TV>::GridType GridType;
  void begin(GridType& g,ImplicitFunc<T2>& func) {}
  void end(GridType& g,ImplicitFunc<T2>& func) {}
};
template <typename T,typename TV>
struct ImplicitFuncTraits<T,T,T,TV> {
  typedef typename GridOp<T,T,TV>::GridType GridType;
  void begin(GridType& g,ImplicitFunc<T>& func) {
    func.beginSampleSet(g);
  }
  void end(GridType& g,ImplicitFunc<T>& func) {
    func.endSampleSet(g);
  }
};
template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::copyVelFromFunc(MACGridType& vel,const VelFunc<T2>& func,bool add)
{
  if(vel.getDim() >= 1)
    copyVelFromFunc(vel.getGu(),func,0,add);
  if(vel.getDim() >= 2)
    copyVelFromFunc(vel.getGv(),func,1,add);
  if(vel.getDim() >= 3)
    copyVelFromFunc(vel.getGw(),func,2,add);
}
template <typename T,typename TI,typename TV>
template<typename T2>
void GridOp<T,TI,TV>::copyVelFromFunc(GridType& vel,const VelFunc<T2>& func,const sizeType& a,bool add)
{
  if(!add)
    vel.init(EigenTraits<T2>::value());
  const Vec3i nrPoint=vel.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++)
      for(sizeType z=0; z<nrPoint.z(); z++)
        vel.get(Vec3i(x,y,z))+=func(vel.getPt(Vec3i(x,y,z)).template cast<T2>())[a];
}
template <typename T,typename TI,typename TV>
template <typename T2>
void GridOp<T,TI,TV>::copyFromImplictFuncCached(GridType& to,ImplicitFunc<T2>& func,bool add)
{
  ImplicitFuncTraits<T,TI,T2,TV>().begin(to,func);
  typedef typename ScalarUtil<T2>::ScalarVec3 Vec3;
  if(!add)
    to.init(EigenTraits<T>::value());
  const Vec3i nrPoint=to.getNrPoint();
  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<nrPoint.x(); x++) {
    //INFOV("%lu",x)
    for(sizeType y=0; y<nrPoint.y(); y++)
      for(sizeType z=0; z<nrPoint.z(); z++) {
        IndexType pt=to.getPt(Vec3i(x,y,z));
        to.get(Vec3i(x,y,z))+=(T)func(Vec3((T2)pt.x(),(T2)pt.y(),(T2)pt.z()));
      }
  }
  ImplicitFuncTraits<T,TI,T2,TV>().end(to,func);
}
template <typename T,typename TI,typename TV>
template <typename T2>
void GridOp<T,TI,TV>::copyFromImplictFunc(GridType& to,const ImplicitFunc<T2>& func,bool add)
{
  typedef typename ScalarUtil<T2>::ScalarVec3 Vec3;
  if(!add)
    to.init(EigenTraits<T>::value());
  const Vec3i nrPoint=to.getNrPoint();
  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<nrPoint.x(); x++) {
    //INFOV("%lu",x)
    for(sizeType y=0; y<nrPoint.y(); y++)
      for(sizeType z=0; z<nrPoint.z(); z++) {
        IndexType pt=to.getPt(Vec3i(x,y,z));
        to.get(Vec3i(x,y,z))+=(T)func(Vec3((T2)pt.x(),(T2)pt.y(),(T2)pt.z()));
      }
  }
}
template <typename T,typename TI,typename TV>
template<typename T2,typename TI2>
void GridOp<T,TI,TV>::copyFromOtherGrid(GridType& to,const Grid<T2,TI2>& from,const Vec3i* warpMode,bool add)
{
  if(!add)
    to.init(EigenTraits<T>::value());
  const Vec3i nrPoint=to.getNrPoint();
  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<nrPoint.x(); x++) {
    for(sizeType y=0; y<nrPoint.y(); y++)
      for(sizeType z=0; z<nrPoint.z(); z++)
        to.get(Vec3i(x,y,z))+=(T)from.sampleSafe(to.getPt(Vec3i(x,y,z)),warpMode);
  }
}
template <typename T,typename TI,typename TV>
template<typename T2,typename TI2>
void GridOp<T,TI,TV>::copyVelFromOtherGrid(MACGridType& to,const MACGrid<T2,TI2>& from,const Vec3c* warpMode,bool add)
{
  for(sizeType d=0; d<from.getDim(); d++) {
    if(!add)
      to.getComp(d).init(EigenTraits<T>::value());
    const Vec3i& nrPoint=to.getComp(d).getNrPoint();
    OMP_PARALLEL_FOR_
    for(sizeType x=0; x<nrPoint.x(); x++) {
      for(sizeType y=0; y<nrPoint.y(); y++)
        for(sizeType z=0; z<nrPoint.z(); z++)
          to.getComp(d).get(Vec3i(x,y,z))+=(T)from.sampleSafe(to.getComp(d).getPt(Vec3i(x,y,z)),warpMode)[d];
    }
  }
}
template <typename T,typename TI,typename TV>
template<typename T2,typename TI2>
void GridOp<T,TI,TV>::copyFromOtherGridOfSameGeometry(GridType& to,const Grid<T2,TI2>& from,bool add)
{
  if(!add)
    to.init(EigenTraits<T>::value());
  const Vec3i nrPoint=to.getNrPoint();
  ASSERT(nrPoint==from.getNrPoint());
  OMP_PARALLEL_FOR_
  for(sizeType x=0; x<nrPoint.x(); x++) {
    for(sizeType y=0; y<nrPoint.y(); y++)
      for(sizeType z=0; z<nrPoint.z(); z++)
        to.get(Vec3i(x,y,z))+=(T)from.get(Vec3i(x,y,z));
  }
}
template <typename T,typename TI,typename TV>
template<typename T2,typename TI2>
void GridOp<T,TI,TV>::copyFromOtherGridOfSameGeometry(MACGridType& to,const MACGrid<T2,TI2>& from,bool add)
{
  if(to.getDim() >= 1)
    copyFromOtherGridOfSameGeometry(to.getGu(),from.getGu(),add);
  if(to.getDim() >= 2)
    copyFromOtherGridOfSameGeometry(to.getGv(),from.getGv(),add);
  if(to.getDim() >= 3)
    copyFromOtherGridOfSameGeometry(to.getGw(),from.getGw(),add);
}
template <typename T,typename TI,typename TV>
template <typename COMPARE,typename ALLOC>
void GridOp<T,TI,TV>::getAllDifferentValues(const GridType& g,std::set<T,COMPARE,ALLOC>& vals)
{
  const Vec3i nrPoint=g.getNrPoint();
  for(sizeType x=0; x<nrPoint.x(); x++)
    for(sizeType y=0; y<nrPoint.y(); y++)
      for(sizeType z=0; z<nrPoint.z(); z++)
        vals.insert(g.get(Vec3i(x,y,z)));
}
template <typename T,typename TI,typename TV>
template <typename COMPARE,typename ALLOC>
void GridOp<T,TI,TV>::getAllDifferentValues(const MACGridType& g,std::set<T,COMPARE,ALLOC>& vals)
{
  if(g.getDim() >= 1)
    getAllDifferentValues(g.getGu(),vals);
  if(g.getDim() >= 2)
    getAllDifferentValues(g.getGv(),vals);
  if(g.getDim() >= 3)
    getAllDifferentValues(g.getGw(),vals);
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromFaceToCenter(const MACGridType& from,Grid<ValueType,T>& to)
{
  ASSERT(from.getNrCell() == to.getNrCell())
  ASSERT(to.isCenter())
  if(to.getDim() == 2)
    fromFaceToCenter2D(from,to);
  else
    fromFaceToCenter3D(from,to);
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromFaceToCenter3D(const MACGridType& from,Grid<ValueType,T>& to)
{
  const Vec3i nrCell=to.getNrCell();
  for(sizeType x=0; x<nrCell.x(); x++)
    for(sizeType y=0; y<nrCell.y(); y++)
      for(sizeType z=0; z<nrCell.z(); z++) {
        ValueType& val=to.get(Vec3i(x,y,z));
        val(0)=(from.getGu().get(Vec3i(x,y,z))+from.getGu().get(Vec3i(x+1,y,z)))*0.5f;
        val(1)=(from.getGv().get(Vec3i(x,y,z))+from.getGv().get(Vec3i(x,y+1,z)))*0.5f;
        val(2)=(from.getGw().get(Vec3i(x,y,z))+from.getGw().get(Vec3i(x,y,z+1)))*0.5f;
      }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromFaceToCenter2D(const MACGridType& from,Grid<ValueType,T>& to)
{
  const Vec3i nrCell=to.getNrCell();
  for(sizeType x=0; x<nrCell.x(); x++)
    for(sizeType y=0; y<nrCell.y(); y++) {
      ValueType& val=to.get(Vec3i(x,y,0));
      val(0)=(from.getGu().get(Vec3i(x,y,0))+from.getGu().get(Vec3i(x+1,y,0)))*0.5f;
      val(1)=(from.getGv().get(Vec3i(x,y,0))+from.getGv().get(Vec3i(x,y+1,0)))*0.5f;
      val(2)=0.0f;
    }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromCenterToFace(const Grid<ValueType,T>& from,MACGridType& to)
{
  ASSERT(from.getNrCell() == to.getNrCell())
  ASSERT(from.isCenter())
  if(to.getDim() == 2)
    fromCenterToFace2D(from,to);
  else
    fromCenterToFace3D(from,to);
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromCenterToFace3D(const Grid<ValueType,T>& from,MACGridType& to)
{
  const Vec3i nrCell=to.getNrCell();
  for(sizeType x=0; x<nrCell.x(); x++)
    for(sizeType y=0; y<nrCell.y(); y++)
      for(sizeType z=0; z<nrCell.z(); z++) {
        //X
        if(x == 0)
          to.getGu().get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z)).x();
        else
          to.getGu().get(Vec3i(x,y,z))=(from.get(Vec3i(x,y,z)).x()+from.get(Vec3i(x-1,y,z)).x())*0.5f;
        if(x == nrCell.x()-1)
          to.getGu().get(Vec3i(x+1,y,z))=from.get(Vec3i(x,y,z)).x();
        else
          to.getGu().get(Vec3i(x+1,y,z))=(from.get(Vec3i(x,y,z)).x()+from.get(Vec3i(x+1,y,z)).x())*0.5f;

        //Y
        if(y == 0)
          to.getGv().get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z)).y();
        else
          to.getGv().get(Vec3i(x,y,z))=(from.get(Vec3i(x,y,z)).y()+from.get(Vec3i(x,y-1,z)).y())*0.5f;
        if(y == nrCell.y()-1)
          to.getGv().get(Vec3i(x,y+1,z))=from.get(Vec3i(x,y,z)).y();
        else
          to.getGv().get(Vec3i(x,y+1,z))=(from.get(Vec3i(x,y,z)).y()+from.get(Vec3i(x,y+1,z)).y())*0.5f;

        //Z
        if(z == 0)
          to.getGw().get(Vec3i(x,y,z))=from.get(Vec3i(x,y,z)).z();
        else
          to.getGw().get(Vec3i(x,y,z))=(from.get(Vec3i(x,y,z)).z()+from.get(Vec3i(x,y-1,z)).z())*0.5f;
        if(z == nrCell.z()-1)
          to.getGw().get(Vec3i(x,y+1,z))=from.get(Vec3i(x,y,z)).z();
        else
          to.getGw().get(Vec3i(x,y+1,z))=(from.get(Vec3i(x,y,z)).z()+from.get(Vec3i(x,y+1,z)).z())*0.5f;
      }
}
template <typename T,typename TI,typename TV>
void GridOp<T,TI,TV>::fromCenterToFace2D(const Grid<ValueType,T>& from,MACGridType& to)
{
  const Vec3i nrCell=to.getNrCell();
  for(sizeType x=0; x<nrCell.x(); x++)
    for(sizeType y=0; y<nrCell.y(); y++) {
      //X
      if(x == 0)
        to.getGu().get(Vec3i(x,y,0))=from.get(Vec3i(x,y,0)).x();
      else
        to.getGu().get(Vec3i(x,y,0))=(from.get(Vec3i(x,y,0)).x()+from.get(Vec3i(x-1,y,0)).x())*0.5f;
      if(x == nrCell.x()-1)
        to.getGu().get(Vec3i(x+1,y,0))=from.get(Vec3i(x,y,0)).x();
      else
        to.getGu().get(Vec3i(x+1,y,0))=(from.get(Vec3i(x,y,0)).x()+from.get(Vec3i(x+1,y,0)).x())*0.5f;

      //Y
      if(y == 0)
        to.getGv().get(Vec3i(x,y,0))=from.get(Vec3i(x,y,0)).y();
      else
        to.getGv().get(Vec3i(x,y,0))=(from.get(Vec3i(x,y,0)).y()+from.get(Vec3i(x,y-1,0)).y())*0.5f;
      if(y == nrCell.y()-1)
        to.getGv().get(Vec3i(x,y+1,0))=from.get(Vec3i(x,y,0)).y();
      else
        to.getGv().get(Vec3i(x,y+1,0))=(from.get(Vec3i(x,y,0)).y()+from.get(Vec3i(x,y+1,0)).y())*0.5f;
    }
}

#define DEF_GRID(T,TI)  \
template class GridOp<T,TI>;  \
template bool GridOp<T,TI>::solveEikonalFM<T>(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T,TI>* extra,T thres,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary);  \
\
template void GridOp<T,TI>::copyVelFromFunc<sizeType>(MACGridType& vel,const VelFunc<sizeType>& func,bool add);  \
template void GridOp<T,TI>::copyVelFromFunc<unsigned char>(MACGridType& vel,const VelFunc<unsigned char>& func,bool add);  \
template void GridOp<T,TI>::copyVelFromFunc<scalarF>(MACGridType& vel,const VelFunc<scalarF>& func,bool add);  \
template void GridOp<T,TI>::copyVelFromFunc<scalarD>(MACGridType& vel,const VelFunc<scalarD>& func,bool add);  \
\
template void GridOp<T,TI>::copyVelFromFunc<sizeType>(GridType& vel,const VelFunc<sizeType>& func,const sizeType& a,bool add); \
template void GridOp<T,TI>::copyVelFromFunc<unsigned char>(GridType& vel,const VelFunc<unsigned char>& func,const sizeType& a,bool add); \
template void GridOp<T,TI>::copyVelFromFunc<scalarF>(GridType& vel,const VelFunc<scalarF>& func,const sizeType& a,bool add); \
template void GridOp<T,TI>::copyVelFromFunc<scalarD>(GridType& vel,const VelFunc<scalarD>& func,const sizeType& a,bool add); \
\
template void GridOp<T,TI>::copyFromImplictFuncCached<sizeType>(GridType& to,ImplicitFunc<sizeType>& func,bool add);  \
template void GridOp<T,TI>::copyFromImplictFuncCached<unsigned char>(GridType& to,ImplicitFunc<unsigned char>& func,bool add);  \
template void GridOp<T,TI>::copyFromImplictFuncCached<scalarF>(GridType& to,ImplicitFunc<scalarF>& func,bool add);  \
template void GridOp<T,TI>::copyFromImplictFuncCached<scalarD>(GridType& to,ImplicitFunc<scalarD>& func,bool add);  \
\
template void GridOp<T,TI>::copyFromImplictFunc<sizeType>(GridType& to,const ImplicitFunc<sizeType>& func,bool add);  \
template void GridOp<T,TI>::copyFromImplictFunc<unsigned char>(GridType& to,const ImplicitFunc<unsigned char>& func,bool add);  \
template void GridOp<T,TI>::copyFromImplictFunc<scalarF>(GridType& to,const ImplicitFunc<scalarF>& func,bool add);  \
template void GridOp<T,TI>::copyFromImplictFunc<scalarD>(GridType& to,const ImplicitFunc<scalarD>& func,bool add);

//instance
DEF_GRID(sizeType,scalarF)
DEF_GRID(unsigned char,scalarF)
DEF_GRID(scalarF,scalarF)
DEF_GRID(scalarD,scalarF)

DEF_GRID(sizeType,scalarD)
DEF_GRID(unsigned char,scalarD)
DEF_GRID(scalarF,scalarD)
DEF_GRID(scalarD,scalarD)

PRJ_END
