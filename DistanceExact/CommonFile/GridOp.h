#ifndef GRID_OP_H
#define GRID_OP_H

#include "GridBasic.h"
#include "ImplicitFuncInterface.h"
#include "ObjMesh.h"
#include <set>

PRJ_BEGIN

template <typename T,typename TI,typename TV=std::vector<T,Eigen::aligned_allocator<T> > >
class GridOp
{
public:
  enum FM_STATE {
    UNKNOWN = 0,
    KNOWN = 1,
    CLOSE = 2,
    BOUNDARY = 4,
    POSITIVE = 16,
    NEGATIVE = 32,
  };
  typedef Grid<T,TI,TV> GridType;
  typedef MACGrid<T,TI,TV> MACGridType;
  typedef Grid<unsigned char,TI> TagGridType;
  typedef Grid<sizeType,TI> MarkGridType;
  typedef typename GridType::ValueType ValueType;
  typedef typename GridType::IndexType IndexType;
  //pde
  static void solveEikonalPDE(GridType& from,const T& dist,bool preserveSign);
  static void smoothEikonalPDE3DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr);
  static void smoothEikonalPDE2DLayer(const sizeType& x,const GridType& from,GridType& smooth,const Vec3i& nrPoint,const T& maxCellSzSqr);
  static void solveEikonalPDE3DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt);
  static void solveEikonalPDE2DLayer(const sizeType& x,const GridType& from,const GridType& smooth,GridType& to,const Vec3i& nrPoint,const IndexType& invCellSz,const T& dt);
  static void solveEikonalPDE3D(GridType& from,GridType& to,GridType& smooth,const T& dist,bool preserveSign);
  static void solveEikonalPDE2D(GridType& from,GridType& to,GridType& smooth,const T& dist,bool preserveSign);
  static FORCE_INLINE T upwindingGrad(T partialPlus,T partialMinus,const T& phi0);
  //fast march
  template<typename VT2> static bool solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& normals,Grid<IndexType,TI>& normalExtra);
  template<typename VT2> static bool solveEikonalFMFromNode(Grid<T,TI>& from,const std::vector<VT2,Eigen::aligned_allocator<VT2> >& nodes);
  template<typename T2> static bool solveEikonalFM(GridType& from,const T& low,const T& high,Grid<T2,TI>* extra=NULL,T thres=(T)-1.0f,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary=NULL) {
    std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> > known;
    const Vec3i nrPoint=from.getNrPoint();
    for(sizeType x=0; x<nrPoint.x(); x++)
      for(sizeType y=0; y<nrPoint.y(); y++)
        for(sizeType z=0; z<nrPoint.z(); z++) {
          const Vec3i id(x,y,z);
          if(from.get(id) > low && from.get(id) < high)
            known.push_back(id);
        }
    return solveEikonalFM(from,known,extra,thres,boundary);
  }
  template<typename T2> static bool solveEikonalFM(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra=NULL,T thres=(T)-1.0f,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary=NULL);
  template<typename T2> static bool solveEikonalFM3D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra=NULL,T thres=(T)-1.0f,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary=NULL);
  template<typename T2> static bool solveEikonalFM2D(GridType& from,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& known,Grid<T2,TI>* extra=NULL,T thres=(T)-1.0f,const std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >* boundary=NULL);
  template<typename T2> static void extrapolate3D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos);
  template<typename T2> static void extrapolate2D(const GridType& from,const TagGridType& tags,Grid<T2,TI>& extra,const Vec3i& pos);
  static void floodFill3D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags);
  static void floodFill2D(GridType& from,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& id,TagGridType& tags);
  static bool updateClose2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close);
  static bool updateClose3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& close);
  static bool updateKnown2D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id);
  static bool updateKnown3D(GridType& result,TagGridType& tags,MarkGridType& heapOffsets,std::vector<Vec3i,Eigen::aligned_allocator<Vec3i> >& heap,const Vec3i& id);
  //reinitialize
  static void reinitialize(GridType& from);
  //smooth
  static void smooth(GridType& from,const TagGridType* tagField=0,unsigned char tag=0);
  static void smooth3D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag);
  static void smooth2D(GridType& from,GridType& tmp,const TagGridType* tagField,unsigned char tag);
  //io
  static void write2DScalarGridObj(ObjMesh& obj,const GridType& grd,bool moveVertex=true);
  static void write2DScalarGridVTK(const std::string& path,const GridType& grd,bool moveVertex=true);
  static void write2DScalarBarChartVTK(const std::string& path,const GridType& grd,bool moveVertex=true,const TI& scale=1.0f);
  static void write2DMarkGridVTK(const std::string& path,const MarkGridType& grd,bool moveVertex=true);
  static void write2DScalarGridGradVTK(const std::string& path,const GridType& grd,const sizeType& sampleInterval,const TI& len);
  static void write2DScalarGridGradMagVTK(const std::string& path,const GridType& grd);
  static void write3DScalarGridVTK(const std::string& path,const GridType& grd);
  static void write2DVectorGridVTK(const std::string& path,const Grid<IndexType,TI>& vel,const TI& len=0.0f);
  static void write2DMACGridVTK(const std::string& path,const MACGridType& mac);
  static void write2DVelocityGridVTK(const std::string& path,const MACGridType& vel,const TI& len=0.0f);
  static void write2DVelocityGridVTK(const std::string& path,const MACGridType& vel,const TI& time,const IndexType& sample,bool jitter);
  static void write3DVelocityGridVTK(const std::string& path,const MACGridType& vel,const TI& time,const IndexType& sample,bool jitter);
  static void advectRK2(const IndexType& from,const TI& time,const MACGridType& vel,
                        std::vector<IndexType,Eigen::aligned_allocator<IndexType> >& vertices,
                        std::vector<Vec2i,Eigen::aligned_allocator<Vec2i> >& indices,
                        std::vector<TI>& colorData);
  //miscellaneous
  template<typename T2>
  static void copyVelFromFunc(MACGridType& vel,const VelFunc<T2>& func,bool add=false);
  template<typename T2>
  static void copyVelFromFunc(GridType& vel,const VelFunc<T2>& func,const sizeType& a,bool add=false);
  template<typename T2>
  static void copyFromImplictFuncCached(GridType& to,ImplicitFunc<T2>& func,bool add=false);
  template<typename T2>
  static void copyFromImplictFunc(GridType& to,const ImplicitFunc<T2>& func,bool add=false);
  template<typename T2,typename TI2>
  static void copyFromOtherGrid(GridType& to,const Grid<T2,TI2>& from,const Vec3i* warpMode=NULL,bool add=false);
  template<typename T2,typename TI2>
  static void copyVelFromOtherGrid(MACGridType& to,const MACGrid<T2,TI2>& from,const Vec3c* warpMode=NULL,bool add=false);
  template<typename T2,typename TI2>
  static void copyFromOtherGridOfSameGeometry(GridType& to,const Grid<T2,TI2>& from,bool add=false);
  template<typename T2,typename TI2>
  static void copyFromOtherGridOfSameGeometry(MACGridType& to,const MACGrid<T2,TI2>& from,bool add=false);
  template <typename COMPARE,typename ALLOC>
  static void getAllDifferentValues(const GridType& g,std::set<T,COMPARE,ALLOC>& vals);
  template <typename COMPARE,typename ALLOC>
  static void getAllDifferentValues(const MACGridType& g,std::set<T,COMPARE,ALLOC>& vals);
  //interpolation
  static void fromFaceToCenter(const MACGridType& from,Grid<ValueType,T>& to);
  static void fromFaceToCenter3D(const MACGridType& from,Grid<ValueType,T>& to);
  static void fromFaceToCenter2D(const MACGridType& from,Grid<ValueType,T>& to);
  static void fromCenterToFace(const Grid<ValueType,T>& from,MACGridType& to);
  static void fromCenterToFace3D(const Grid<ValueType,T>& from,MACGridType& to);
  static void fromCenterToFace2D(const Grid<ValueType,T>& from,MACGridType& to);
};

PRJ_END

#endif
