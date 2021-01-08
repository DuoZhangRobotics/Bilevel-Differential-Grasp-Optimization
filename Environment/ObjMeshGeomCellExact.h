#ifndef OBJ_MESH_GEOM_CELL_EXACT_H
#define OBJ_MESH_GEOM_CELL_EXACT_H

#include <CommonFile/geom/ObjMeshGeomCell.h>
#include "TriangleExact.h"
#include "BBoxExact.h"
#include "MPQZIO.h"

PRJ_BEGIN

struct ALIGN_16 ObjMeshGeomCellExact : public SerializableBase
{
  typedef mpq_class T;
  typedef Eigen::Matrix<T,3,1> PT;
  typedef Eigen::Matrix<T,3,3> MAT3;
  ObjMeshGeomCellExact();
  ObjMeshGeomCellExact(const ObjMeshGeomCell& exact);
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  const BBoxExact& getBB() const;
  bool empty() const;
  virtual bool closest(const PT& pt,PT& n,PT& normal,MAT3& hessian,Vec2i& feat,bool cache=false,std::vector<PT,Eigen::aligned_allocator<PT>>* history=NULL) const;
  template <typename T2>
  T2 closest(const Eigen::Matrix<T2,3,1>& pt,Eigen::Matrix<T2,3,1>& n,Eigen::Matrix<T2,3,1>& normal,Eigen::Matrix<T2,3,3>& hessian,Vec2i& feat,bool cache=false,std::vector<Eigen::Matrix<T2,3,1>,Eigen::aligned_allocator<Eigen::Matrix<T2,3,1>>>* history=NULL) const
  {
    bool inside;
    MAT3 hessianR;
    PT ptR=pt.unaryExpr([&](const T2& t) {
      T ret;
      castRational(t,ret);
      return ret;
    }),nR,normalR;
    std::vector<PT,Eigen::aligned_allocator<PT>> historyT2;
    inside=closest(ptR,nR,normalR,hessianR,feat,cache,history?&historyT2:NULL);
    if(history) {
      history->resize(historyT2.size());
      for(sizeType i=0; i<(sizeType)history->size(); i++)
        history->at(i)=castRational<Eigen::Matrix<T2,3,1>,PT>(historyT2[i]);
    }
    inside=closest(ptR,nR,normalR,hessianR,feat,cache,history?&historyT2:NULL);
    //cast
    n=castRational<Eigen::Matrix<T2,3,1>,PT>(nR);
    normal=castRational<Eigen::Matrix<T2,3,1>,PT>(normalR);
    hessian=castRational<Eigen::Matrix<T2,3,3>,MAT3>(hessianR);
    //post process
    T2 nLen=std::sqrt(n.squaredNorm());
    hessian/=std::max(std::numeric_limits<T2>::epsilon(),nLen);
    if(nR[0]==0 && nR[1]==0 && nR[2]==0) {
      nLen=0;
      normal/=std::max(std::numeric_limits<T2>::epsilon(),std::sqrt(normal.squaredNorm()));
    } else {
      normal=n/std::max(std::numeric_limits<T2>::epsilon(),nLen);
      if(inside)
        nLen*=-1;
      else {
        normal*=-1;
        hessian*=-1;
      }
    }
    return nLen;
  }
  virtual void writePointDistVTK(const std::string& path,sizeType res=10) const;
  virtual void debugPointDist(sizeType nrIter=100) const;
  ObjMesh getMesh() const;
protected:
  ALIGN_16 std::vector<TriangleExact> _tss;
  ALIGN_16 std::vector<PT,Eigen::aligned_allocator<PT>> _vss;
  ALIGN_16 std::vector<Vec3i,Eigen::aligned_allocator<Vec3i>> _iss;
  ALIGN_16 std::vector<Node<sizeType,BBoxExact>> _bvh;
};

PRJ_END

#endif
