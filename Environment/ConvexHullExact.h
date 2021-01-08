#ifndef CONVEX_HULL_EXACT_H
#define CONVEX_HULL_EXACT_H

#include "ObjMeshGeomCellExact.h"
#include "EdgeExact.h"
extern "C"
{
#include "openGJK.h"
}

PRJ_BEGIN

struct ALIGN_16 ConvexHullExact : public ObjMeshGeomCellExact
{
  ConvexHullExact();
  ConvexHullExact(const ObjMeshGeomCell& exact);
  virtual ~ConvexHullExact();
  virtual bool read(std::istream& is,IOData* dat) override;
  virtual bool write(std::ostream& os,IOData* dat) const override;
  virtual std::shared_ptr<SerializableBase> copy() const override;
  virtual std::string type() const override;
  void parityCheck() const;
  bool closest(const PT& pt,PT& n,PT& normal,MAT3& hessian,Vec2i& feat,bool cache=false,std::vector<PT,Eigen::aligned_allocator<PT>>* history=NULL) const override;
  template <typename T2>
  T2 closestGJK(const Eigen::Matrix<T2,3,1>& pt,Eigen::Matrix<T2,3,1>& n,Eigen::Matrix<T2,3,1>& normal) const
  {
    double* coord0=new double[3];

    bd bd0;
    bd0.numpoints=1;
    bd0.coord=&coord0;
    for(sizeType d=0; d<3; d++)
      bd0.s[d]=coord0[d]=std::to_double(pt[d]);

    bd bd1;
    bd1.numpoints=_bd.numpoints;
    bd1.coord=_bd.coord;

    simplex s;
    T2 d=gjk(bd0,bd1,&s);

    n.setZero();
    for(sizeType i=0; i<s.nvrtx; i++)
      for(sizeType d=0; d<3; d++)
        n[d]+=s.lambdas[i]*(s.q[i][d]-s.p[i][d]);
    normal=-n;
    normal/=std::sqrt(n.squaredNorm());

    delete [] coord0;
    return d;
  }
  void writeHistoryDistVTK(const std::string& path,sizeType res=10) const;
  void compare(const ObjMeshGeomCellExact& mExact,sizeType res=10) const;
private:
  void buildBD();
  void removeBD();
  ALIGN_16 std::vector<EdgeExact> _ess;
  ALIGN_16 std::vector<std::vector<sizeType>> _eNss;
  bd _bd;
};

PRJ_END

#endif
