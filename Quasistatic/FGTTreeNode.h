#ifndef FGT_TREE_NODE_H
#define FGT_TREE_NODE_H

#include <CommonFile/MathBasic.h>
#include <Utils/SparseUtils.h>

PRJ_BEGIN

#include <Utils/MapTypePragma.h>
template <typename T>
struct FGTTreeNode
{
public:
  DECL_MAP_FUNCS
  DECL_MAP_TYPES_T
  FGTTreeNode();
  FGTTreeNode(Vec* Sy,Mat3XT& xy,const Vec2i& range,sizeType leafThres,Mat3XT* xyn=NULL);
  void parityCheck(const Vec* Sy,const Mat3XT& xy,T thres=1e-6f) const;
  void transform(const Mat3T& R,const Vec3T& t,const Mat3XT& xy,bool xyTransformed);
  void transform(const Mat3T& R,const Vec3T& t);
  Vec3T variance(Mat3XT& y) const;
  Sphere<T> computeSphere(const Mat3XT& y,const Vec3T& ctr) const;
  Sphere<T> mergeSphere(const Sphere<T>& l,const Sphere<T>& r) const;
  T distTo(const FGTTreeNode<T>& other) const;
  sizeType size() const;
  //FGT
  static void closestYNode(const FGTTreeNode<T>** minLeaf,T& minDist,const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode);
  static void initErrorBound(const Vec* Sy,const Mat3XT& y,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr);
  static void FGT(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr,T eps,Vec3i* profile=NULL);
  static void contrib(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T F,T invHSqr,T eps,Vec3i* profile=NULL);
  static sizeType cost(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,sizeType& p,T invHSqr,T errBound);
  static sizeType costChildren(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,T invHSqr,T errBound);
  //mean
  static void mean(Vec& G,MatX4T* DGDT,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,const Vec2T& minMax,T F,T invHSqr,T errMean,T eps);
  static T errorMean(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,Vec2T& minMax,T invHSqr);
  static sizeType costMean(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode);
  //mean-ctr
  static void meanCtr(Vec& G,MatX4T* DGDT,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T F,T invHSqr,T errMean,T eps);
  static T errorMeanCtr(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,Vec2T& minMax,T invHSqr);
  static sizeType costMeanCtr(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode);
  //direct
  static void direct(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr);
  static sizeType costDirect(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode);
  //taylor
  static void taylor(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr,T eps);
  static void taylor(Vec& G,MatX4T* DGDT,const Vec* Sy,const Mat3XT& y,const Mat3XT* yl,const Mat3XT& x,const FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,T invHSqr,sizeType p);
  static sizeType costTaylor(const FGTTreeNode<T>& yNode,const FGTTreeNode<T>& xNode,sizeType& p,T invHSqr,T errBound);
  //helper
  static void debugSwapId(sizeType base,sizeType N);
  static void swapId(std::vector<sizeType> id,std::function<void(sizeType,sizeType)> f);
  static void randomTree(Vec* Sy,Mat3XT& y,Mat3XT& yl,Mat3XT& x,FGTTreeNode<T>& yNode,FGTTreeNode<T>& xNode,sizeType N,sizeType leafThres,bool random);
  static void debugTaylor(sizeType N,sizeType leafThres,bool random,T invHSqr,T eps,bool useSy);
  static void debugFGT(sizeType N,sizeType leafThres,bool random,T invHSqr,T eps,bool useSy);
  static void debugTree(sizeType N,sizeType leafThres,bool random,bool useSy);
  static sizeType combination(sizeType n,sizeType p);
  static sizeType factorial(sizeType p);
private:
  //tmp
  T _tildeGMinInit;
  T _tildeGMin;
  T _Fs,_FtSave;
  //data
  std::shared_ptr<FGTTreeNode> _l,_r;
  Sphere<T> _spherel,_sphere;
  BBox<T> _bbl,_bb;
  Vec2i _range;
};

PRJ_END

#endif
