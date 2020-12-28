#ifndef ARTICULATED_BODY_PRAGMA_H
#define ARTICULATED_BODY_PRAGMA_H

//basic rigid body transformation
#define TRANSI(JSS,I) JSS.template block<3,4>(0,(I)*4)
#define TRANSI6(JSS,I) JSS.template block<6,6>(0,(I)*6)
#define GETT(TID,JSS,I) Eigen::Block<Mat3Xd,3,4> TID=TRANSI(JSS,I);
#define GETTC(TID,JSS,I) Eigen::Block<const Mat3Xd,3,4> TID=TRANSI(JSS,I);
#define GETTM(TID,JSS,I) Eigen::Block<Eigen::Map<Mat3Xd,0,Eigen::OuterStride<> >,3,4> TID=TRANSI(JSS,I);
#define GETTCM(TID,JSS,I) Eigen::Block<const Eigen::Map<Mat3Xd,0,Eigen::OuterStride<> >,3,4> TID=TRANSI(JSS,I);
#define GETT_T(TID,JSS,I) Eigen::Block<Mat3XT,3,4> TID=TRANSI(JSS,I);
#define GETTC_T(TID,JSS,I) Eigen::Block<const Mat3XT,3,4> TID=TRANSI(JSS,I);
#define GETTM_T(TID,JSS,I) Eigen::Block<Eigen::Map<Mat3XT,0,Eigen::OuterStride<> >,3,4> TID=TRANSI(JSS,I);
#define GETTCM_T(TID,JSS,I) Eigen::Block<const Eigen::Map<Mat3XT,0,Eigen::OuterStride<> >,3,4> TID=TRANSI(JSS,I);
#define ROT(A) (A).template block<3,3>(0,0)
#define CTR(A) (A).template block<3,1>(0,3)
#define ROTI(A,I) (A).template block<3,3>(0,(I)*4+0)
#define CTRI(A,I) (A).template block<3,1>(0,(I)*4+3)
#define DWDLI(A,I) (A).col(((I)<<1)+0)
#define DTDLI(A,I) (A).col(((I)<<1)+1)

//6x6 spatial block
#define SPATIAL_BLK00(A) A.template block<3,3>(0,0)
#define SPATIAL_BLK10(A) A.template block<3,3>(3,0)
#define SPATIAL_BLK01(A) A.template block<3,3>(0,3)
#define SPATIAL_BLK11(A) A.template block<3,3>(3,3)

//multiply transformation
#define APPLY_TRANS(C,A,B)  \
C=(ROT(A)*(B)).eval();  \
CTR(C)+=CTR(A);

//invert the transformation
#define INV(IA,A) \
ROT(IA)=ROT(A).transpose().eval();  \
CTR(IA)=-(ROT(IA)*CTR(A)).eval();

//degree to radius
#define D2R(X) (X*M_PI/180.0f)

#define DECL_MAP_TYPES_T  \
typedef typename SparseTraits<T>::DMat DMat;   \
typedef typename SparseTraits<T>::SMat SMat;   \
typedef typename SparseTraits<T>::STrip STrip;  \
typedef typename SparseTraits<T>::STrips STrips; \
typedef Eigen::Matrix<T,-1,1> Vec;  \
typedef Eigen::Matrix<T,2,1> Vec2T;  \
typedef Eigen::Matrix<T,3,1> Vec3T;  \
typedef Eigen::Matrix<T,4,1> Vec4T;  \
typedef Eigen::Matrix<T,-1,-1> MatT;  \
    \
typedef Eigen::Matrix<T,2,2> Mat2T;  \
typedef Eigen::Matrix<T,3,3> Mat3T;  \
typedef Eigen::Matrix<T,2,-1> Mat2XT;  \
typedef Eigen::Matrix<T,3,-1> Mat3XT;  \
typedef Eigen::Matrix<T,-1,3> MatX3T;  \
typedef Eigen::Matrix<T,4,4> Mat4T;  \
typedef Eigen::Matrix<T,4,-1> Mat4XT;  \
typedef Eigen::Matrix<T,-1,4> MatX4T;  \
typedef Eigen::Matrix<T,3,4> Mat3X4T;  \
typedef Eigen::Matrix<T,6,3> Mat6X3T;  \
typedef Eigen::Matrix<T,6,-1> Mat6XT;  \
typedef Eigen::Matrix<T,-1,6> MatX6T;  \
typedef Eigen::Matrix<T,6,6> Mat6T;  \
typedef Eigen::Matrix<T,6,1> Vec6T;  \
typedef Eigen::Matrix<T,9,9> Mat9T;  \
typedef Eigen::Matrix<T,9,12> Mat9X12T;  \
typedef Eigen::Matrix<T,12,-1> Mat12XT;  \
    \
typedef Eigen::Map<Vec> VecM;  \
typedef Eigen::Map<Coli> ColiM;  \
typedef Eigen::Map<const Coli> ColiCM;  \
typedef Eigen::Map<Vec,0,Eigen::InnerStride<> > VecMI;  \
typedef Eigen::Map<const Vec> VecCM;  \
typedef Eigen::Map<const Vec,0,Eigen::InnerStride<> > VecCMI;  \
typedef Eigen::Map<Vec3T> Vec3TM;  \
typedef Eigen::Map<const Vec3T> Vec3TCM;  \
typedef Eigen::Map<Vec4T> Vec4TM;  \
typedef Eigen::Map<const Vec4T> Vec4TCM;  \
typedef Eigen::Matrix<T,-1,-1,Eigen::RowMajor> MatTR; \
typedef Eigen::Map<MatTR,0,Eigen::OuterStride<> > MatTRM;  \
typedef Eigen::Map<const MatTR,0,Eigen::OuterStride<> > MatTRCM;  \
typedef Eigen::Map<MatT,0,Eigen::OuterStride<> > MatTM;  \
typedef Eigen::Map<const MatT,0,Eigen::OuterStride<> > MatTCM;  \
typedef Eigen::Map<Mat2XT,0,Eigen::OuterStride<> > Mat2XTM;  \
typedef Eigen::Map<const Mat2XT,0,Eigen::OuterStride<> > Mat2XTCM;  \
typedef Eigen::Map<Mat3XT,0,Eigen::OuterStride<> > Mat3XTM;  \
typedef Eigen::Map<const Mat3XT,0,Eigen::OuterStride<> > Mat3XTCM;  \
typedef Eigen::Map<Mat3X4T,0,Eigen::OuterStride<> > Mat3X4TM;  \
typedef Eigen::Map<const Mat3X4T,0,Eigen::OuterStride<> > Mat3X4TCM;  \
typedef Eigen::Map<MatX3T,0,Eigen::OuterStride<> > MatX3TM;  \
typedef Eigen::Map<const MatX3T,0,Eigen::OuterStride<> > MatX3TCM;  \
typedef Eigen::Map<Mat3T,0,Eigen::OuterStride<> > Mat3TM;  \
typedef Eigen::Map<const Mat3T,0,Eigen::OuterStride<> > Mat3TCM;    \
typedef Eigen::Map<Mat4XT,0,Eigen::OuterStride<> > Mat4XTM;  \
typedef Eigen::Map<const Mat4XT,0,Eigen::OuterStride<> > Mat4XTCM;  \
typedef Eigen::Map<MatX4T,0,Eigen::OuterStride<> > MatX4TM;  \
typedef Eigen::Map<const MatX4T,0,Eigen::OuterStride<> > MatX4TCM;  \
typedef Eigen::Map<Mat4T,0,Eigen::OuterStride<> > Mat4TM;  \
typedef Eigen::Map<const Mat4T,0,Eigen::OuterStride<> > Mat4TCM;    \
typedef Eigen::Map<Vec6T> Vec6TM;  \
typedef Eigen::Map<const Vec6T> Vec6TCM;  \
typedef Eigen::Map<Mat6X3T,0,Eigen::OuterStride<> > Mat6X3TM;  \
typedef Eigen::Map<const Mat6X3T,0,Eigen::OuterStride<> > Mat6X3TCM;  \
typedef Eigen::Map<Mat6XT,0,Eigen::OuterStride<> > Mat6XTM;  \
typedef Eigen::Map<const Mat6XT,0,Eigen::OuterStride<> > Mat6XTCM;  \
typedef Eigen::Map<MatX6T,0,Eigen::OuterStride<> > MatX6TM;  \
typedef Eigen::Map<const MatX6T,0,Eigen::OuterStride<> > MatX6TCM;  \
typedef Eigen::Map<Mat6T,0,Eigen::OuterStride<> > Mat6TM;  \
typedef Eigen::Map<const Mat6T,0,Eigen::OuterStride<> > Mat6TCM;

#define DECL_MAP_FUNCS  \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<T2,0,Eigen::OuterStride<> > mapM(T2& m) {   \
  return Eigen::Map<T2,0,Eigen::OuterStride<> >(m.data(),m.rows(),m.cols(),Eigen::OuterStride<>(m.rows()));  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<T2> mapV(T2& m) {   \
  return Eigen::Map<T2>(m.data(),m.rows());  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<T2,0,Eigen::OuterStride<> > mapM(T2* m) {   \
  static const int rows = T2::RowsAtCompileTime>0?T2::RowsAtCompileTime:0;  \
  static const int cols = T2::ColsAtCompileTime>0?T2::ColsAtCompileTime:0;  \
  return m?Eigen::Map<T2,0,Eigen::OuterStride<> >(m->data(),m->rows(),m->cols(),Eigen::OuterStride<>(m->rows())):   \
           Eigen::Map<T2,0,Eigen::OuterStride<> >(NULL,rows,cols,Eigen::OuterStride<>(rows));  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<T2> mapV(T2* m) {   \
  return m?Eigen::Map<T2>(m->data(),m->rows()):Eigen::Map<T2>(NULL,0);  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<const T2,0,Eigen::OuterStride<> > mapCM(const T2& m) {   \
  return Eigen::Map<const T2,0,Eigen::OuterStride<> >(m.data(),m.rows(),m.cols(),Eigen::OuterStride<>(m.rows()));  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<const T2> mapCV(const T2& m) {   \
  return Eigen::Map<const T2>(m.data(),m.rows());  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<const T2,0,Eigen::OuterStride<> > mapCM(const T2* m) {   \
  static const int rows = T2::RowsAtCompileTime>0?T2::RowsAtCompileTime:0;  \
  static const int cols = T2::ColsAtCompileTime>0?T2::ColsAtCompileTime:0;  \
  return m?Eigen::Map<const T2,0,Eigen::OuterStride<> >(m->data(),m->rows(),m->cols(),Eigen::OuterStride<>(m->rows())):   \
           Eigen::Map<const T2,0,Eigen::OuterStride<> >(NULL,rows,cols,Eigen::OuterStride<>(rows));  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<const T2> mapCV(const T2* m) {   \
  return m?Eigen::Map<const T2>(m->data(),m->rows()):Eigen::Map<const T2>(NULL,0);  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<const T2,0,Eigen::OuterStride<> > mapM2CM(Eigen::Map<T2,0,Eigen::OuterStride<> > m) {   \
  return Eigen::Map<const T2,0,Eigen::OuterStride<> >(m.data(),m.rows(),m.cols(),Eigen::OuterStride<>(m.outerStride()));  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<T2,0,Eigen::OuterStride<> > mapCM2M(Eigen::Map<const T2,0,Eigen::OuterStride<> > m) {   \
  return Eigen::Map<T2,0,Eigen::OuterStride<> >(m.data(),m.rows(),m.cols(),Eigen::OuterStride<>(m.outerStride()));  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<const T2> mapV2CV(Eigen::Map<T2> m) {   \
  return Eigen::Map<const T2>(m.data(),m.rows());  \
}   \
template <typename T2>   \
static FORCE_INLINE Eigen::Map<T2> mapCV2V(Eigen::Map<const T2> m) {   \
  return Eigen::Map<T2>(m.data(),m.rows());  \
}

#endif
