#ifndef SVD_H
#define SVD_H

#include <CommonFile/MathBasic.h>
#include "DebugGradient.h"
#include <Eigen/Eigen>

PRJ_BEGIN

//derivative of SVD decomposition
template <typename T>
void adjustF3D(typename ScalarUtil<T>::ScalarMat3& F,
               typename ScalarUtil<T>::ScalarMat3& U,
               typename ScalarUtil<T>::ScalarMat3& V,
               typename ScalarUtil<T>::ScalarVec3& S)
{
  Mat3d FTF=(F.transpose()*F).unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  });
  Eigen::SelfAdjointEigenSolver<Mat3d> eig(FTF);
  //F Hat
  S=eig.eigenvalues().cwiseAbs().cwiseSqrt().cast<T>();
  //adjust V
  V=eig.eigenvectors().cast<T>();
  if(V.determinant() < 0.0f)
    V.col(0)*=-1.0f;
  //adjust U
  sizeType good[3];
  sizeType nrGood=0;
  for(sizeType v=0; v<3; v++)
    if(S[v] > ScalarUtil<T>::scalar_eps()) {
      U.col(v)=F*V.col(v)/S[v];
      good[nrGood++]=v;
    }
  //the tet is a point
  if(nrGood == 0)
    U.setIdentity();
  //set columns of U to be orthogonal
  else for(sizeType v=0; v<3; v++)
      if(S[v] <= ScalarUtil<T>::scalar_eps()) {
        if(nrGood == 1) {
          typename ScalarUtil<T>::ScalarVec3 d;
          T nd;
          while(true) {
            d=U.col(good[0]).cross(ScalarUtil<T>::ScalarVec3::Random());
            nd=std::sqrt(d.squaredNorm());
            if(nd > 1E-3f) {
              U.col(v)=d/nd;
              break;
            }
          }
        } else {
          ASSERT(nrGood == 2)
          U.col(v)=U.col(good[0]).cross(U.col(good[1]));
        }
        good[nrGood++]=v;
      }
  //negate negative U columns
  if(U.determinant() < 0.0f) {
    sizeType id;
    S.minCoeff(&id);
    S[id]*=-1.0f;
    U.col(id)*=-1.0f;
  }
  //rebuild F
  typename ScalarUtil<T>::ScalarMat3 Sigma=ScalarUtil<T>::ScalarMat3::Zero();
  Sigma.diagonal()=S;
  F=U*Sigma*V.transpose();
}
template <typename T>
void adjustF2D(typename ScalarUtil<T>::ScalarMat3& F,
               typename ScalarUtil<T>::ScalarMat3& U,
               typename ScalarUtil<T>::ScalarMat3& V,
               typename ScalarUtil<T>::ScalarVec3& S)
{
  Mat2d FTF=(F.transpose()*F).unaryExpr([&](const T& in) {
    return (scalarD)std::to_double(in);
  }).template block<2,2>(0,0);
  Eigen::SelfAdjointEigenSolver<Mat2d> eig(FTF);
  //F Hat
  S.setZero();
  S.template block<2,1>(0,0)=eig.eigenvalues().cwiseAbs().cwiseSqrt().cast<T>();
  //adjust V
  V.setZero();
  V.template block<2,2>(0,0)=eig.eigenvectors().cast<T>();
  if(V.template block<2,2>(0,0).determinant() < 0.0f)
    V.col(0)*=-1.0f;
  //adjust U
  U.setZero();
  sizeType good[2];
  sizeType nrGood=0;
  for(sizeType v=0; v<2; v++)
    if(S[v] > ScalarUtil<T>::scalar_eps()) {
      U.col(v)=F*V.col(v)/S[v];
      good[nrGood++]=v;
    }
  //the tet is a point
  if(nrGood == 0)
    U.template block<2,2>(0,0).setIdentity();
  //set columns of U to be orthogonal
  else for(sizeType v=0; v<2; v++)
      if(S[v] <= ScalarUtil<T>::scalar_eps()) {
        U.col(v)=typename ScalarUtil<T>::ScalarVec3(-U(1,good[0]),U(0,good[0]),0.0f);
        good[nrGood++]=v;
      }
  //negate negative U columns
  if(U.template block<2,2>(0,0).determinant() < 0.0f) {
    sizeType id;
    S.template block<2,1>(0,0).minCoeff(&id);
    S[id]*=-1.0f;
    U.col(id)*=-1.0f;
  }
  //rebuild F
  typename ScalarUtil<T>::ScalarMat3 Sigma=ScalarUtil<T>::ScalarMat3::Zero();
  Sigma.diagonal()=S;
  F=U*Sigma*V.transpose();
  F(2,2)=1;
}
template <typename T>
void derivSVD(typename ScalarUtil<T>::ScalarMat3& derivU,
              typename ScalarUtil<T>::ScalarVec3& derivD,
              typename ScalarUtil<T>::ScalarMat3& derivV,
              const typename ScalarUtil<T>::ScalarMat3& LHS,
              const typename ScalarUtil<T>::ScalarMat3& U,
              const typename ScalarUtil<T>::ScalarVec3& D,
              const typename ScalarUtil<T>::ScalarMat3& V)
{
#define SOLVE_SAFE(a,b,c) \
tmp << D(c),D(b),D(b),D(c); \
if(std::abs(tmp(0,0)-tmp(0,1)) < 1E-6f)	\
  tmpInv=(tmp*tmp+reg).inverse()*tmp; \
else  \
  tmpInv=tmp.inverse(); \
a=tmpInv*typename ScalarUtil<T>::ScalarVec2(LHS(b,c),-LHS(c,b));

  derivD=typename ScalarUtil<T>::ScalarVec3(LHS(0,0),LHS(1,1),LHS(2,2));

  //21,31,32
  typename ScalarUtil<T>::ScalarMat2 reg=ScalarUtil<T>::ScalarMat2::Identity()*1E-6f,tmp,tmpInv;
  typename ScalarUtil<T>::ScalarVec2 omg10,omg20,omg21;
  SOLVE_SAFE(omg10,1,0)
  SOLVE_SAFE(omg20,2,0)
  SOLVE_SAFE(omg21,2,1)
#undef SOLVE_SAFE

  typename ScalarUtil<T>::ScalarMat3 omgU;
  omgU<<      0,-omg10(0),-omg20(0),
       omg10(0),        0,-omg21(0),
       omg20(0), omg21(0),        0;
  derivU=U*omgU;

  typename ScalarUtil<T>::ScalarMat3 omgV;
  omgV<<      0,-omg10(1),-omg20(1),
       omg10(1),        0,-omg21(1),
       omg20(1), omg21(1),        0;
  derivV=-V*omgV;
}
template <typename T>
void derivSVDDir(typename ScalarUtil<T>::ScalarMat3& derivU,
                 typename ScalarUtil<T>::ScalarVec3& derivD,
                 typename ScalarUtil<T>::ScalarMat3& derivV,
                 const typename ScalarUtil<T>::ScalarMat3& dir,
                 const typename ScalarUtil<T>::ScalarMat3& U,
                 const typename ScalarUtil<T>::ScalarVec3& D,
                 const typename ScalarUtil<T>::ScalarMat3& V)
{
  typename ScalarUtil<T>::ScalarMat3 LHS=ScalarUtil<T>::ScalarMat3::Zero();
  for(sizeType r=0; r<3; r++)
    for(sizeType c=0; c<3; c++)
      LHS+=U.row(r).transpose()*V.row(c)*dir(r,c);
  derivSVD<T>(derivU,derivD,derivV,LHS,U,D,V);
}
template <typename T>
void derivRDF(typename ScalarUtil<T>::ScalarMat9& DRDF,
              const typename ScalarUtil<T>::ScalarMat3& F,
              const typename ScalarUtil<T>::ScalarMat3& U,
              const typename ScalarUtil<T>::ScalarVec3& D,
              const typename ScalarUtil<T>::ScalarMat3& V)
{
  typename ScalarUtil<T>::ScalarVec3 derivD;
  typename ScalarUtil<T>::ScalarMat3 derivU,derivV;
  if(D.minCoeff() < ScalarUtil<T>::scalar_eps())
    DRDF.setZero();
  for(sizeType r=0; r<3; r++)
    for(sizeType c=0; c<3; c++) {
      derivSVD<T>(derivU,derivD,derivV,U.row(r).transpose()*V.row(c),U,D,V);
      Eigen::Map<typename ScalarUtil<T>::ScalarMat3>(DRDF.col(r+c*3).data())=(derivU*V.transpose()+U*derivV.transpose()).template cast<T>();
    }
}
template <typename T>
void debugDerivSVD(bool debugSVD,bool debugRDF,int maxIt)
{
  DEFINE_NUMERIC_DELTA
  for(int it=0; it<maxIt; it++) {
    //debug derivative of SVD decomposition
    if(debugSVD) {
      typename ScalarUtil<T>::ScalarVec3 dD;
      typename ScalarUtil<T>::ScalarMat3 F=ScalarUtil<T>::ScalarMat3::Random(),dF=ScalarUtil<T>::ScalarMat3::Random(),dU,dV;
      Eigen::JacobiSVD<typename ScalarUtil<T>::ScalarMat3> svd(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
      derivSVDDir<T>(dU,dD,dV,dF,svd.matrixU(),svd.singularValues(),svd.matrixV());
      Eigen::JacobiSVD<typename ScalarUtil<T>::ScalarMat3> svd2(F+dF*DELTA,Eigen::ComputeFullU|Eigen::ComputeFullV);
      INFOV("DerivSVD: %f %f %f Err: %f %f %f",dU.norm(),dD.norm(),dV.norm(),
            (dU-(svd2.matrixU()-svd.matrixU())/DELTA).norm(),
            (dD-(svd2.singularValues()-svd.singularValues())/DELTA).norm(),
            (dV-(svd2.matrixV()-svd.matrixV())/DELTA).norm())
    }

    //debug derivative of RS (polar) decomposition
    if(debugRDF) {
      typename ScalarUtil<T>::ScalarMat9 DRDF;
      typename ScalarUtil<T>::ScalarMat3 F=ScalarUtil<T>::ScalarMat3::Random(),dF=ScalarUtil<T>::ScalarMat3::Random();
      Eigen::JacobiSVD<typename ScalarUtil<T>::ScalarMat3> svd(F,Eigen::ComputeFullU|Eigen::ComputeFullV);
      derivRDF<T>(DRDF,F,svd.matrixU(),svd.singularValues(),svd.matrixV());
      Eigen::JacobiSVD<typename ScalarUtil<T>::ScalarMat3> svd2(F+dF*DELTA,Eigen::ComputeFullU|Eigen::ComputeFullV);
      typename ScalarUtil<T>::ScalarMat3 R=svd.matrixU()*svd.matrixV().transpose();
      typename ScalarUtil<T>::ScalarMat3 R2=svd2.matrixU()*svd2.matrixV().transpose(),DRDFN=(R2-R)/DELTA;
      typename ScalarUtil<T>::ScalarCol DRDFA=DRDF.template cast<T>()*Eigen::Map<typename ScalarUtil<T>::ScalarCol>(dF.data(),9);
      INFOV("DerivRS: %f Err: %f",DRDFA.norm(),(DRDFA-Eigen::Map<typename ScalarUtil<T>::ScalarCol>(DRDFN.data(),9)).norm())
    }
  }
}

PRJ_END

#endif
