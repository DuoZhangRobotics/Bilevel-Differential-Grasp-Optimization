#ifndef INTERP_H
#define INTERP_H

#include "MathBasic.h"

PRJ_BEGIN

//linear interpolation
template <typename T,typename TF>
T interp1D(const T& v0,const T& v1,const TF& px)
{
  return v0*(typename EigenTraits<T>::ScalarType)(1.0f-px)+
         v1*(typename EigenTraits<T>::ScalarType)px;
}
template <typename T,typename TF>
T interp2D(const T& v0,const T& v1,const T& v2,const T& v3,const TF& px,const TF& py)
{
  return interp1D(interp1D(v0,v1,px),
                  interp1D(v2,v3,px),py);
}
template <typename T,typename TF>
T interp3D(const T& v0,const T& v1,const T& v2,const T& v3,
           const T& v4,const T& v5,const T& v6,const T& v7,
           const TF& px,const TF& py,const TF& pz)
{
  return interp1D(interp2D(v0,v1,v2,v3,px,py),
                  interp2D(v4,v5,v6,v7,px,py),pz);
}
//linear interpolation stencil
template <typename TF>
sizeType stencil1D(TF* coefs,const TF& px)
{
  coefs[0]=1.0f-px;
  coefs[1]=px;
  return 2;
}
template <typename TF>
sizeType stencil2D(TF* coefs,const TF& px,const TF& py)
{
  coefs[0]=(1.0f-px)*(1.0f-py);
  coefs[1]=px*(1.0f-py);
  coefs[2]=(1.0f-px)*py;
  coefs[3]=px*py;
  return 4;
}
template <typename TF>
sizeType stencil3D(TF* coefs,const TF& px,const TF& py,const TF& pz)
{
  coefs[0]=(1.0f-px)*(1.0f-py)*(1.0f-pz);
  coefs[1]=px*(1.0f-py)*(1.0f-pz);
  coefs[2]=(1.0f-px)*py*(1.0f-pz);
  coefs[3]=px*py*(1.0f-pz);
  coefs[4]=(1.0f-px)*(1.0f-py)*pz;
  coefs[5]=px*(1.0f-py)*pz;
  coefs[6]=(1.0f-px)*py*pz;
  coefs[7]=px*py*pz;
  return 8;
}
//linear interpolation with minmax
template <typename T,typename TF>
T interp1D(const T& v0,const T& v1,const TF& px,T& minV,T& maxV)
{
  minV=compMin(v0,v1);
  maxV=compMax(v0,v1);
  return v0*(typename EigenTraits<T>::ScalarType)(1.0f-px)+
         v1*(typename EigenTraits<T>::ScalarType)px;
}
template <typename T,typename TF>
T interp2D(const T& v0,const T& v1,const T& v2,const T& v3,const TF& px,const TF& py,T& minV,T& maxV)
{
  T minV2,maxV2;
  T ret=interp1D(interp1D(v0,v1,px,minV,maxV),
                 interp1D(v2,v3,px,minV2,maxV2),py);
  minV=compMin(minV,minV2);
  maxV=compMax(maxV,maxV2);
  return ret;
}
template <typename T,typename TF>
T interp3D(const T& v0,const T& v1,const T& v2,const T& v3,
           const T& v4,const T& v5,const T& v6,const T& v7,
           const TF& px,const TF& py,const TF& pz,T& minV,T& maxV)
{
  T minV2,maxV2;
  T ret=interp1D(interp2D(v0,v1,v2,v3,px,py,minV,maxV),
                 interp2D(v4,v5,v6,v7,px,py,minV2,maxV2),pz);
  minV=compMin(minV,minV2);
  maxV=compMax(maxV,maxV2);
  return ret;
}
//cubic interpolation
template <typename T>
Eigen::Matrix<T,4,4> cubicHPosition()
{
  Eigen::Matrix<T,4,4> ret;
  ret <<
      0,-1, 2,-1,
      2, 0,-5, 3,
      0, 1, 4,-3,
      0, 0,-1, 1;
  return ret/2;
}
template <typename T>
Eigen::Matrix<T,4,4> cubicHTangent()
{
  Eigen::Matrix<T,4,4> ret;
  ret <<
      1,0,-3, 2,
      0,1,-2, 1,
      0,0, 3,-2,
      0,0,-1, 1;
  return ret;
}
template <typename T>
Eigen::Matrix<T,4,1> cubicHPosition(T px)
{
  return cubicHPosition<T>()*Eigen::Matrix<T,4,1>(1,px,px*px,px*px*px);
}
template <typename T>
Eigen::Matrix<T,4,1> cubicHTangent(T px)
{
  return cubicHTangent<T>()*Eigen::Matrix<T,4,1>(1,px,px*px,px*px*px);
}
template <typename T,typename TF>
T interp1DCubic(const T a[4],TF px,T* valLinear,T* minV,T* maxV,bool linear=false)
{
  if(linear)
    return interp1D(a[1],a[2],px);
  Eigen::Matrix<TF,4,1> coef=cubicHPosition<TF>(px);
  T ret=EigenTraits<T>::value();
  for(sizeType r=0; r<4; r++)
    ret+=a[r]*(typename EigenTraits<T>::ScalarType)coef[r];
  if(valLinear) {
    *valLinear=interp1D<T,TF>(a[1],a[2],px);
    *minV=compMin(compMin(*minV,a[1]),a[2]);
    *maxV=compMax(compMax(*maxV,a[1]),a[2]);
  }
  return ret;
}
template <typename T,typename TF,typename GRID>
T interp1DCubic(const sizeType index[4],const GRID& grid,const TF& px,T* valLinear,T* minV,T* maxV,bool linear=false)
{
  const T a[4]= {grid[index[0]],grid[index[1]],grid[index[2]],grid[index[3]]};
  return interp1DCubic<T,TF>(a,px,valLinear,minV,maxV,linear);
}
template <typename T,typename TF,typename GRID>
T interp2DCubic(const sizeType index[4][4],const GRID& grid,const TF& px,const TF& py,T* valLinear,T* minV,T* maxV,bool linear=false)
{
  T aL[2];
  const T a[4]= {
    interp1DCubic<T,TF,GRID>(index[0],grid,px,NULL,NULL,NULL),
    interp1DCubic<T,TF,GRID>(index[1],grid,px,valLinear ? aL+0 : NULL,minV,maxV,linear),
    interp1DCubic<T,TF,GRID>(index[2],grid,px,valLinear ? aL+1 : NULL,minV,maxV,linear),
    interp1DCubic<T,TF,GRID>(index[3],grid,px,NULL,NULL,NULL)
  };
  if(linear)
    return interp1D(a[1],a[2],py);
  if(valLinear)
    *valLinear=interp1D(aL[0],aL[1],py);
  return interp1DCubic<T,TF>(a,py,NULL,NULL,NULL);
}
template <typename T,typename TF,typename GRID>
T interp3DCubic(const sizeType index[4][4][4],const GRID& grid,const TF& px,const TF& py,const TF& pz,T* valLinear,T* minV,T* maxV)
{
  T aL[2];
  const T a[4]= {
    interp2DCubic<T,TF,GRID>(index[0],grid,px,py,NULL,NULL,NULL),
    interp2DCubic<T,TF,GRID>(index[1],grid,px,py,valLinear ? aL+0 : NULL,minV,maxV),
    interp2DCubic<T,TF,GRID>(index[2],grid,px,py,valLinear ? aL+1 : NULL,minV,maxV),
    interp2DCubic<T,TF,GRID>(index[3],grid,px,py,NULL,NULL,NULL)
  };
  if(valLinear)
    *valLinear=interp1D(aL[0],aL[1],pz);
  return interp1DCubic<T,TF>(a,pz,NULL,NULL,NULL);
}
//cubic interpolation stencil
template <typename TF>
void stencil1DCubic(TF coefs[4],const TF& px,bool linear)
{
  if(linear) {
    coefs[0]=0;
    coefs[1]=1.0f-px;
    coefs[2]=px;
    coefs[3]=0;
  } else {
    Eigen::Matrix<TF,4,1> c=cubicHPosition<TF>(px);
    for(sizeType i=0;i<4;i++)
      coefs[i]=c[i];
  }
}
template <typename TF>
void stencil2DCubic(TF coefs[4][4],const TF& px,const TF& py,bool linear)
{
  stencil1DCubic(coefs[0],px,linear);
  stencil1DCubic(coefs[1],px,linear);
  stencil1DCubic(coefs[2],px,linear);
  stencil1DCubic(coefs[3],px,linear);
  TF coefTmp[4];
  stencil1DCubic(coefTmp,py,linear);
  for(sizeType x=0; x<4; x++)
    for(sizeType y=0; y<4; y++)
      coefs[x][y]*=coefTmp[x];
}
template <typename TF>
void stencil3DCubic(TF coefs[4][4][4],const TF& px,const TF& py,const TF& pz,bool linear)
{
  stencil2DCubic(coefs[0],px,py,linear);
  stencil2DCubic(coefs[1],px,py,linear);
  stencil2DCubic(coefs[2],px,py,linear);
  stencil2DCubic(coefs[3],px,py,linear);
  TF coefTmp[4];
  stencil1DCubic(coefTmp,pz,linear);
  for(sizeType x=0; x<4; x++)
    for(sizeType y=0; y<4; y++)
      for(sizeType z=0; z<4; z++)
        coefs[x][y][z]*=coefTmp[x];
}
//cubic interpolation gradient
template <typename T,typename TF>
T interp1DGradCubic(const T a[4],TF px,bool linear)
{
  if(linear)
    return a[2]-a[1];
  else {
    Eigen::Matrix<TF,4,1> coef=cubicHPosition<TF>(px);
    T ret=EigenTraits<T>::value();
    for(sizeType r=0; r<4; r++)
      ret+=a[r]*(typename EigenTraits<T>::ScalarType)coef[r];
    return ret;
  }
}
template <typename T,typename TF,typename GRID>
Eigen::Matrix<T,1,1> interp1DGradCubic(const sizeType index[4],const GRID& grid,const TF& px,bool linear)
{
  const T a[4]= {
    grid[index[0]],
    grid[index[1]],
    grid[index[2]],
    grid[index[3]]
  };
  Eigen::Matrix<T,1,1> ret;
  ret[0]=interp1DGradCubic(a,px,linear);
  return ret;
}
template <typename T,typename TF,typename GRID>
Eigen::Matrix<T,2,1> interp2DGradCubic(const sizeType index[4][4],const GRID& grid,const TF& px,const TF& py,bool linear)
{
  const T aGrad[4]= {
    interp1DGradCubic<T,TF,GRID>(index[0],grid,px,linear)[0],
    interp1DGradCubic<T,TF,GRID>(index[1],grid,px,linear)[0],
    interp1DGradCubic<T,TF,GRID>(index[2],grid,px,linear)[0],
    interp1DGradCubic<T,TF,GRID>(index[3],grid,px,linear)[0],
  };
  const T a[4]= {
    interp1DCubic<T,TF,GRID>(index[0],grid,px,NULL,NULL,NULL,linear),
    interp1DCubic<T,TF,GRID>(index[1],grid,px,NULL,NULL,NULL,linear),
    interp1DCubic<T,TF,GRID>(index[2],grid,px,NULL,NULL,NULL,linear),
    interp1DCubic<T,TF,GRID>(index[3],grid,px,NULL,NULL,NULL,linear),
  };
  Eigen::Matrix<T,2,1> ret;
  ret[0]=interp1DCubic<T,TF>(aGrad,py,NULL,NULL,NULL,linear);
  ret[1]=interp1DGradCubic<T,TF>(a,py,linear);
  return ret;
}
template <typename T,typename TF,typename GRID>
Eigen::Matrix<T,3,1> interp3DGradCubic(const sizeType index[4][4][4],const GRID& grid,const TF& px,const TF& py,const TF& pz,bool linear)
{
  Eigen::Matrix<T,4,2> aGrad;
  aGrad.row(0)=interp2DGradCubic<T,TF,GRID>(index[0],grid,px,py,linear);
  aGrad.row(1)=interp2DGradCubic<T,TF,GRID>(index[1],grid,px,py,linear);
  aGrad.row(2)=interp2DGradCubic<T,TF,GRID>(index[2],grid,px,py,linear);
  aGrad.row(3)=interp2DGradCubic<T,TF,GRID>(index[3],grid,px,py,linear);
  const T a[4]= {
    interp2DCubic<T,TF,GRID>(index[0],grid,px,py,NULL,NULL,NULL,linear),
    interp2DCubic<T,TF,GRID>(index[1],grid,px,py,NULL,NULL,NULL,linear),
    interp2DCubic<T,TF,GRID>(index[2],grid,px,py,NULL,NULL,NULL,linear),
    interp2DCubic<T,TF,GRID>(index[3],grid,px,py,NULL,NULL,NULL,linear),
  };
  Eigen::Matrix<T,3,1> ret;
  ret[0]=interp1DCubic<T,TF>(aGrad.data(),pz,NULL,NULL,NULL,linear);
  ret[1]=interp1DCubic<T,TF>(aGrad.data()+4,pz,NULL,NULL,NULL,linear);
  ret[2]=interp1DGradCubic<T,TF>(a,pz,linear);
  return ret;
}

PRJ_END

#endif
