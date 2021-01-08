#include <Deformable/FEMEnergy.h>

USE_PRJ_NAMESPACE

template <typename T>
void debug()
{
  FEMConstantForceEnergy<T>::debugDEDX(3);
  FEMConstantForceEnergy<T>::debugDEDX(4);
  FEMConstantForceEnergy<T>::debugDEDX(8);

  FEMPressureEnergy<T>::debugDEDX();

  FEMEmbeddedPressureEnergy<T>::debugDEDX(3);
  FEMEmbeddedPressureEnergy<T>::debugDEDX(4);
  FEMEmbeddedPressureEnergy<T>::debugDEDX(8);

  FEMLinearElasticEnergy<T>::debugDEDX(3,0.1f);
  FEMLinearElasticEnergy<T>::debugDEDX(4,0.1f);
  FEMLinearElasticEnergy<T>::debugDEDX(8,0.1f);
  FEMLinearElasticEnergy<T>::debugDEDX(8,0.1f,0);

  FEMStVKElasticEnergy<T>::debugDEDX(3,0.1f);
  FEMStVKElasticEnergy<T>::debugDEDX(4,0.1f);
  FEMStVKElasticEnergy<T>::debugDEDX(8,0.1f);
  FEMStVKElasticEnergy<T>::debugDEDX(8,0.1f,0);

  FEMCorotatedElasticEnergy<T>::debugDEDX(3,0.1f);
  FEMCorotatedElasticEnergy<T>::debugDEDX(4,0.1f);
  FEMCorotatedElasticEnergy<T>::debugDEDX(8,0.1f);
  FEMCorotatedElasticEnergy<T>::debugDEDX(8,0.1f,0);

  //FEMNonHookeanElasticEnergy<T>::debugDEDX(3,0.1f);
  FEMNonHookeanElasticEnergy<T>::debugDEDX(4,0.1f);
  FEMNonHookeanElasticEnergy<T>::debugDEDX(8,0.1f);
  FEMNonHookeanElasticEnergy<T>::debugDEDX(8,0.1f,0);
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debug<double>();
  return 0;
}
