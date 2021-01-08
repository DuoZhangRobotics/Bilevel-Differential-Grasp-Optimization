#include <Deformable/FEMCell.h>

USE_PRJ_NAMESPACE

template <typename T>
void debug()
{
  FEMCell<T>::debugDFDX(3);
  FEMCell<T>::debugDFDX(4);
  FEMCell<T>::debugDFDX(8);
  FEMCell<T>::debugDFDX(8,0);
}
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debug<double>();
  return 0;
}
