#include <Articulated/MultiPrecisionSeparatingPlane.h>

USE_PRJ_NAMESPACE

typedef double T;
int main()
{
  typedef MultiPrecisionSeparatingPlane<T>::Vec4T Vec4T;
  typedef MultiPrecisionSeparatingPlane<T>::Vec3T Vec3T;
  typedef MultiPrecisionSeparatingPlane<T>::PSS PSS;
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  Options ops;
  PSS pss;
  for(sizeType i=0; i<100; i++)
    pss.push_back(Vec3T::Random());
  Vec4T plane=Vec4T::Random();
  bool succ;
  MultiPrecisionSeparatingPlane<T> sol(ops,plane,5.0f);
  sol.resetPoints(pss);
  sol.debugGradient();
  sol.solve(succ);
  return 0;
}

