#include <Utils/ArticulatedBodyPragma.h>
#include <Utils/SpatialRotationUtil.h>

USE_PRJ_NAMESPACE

typedef mpfr::mpreal T;
int main()
{
  mpfr_set_default_prec(1024U);
  RandEngine::useDeterministic();
  RandEngine::seed(0);
  debugCross<T>();
  debugSpatial<T>();
  debugRotation<T>();
  debugSpatialRotation<T>();
  return 0;
}

