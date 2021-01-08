#include <Deformable/NNHTP.h>

USE_PRJ_NAMESPACE

int main()
{
  NNHTP<scalarD>::MatT A;
  NNHTP<scalarD>::Vec b;
  A.setRandom(1000,100);
  b.setRandom(1000);

  NNHTP<scalarD> sol(A,b);
  sol.solve(10,1e-3f,true);
  return 0;
}

