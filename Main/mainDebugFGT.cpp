#include <Quasistatic/FGTTreeNode.h>

USE_PRJ_NAMESPACE

typedef double T;
int main()
{
  //FGTTreeNode<T>::debugSwapId(10,1000);
  FGTTreeNode<T>::debugTree(10,2,true,true);
  FGTTreeNode<T>::debugTree(10,2,true,false);
  //FGTTreeNode<T>::debugTaylor(10,2,true,0.1f,1e-6f,true);
  //FGTTreeNode<T>::debugTaylor(10,2,true,0.1f,1e-6f,false);
  FGTTreeNode<T>::debugFGT(10,2,true,0.01f,1e-3f,true);
  FGTTreeNode<T>::debugFGT(10,2,true,0.01f,1e-3f,false);
  FGTTreeNode<T>::debugFGT(10,2,true,0.1f,1e-3f,true);
  FGTTreeNode<T>::debugFGT(10,2,true,0.1f,1e-3f,false);
  FGTTreeNode<T>::debugFGT(10,2,true,1.0f,1e-3f,true);
  FGTTreeNode<T>::debugFGT(10,2,true,1.0f,1e-3f,false);
  //std::cout << FGTTreeNode<T>::combination(10,3) << std::endl;
  return 0;
}
