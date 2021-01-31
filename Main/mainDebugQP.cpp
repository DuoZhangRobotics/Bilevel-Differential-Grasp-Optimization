#include <Optimizer/QCQPSolverMosek.h>
#include <Optimizer/QCQPSolverGurobi.h>

USE_PRJ_NAMESPACE

typedef double T;
typedef Eigen::Matrix<T,-1,1> Vec;
typedef Eigen::Matrix<T,3,1> Vec3T;
typedef Eigen::Matrix<T,-1,-1> MatT;
typedef Eigen::SparseMatrix<T,0,sizeType> SMat;
template <typename MAT>
void debugQP(QCQPSolver<T>& sol,T rho,bool cone)
{
  int N=10,M=4;

  Vec g=Vec::Random(N);
  Vec x=Vec::Random(N),d,d2;
  MAT H=MatT::Random(N,N).sparseView();
  H=H*H.transpose();

  MAT cjac=MatT::Random(M,N).sparseView();
  Vec l=Vec::Constant(N,-0.1f);
  Vec u=Vec::Constant(N, 0.1f);
  Vec lA=Vec::Constant(M,-0.1f);
  Vec uA=Vec::Constant(M, 0.1f);
  l[0]=0;

  std::vector<Coli,Eigen::aligned_allocator<Coli>> QCones;
  if(!sol.solveQP(d=x,H,g,&cjac,&l,&u,&lA,&uA,QCones)) {
    WARNING("qp solver failed")
  } else {
    std::cout << "normSqr(d): " << d.squaredNorm() << std::endl;
    std::cout << "d: ";
    for(sizeType i=0; i<d.size(); i++)
      std::cout << d[i] << " ";
    std::cout << std::endl;

    std::cout << "d-l: ";
    for(sizeType i=0; i<d.size(); i++)
      std::cout << (d-l)[i] << " ";
    std::cout << std::endl;

    std::cout << "u-d: ";
    for(sizeType i=0; i<d.size(); i++)
      std::cout << (u-d)[i] << " ";
    std::cout << std::endl;

    std::cout << "(cjac*d-lA): ";
    for(sizeType i=0; i<cjac.rows(); i++)
      std::cout << (cjac*d-lA)[i] << " ";
    std::cout << std::endl;

    std::cout << "(uA-cjac*d): ";
    for(sizeType i=0; i<cjac.rows(); i++)
      std::cout << (uA-cjac*d)[i] << " ";
    std::cout << std::endl;
  }

  if(cone) {
    QCones.push_back(Vec3i(3,2,1));
    if(!sol.solveQP(d=x,H,g,NULL,NULL,NULL,NULL,NULL,QCones)) {
      WARNING("qp solver failed")
    } else {
      T val=std::pow((d+x)[QCones[0][0]],2);
      for(sizeType i=1; i<QCones[0].size(); i++)
        val-=std::pow((d+x)[QCones[0][i]],2);
      std::cout << "QConeErr: " << val << std::endl;
    }
    QCones.clear();
  }

  T TR=1;
  if(!sol.solveL1QP(d2=x,H,g,&cjac,&l,&u,&lA,&uA,TR,rho,QCones)) {
    WARNING("l1qp solver failed")
  } else {
    std::cout << "normSqr(d): " << d2.squaredNorm() << std::endl;
    d2=d2.segment(0,x.size()).eval();
    std::cout << "d-d2: ";
    for(sizeType i=0; i<d.size(); i++)
      std::cout << (d-d2)[i] << " ";
    std::cout << std::endl;
  }

  TR=0.01f;
  if(!sol.solveL1QP(d2=x,H,g,NULL,NULL,NULL,NULL,NULL,TR,rho,QCones)) {
    WARNING("l1qp solver failed")
  } else {
    std::cout << "normSqr(d): " << d2.squaredNorm() << std::endl;
  }
}
int main()
{
  {
    QCQPSolverMosek<T> sol;
    INFO("-------------------------------------------------------Mosek Dense")
    srand(0);
    debugQP<MatT>(sol,10,false);
    INFO("-------------------------------------------------------Mosek Sparse")
    srand(0);
    debugQP<SMat>(sol,10,false);
  }
  {
    QCQPSolverGurobi<T> sol;
    INFO("-------------------------------------------------------Gurobi Dense")
    srand(0);
    debugQP<MatT>(sol,10,true);
    INFO("-------------------------------------------------------Gurobi Sparse")
    srand(0);
    debugQP<SMat>(sol,10,true);
  }
  return 0;
}
