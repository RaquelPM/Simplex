#include <iostream>
#include <cstdlib>
#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

#include <numeric>
#include <ostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <map>

#include <stdio.h>
#include <umfpack.h>

using Eigen::MatrixXd;
using namespace std;

#include "Simplex.h"

// gerador de matrizes inversiveis B aleatorias (it_max e o grau de aleatoriedade)
Eigen::MatrixXd gen_random_non_singular_mat(int n, int it_max)
{

  Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);

  std::vector<int> idx_vec(n);
  std::iota(idx_vec.begin(), idx_vec.end(), 1);

  for (int k = 0; k < it_max; k++)
  {
    int r1 = idx_vec[rand() % n] - 1;
    std::swap(idx_vec.back(), idx_vec[r1]);
    int r2 = idx_vec[rand() % (n - 1)] - 1;

    I.row(r1) += ((double)rand() / RAND_MAX) * (rand() % 2 ? -1 : 1) * I.row(r2);
  }

  return I;
}

// cria matriz E aleatoria
std::pair<int, Eigen::VectorXd> gen_random_eta_mat(int n)
{
  int p = rand() % n;
  Eigen::VectorXd d = Eigen::VectorXd::Random(n);

  while (std::abs(d[p]) < 0.000001)
  {
    d = Eigen::VectorXd::Random(n);
  }

  return std::make_pair(p, d);
}

int main(int argc, char **argv)
{
  int seed = std::atoi(argv[2]);
  // srand(time(NULL));
  srand(seed);

  int n = std::atoi(argv[1]);

  // Eigen::MatrixXd B_dense = gen_random_non_singular_mat(n, 20);

  // commprimindo a matriz, convertendo-a para uma matriz esparsa (seboso, nao faça)
  // Eigen::SparseMatrix<double> B = B_dense.sparseView();

  Eigen::SparseMatrix<double> B(n, n);
  B.setIdentity();
  vector<pair<int, Eigen::VectorXd>> Ek;

  VectorXd Et(n);
  Et << 1, 1, 3;
  pair<int, VectorXd> p_d = make_pair(1, Et);
  Ek.push_back(p_d);

  Simplex s = Simplex(B, n, Ek);

  // Criando decomposicao LU para a matriz esparsa B usando UMFPACK
  double *null = (double *)NULL;
  void *Symbolic, *Numeric;

  (void)umfpack_di_symbolic(n, n, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, null, null);
  (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, null, null);

  // Resolvendo o sistema Bd = a varias vezes, reutilizando a mesma decomposicao LU de B para varios vetores a diferentes
  Eigen::VectorXd d(n);

  for (int i = 0; i < 100; i++)
  {
    Eigen::VectorXd a = Eigen::VectorXd::Random(n);

    (void)umfpack_di_solve(UMFPACK_A, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), d.data(), a.data(), Numeric, null, null);

    // verificando se de fato Bd = a
    if ((B * d - a).norm() > 0.0000001)
    {
      std::cout << "Deu errado" << std::endl;
      exit(0);
    }
  }

  umfpack_di_free_symbolic(&Symbolic);
  umfpack_di_free_numeric(&Numeric);

  // gerando matriz E aleatoria (apenas a coluna e sua localizacao sao necessarias para representa-la)
  auto [eta_idx, eta_col] = gen_random_eta_mat(n);

  // matriz E representada na tora, apenas para fins ilustrativos
  Eigen::MatrixXd E = Eigen::MatrixXd::Identity(n, n);
  E.col(eta_idx) = eta_col;

  cout << E << endl;
  cout << eta_idx << endl;

  std::cout << "Antiga matriz B:\n"
            << B.toDense() << std::endl
            << std::endl;
  std::cout << "Nova matriz B:\n"
            << B * E << std::endl;

  return 0;
}

// n = 5 ;
// int Ap [ ] = {0, 2, 5, 9, 10, 12} ;
// int Ai [ ] = { 0, 1, 0, 2, 4, 1, 2, 3, 4, 2, 1, 4} ;
// double Ax [ ] = {2., 3., 3., -1., 4., 4., -3., 1., 2., 2., 6., 1.} ;
// double b[ ] = {8., 45., -3., 3., 19.} ;
// double x[n] ;

// double *null = (double *) NULL ;
// void *Symbolic, *Numeric ;
// (void) umfpack_di_symbolic (n, n, Ap, Ai, Ax, &Symbolic, null, null) ;
// (void) umfpack_di_numeric (Ap, Ai, Ax, Symbolic, &Numeric, null, null) ;
// umfpack_di_free_symbolic (&Symbolic) ;
// (void) umfpack_di_solve (UMFPACK_A, Ap, Ai, Ax, x, b, Numeric, null, null) ;
// umfpack_di_free_numeric (&Numeric) ;
// for (i = 0 ; i < n ; i++) printf ("x [%d] = %g\n", i, x [i]) ;
// return (0) ;
