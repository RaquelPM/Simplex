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
#include "GS.h"

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
  srand(seed);
  int n = std::atoi(argv[1]);

  // Eigen::MatrixXd B_dense = gen_random_non_singular_mat(n, 20);

  // commprimindo a matriz, convertendo-a para uma matriz esparsa (seboso, nao faça)
  // Eigen::SparseMatrix<double> B = B_dense.sparseView();

  // exemplo do livro Chvatel
  Eigen::SparseMatrix<double> B(n, n);
  B.setIdentity();

  VectorXd Et(n);
  Et << 1, 1, 3;
  pair<int, VectorXd> p_d = make_pair(1, Et);

  GS *g = new GS(B, n);
  g->addEk(p_d);

  VectorXd a(n);
  a << 3, 1, 4;
  VectorXd d = g->FTRAN(a);
  cout << d << endl;
  VectorXd c(n);
  c << 0, 12, 0;
  VectorXd y = g->BTRAN(c);
  cout << y << endl;

  g->refatorar();

  d = VectorXd::Zero(n);
  y = VectorXd::Zero(n);

  d = g->FTRAN(a);
  cout << d << endl;

  y = g->BTRAN(c);
  cout << y << endl;


  delete g;


  // // gerando matriz E aleatoria (apenas a coluna e sua localizacao sao necessarias para representa-la)
  // auto [eta_idx, eta_col] = gen_random_eta_mat(n);

  // // matriz E representada na tora, apenas para fins ilustrativos
  // Eigen::MatrixXd E = Eigen::MatrixXd::Identity(n, n);
  // E.col(eta_idx) = eta_col;

  return 0;
}