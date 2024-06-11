#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <map>
#include <limits>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

#include "Simplex.h"
#include "GS.h"
#include "Data.h"

using namespace std;

double pInf = numeric_limits<double>::infinity();
double nInf = -numeric_limits<double>::infinity();

int main(int argc, char **argv)
{
  // int seed = std::atoi(argv[2]);
  // srand(seed);
  // int n = std::atoi(argv[1]);

  MatrixXd A(2, 12);
  A.row(0) << 3, 1, 5, 6, 9, 4, 3, 4, 7, 6, 4, 5;
  A.row(1) << 1, 0, 9, 5, 8, 1, 2, 7, 8, 7, 9, 1;

  VectorXd b(2);
  b << 72, 62;

  VectorXd c(12);
  c << 2, 1, -2, -2, 3, 2, 3, -4, 0, -2, -3, 3;

  VectorXd u(12);
  u << pInf, 3, -2, 3, 5, 1, pInf, pInf, 0, 5, pInf, pInf;

  VectorXd l(12);
  l << -5, nInf, -4, -2, 2, 0, 0, 3, nInf, nInf, nInf, nInf;

  vector<int> B = {0, 1};
  vector<int> N = {2, 3, 4, 5, 6, 7, 8, 9, 10, 11};

  Data d(A, b, c, u, l, B, N, 2, 12);

  MatrixXd B_inicial(2, 2);
  B_inicial.col(0) = A.col(0);
  B_inicial.col(1) = A.col(1);
  Eigen::SparseMatrix<double> B_sparse = B_inicial.sparseView();

  GS g(B_sparse, 2);

  Simplex s(d, g);

  s.findInitialSolution();
  pair<int, int> variable = s.chooseEnteringVariable();
  cout << variable.first << endl;
  cout << variable.second << endl;

  return 0;
}