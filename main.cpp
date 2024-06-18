#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <map>
#include <limits>
#include <unistd.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

#include "Simplex.h"
#include "GS.h"
#include "Data.h"
#include "mpsReader.h"

using namespace std;

double pInf = numeric_limits<double>::infinity();
double nInf = -numeric_limits<double>::infinity();
double EPSILON_1 = 1e-5;

int main(int argc, char **argv)
{
  mpsReader mps("SIMPLE.QPS");
  // int seed = std::atoi(argv[2]);
  // srand(seed);
  // int n = std::atoi(argv[1]);

  /*
    MatrixXd A(2, 4);
    A.row(0) << 2, 1, 1, 0;
    A.row(1) << 1, 3, 0, 1;

    VectorXd b(2);
    b << 2, 3;

    VectorXd c(4);
    c << 1, 1, 0, 0;

    VectorXd u(4);
    u << pInf, pInf, pInf, pInf;

    VectorXd l(4);
    l << 0, 0, 0, 0;

    vector<int> B = {2, 3};
    vector<int> N = {0, 1};

    int m = 2;
    int n = 4;

    MatrixXd B_inicial = MatrixXd::Identity(2, 2);
    Eigen::SparseMatrix<double> B_sparse = B_inicial.sparseView();

    */

  /*
    MatrixXd A(3, 7);
    A.row(0) << 3, 2, 1, 2, 1, 0, 0;
    A.row(1) << 1, 1, 1, 1, 0, 1, 0;
    A.row(2) << 4, 3, 3, 4, 0, 0, 1;

    VectorXd b(3);
    b << 225, 117, 420;

    VectorXd c(7);
    c << 19, 13, 12, 17, 0, 0, 0;

    VectorXd u(7);
    u << pInf, pInf, pInf, pInf, pInf, pInf, pInf;

    VectorXd l(7);
    l << 0, 0, 0, 0, 0, 0, 0;

    vector<int> B = {4, 5, 6};
    vector<int> N = {0, 1, 2, 3};

    int m = 3;
    int n = 7;

    MatrixXd B_inicial = MatrixXd::Identity(3, 3);
    Eigen::SparseMatrix<double> B_sparse = B_inicial.sparseView();
  */

  /*

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

   int m = 2;
   int n = 12;

   MatrixXd B_inicial(2, 2);
   B_inicial.col(0) = A.col(0);
   B_inicial.col(1) = A.col(1);
   Eigen::SparseMatrix<double> B_sparse = B_inicial.sparseView();

   // classe data armazenar as informações da instância e as matrizes B e N
   Data d(A, b, c, u, l, B, N, m, n);
   // classe GS para resolver sistemas lineares com a matriz básica
   GS g(B_sparse, m);

   // inicializando o simplex
   Simplex s(d, g);
   s.findInitialSolution();

   // Simplex loop
   while (true)
   {

     // escolhendo a variavel que vai entrar na base
     pair<int, int> variable = s.chooseEnteringVariable();
     cout << "variavel de entrada " << variable.first << " t_sign: " << variable.second << endl;

     // caso nenhuma variável aumente o custo (problema de maximazação) a solução é otima
     if (variable.first == -1)
     {
       cout << "Optimal: " << s.objectiveFunction() << endl;
       return 0;
     }

     // escolhendo a variável de saída
     pair<int, double> leavingVariable = s.chooseLeavingVariable(variable);

     // caso a variável de saída não possua limitante, solução unbouded
     if (leavingVariable.second >= pInf)
     {
       cout << "Unbounded" << endl;
       return 0;
     }

     cout << "variavel de saida " << leavingVariable.first << " max_t: " << leavingVariable.second << endl;

     // atualizando a base
     s.updateBasis(variable, leavingVariable);

     sleep(1);
   }
   */

  return 0;
}