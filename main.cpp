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
  string filename = argv[1];
  string fo = argv[2];
  // leitor de instâncias mps
  mpsReader mps(filename);

  cout << "leitor" << endl;

  // Matriz A esparsa
  Eigen::SparseMatrix<double> A = mps.A.sparseView();

  // classe data armazenar as informações da instância e as matrizes B e N
  Data d(A, mps.b, mps.c, mps.ub, mps.lb, mps.n_rows_eq + mps.n_rows_inq, mps.n_cols + mps.n_rows_inq + mps.n_rows_eq);

  cout << "lb: " << d.l.transpose() << endl;
  cout << "ub: " << d.u.transpose() << endl;
  cout << "c: " << d.c.transpose() << endl;
  cout << "A: " << endl;
  cout << d.A << endl;
  cout << "b: " << d.b.transpose() << endl;

  if (fo == "min")
    d.c = -d.c;

  // Inicializando B
  Eigen::SparseMatrix<double> B_sparse(d.m, d.m);
  for (size_t i = 0; i < d.B.size(); i++)
  {
    int xi = d.B[i];
    B_sparse.col(i) = d.A.col(xi);
  }

  // classe GS para resolver sistemas lineares com a matriz básica
  GS g(B_sparse, A, d.m);

  // inicializando o simplex
  Simplex s(d, g);
  s.findInitialSolution();
  // verificando se solução inicial é infeasible
  bool phase = s.computeInfeasibility();

  cout << "fase: " << phase << endl;

  // Simplex loop
  while (true)
  {
    // escolhendo a variavel que vai entrar na base
    pair<int, int> variable = s.chooseEnteringVariable(phase);
    cout << "variavel de entrada " << variable.first << " t_sign: " << variable.second << endl;

    // caso nenhuma variável aumente o custo (problema de maximazação) a solução é otima
    if (variable.first == -1)
    {
      cout << "Optimal: " << s.objectiveFunction() << endl;
      return 0;
    }

    // escolhendo a variável de saída
    pair<int, double> leavingVariable = s.chooseLeavingVariable(variable, phase);

    // caso a variável de saída não possua limitante, solução unbouded
    if (leavingVariable.second >= pInf)
    {
      cout << "Unbounded" << endl;
      return 0;
    }

    cout << "variavel de saida " << leavingVariable.first << " max_t: " << leavingVariable.second << endl;

    // atualizando a base
    s.updateBasis(variable, leavingVariable);

    // verificando se a solução continua inviavel
    if (phase)
      phase = s.computeInfeasibility();

    //sleep(1);
  }

  return 0;
}