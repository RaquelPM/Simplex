#include <iostream>
#include <stdio.h>
#include <cstdlib>
#include <vector>
#include <algorithm>
#include <map>
#include <limits>
#include <unistd.h>
#include <fstream>
#include <string>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

#include "Simplex.h"
#include "GS.h"
#include "Data.h"
#include "mpsReader.h"
#include "Scaling.h"

using namespace std;

double pInf = numeric_limits<double>::infinity();
double nInf = -numeric_limits<double>::infinity();
double EPSILON_1 = 1e-5;

int main(int argc, char **argv)
{
  string filename = argv[1];
  string fo = argv[2];
  int pp = atoi(argv[3]);
  int refactor = atoi(argv[4]);

  // scaling para normalizar a matriz A
  Scaling sa;

  // variaveis para armazenar as informações da instância
  MatrixXd A_dense;
  VectorXd b;
  VectorXd l;
  VectorXd u;
  VectorXd c;

  int m, n;

  // leitor de instâncias mps
  mpsReader mps;

  if (fo != "mps")
  {
    ifstream readFile(filename);
    readFile >> m >> n;

    l = VectorXd::Zero(n);
    u = VectorXd::Zero(n);
    A_dense = MatrixXd::Zero(m, n);
    b = VectorXd::Zero(m);
    c = VectorXd::Zero(n);

    readFile.ignore(numeric_limits<streamsize>::max(), '\n');

    for (int i = 0; i < m; i++)
    {
      for (int j = 0; j < n; j++)
      {
        readFile >> A_dense(i, j);
      }
      readFile.ignore(numeric_limits<streamsize>::max(), '\n');
    }

    for (int i = 0; i < n; i++)
    {
      readFile >> c(i);
    }
    readFile.ignore(numeric_limits<streamsize>::max(), '\n');

    string str;

    for (int i = 0; i < n; i++)
    {
      readFile >> str;
      if (!str.compare("inf"))
        l(i) = pInf;
      else if (!str.compare("-inf"))
        l(i) = nInf;
      else
        l(i) = stof(str);
    }

    readFile.ignore(numeric_limits<streamsize>::max(), '\n');

    for (int i = 0; i < n; i++)
    {
      readFile >> str;
      if (!str.compare("inf"))
        u(i) = pInf;
      else if (!str.compare("-inf"))
        u(i) = nInf;
      else
        u(i) = stof(str);
    }

    readFile.ignore(numeric_limits<streamsize>::max(), '\n');

    for (int i = 0; i < m; i++)
    {
      readFile >> b(i);
    }
  }
  else
  {
    mps.read(filename, pp);

    l = mps.lb;
    u = mps.ub;
    A_dense = mps.A;
    b = mps.b;
    c = mps.c;
    m = mps.n_rows_eq + mps.n_rows_inq;
    n = mps.n_cols + mps.n_rows_inq + mps.n_rows_eq;
  }

  // Matriz A esparsa
  Eigen::SparseMatrix<double> A = A_dense.sparseView();

  // classe data armazenar as informações da instância e as matrizes B e N
  Data d(A, b, c, u, l, m, n);

  if (fo == "mps")
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
  Simplex s(d, g, refactor);
  s.findInitialSolution();

  // verificando se solução inicial é infeasible
  bool phase = s.computeInfeasibility();

  cout << "fase: " << phase << endl;
  int count = 0;
  // Simplex loop
  while (true)
  {
    count++;
    // escolhendo a variavel que vai entrar na base
    pair<int, int> variable = s.chooseEnteringVariable(phase);

    // cout << "variavel de entrada " << variable.first << " t_sign: " << variable.second << endl;

    // caso nenhuma variável aumente o custo (problema de maximazação) a solução é otima
    if (variable.first == INT_MAX)
    {
      cout << "Optimal: " << s.objectiveFunction() << " interações: " << count << endl;
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

    // cout << "variavel de saida " << leavingVariable.first << " max_t: " << leavingVariable.second << endl;

    // atualizando a base
    s.updateBasis(variable, leavingVariable);

    // verificando se a solução continua inviavel
    if (phase)
      phase = s.computeInfeasibility();

    cout << "cost: " << s.objectiveFunction() << " fase: " << phase << endl;
  }

  return 0;
}