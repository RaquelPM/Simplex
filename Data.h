#ifndef DATA_H
#define DATA_H

#include <stdio.h>
#include <iostream>
#include <vector>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class Data
{
public:
    Eigen::SparseMatrix<double> A;
    VectorXd b;
    VectorXd c;
    VectorXd u;
    VectorXd l;
    vector<int> B;
    vector<int> N;
    int m;
    int n;

    Data(Eigen::SparseMatrix<double> &A, VectorXd &b, VectorXd &c, VectorXd &u, VectorXd &l, vector<int> B, vector<int> N, int m, int n) : A(A), b(b), c(c), u(u), l(l), B(B), N(N), m(m), n(n){};
    Data(Eigen::SparseMatrix<double> &A, VectorXd &b, VectorXd &c, VectorXd &u, VectorXd &l, int m, int n);
    Data(int m, int n);
    // gerador de matrizes inversiveis B aleatorias (it_max e o grau de aleatoriedade)
    MatrixXd gen_random_non_singular_mat(int n, int it_max);
    // cria matriz E aleatoria
    std::pair<int, VectorXd> gen_random_eta_mat(int n);

private:
};

#endif // DATA_H