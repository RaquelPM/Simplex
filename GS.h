#ifndef GS_H
#define GS_H

#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <vector>
#include <umfpack.h>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

class GS
{
private:
    Eigen::SparseMatrix<double> B;
    vector<pair<int, VectorXd>> Ek;
    int n;

    double *null;
    void *Symbolic, *Numeric;

public:
    GS(Eigen::SparseMatrix<double> B, int n);
    ~GS();

    void addEk(pair<int, VectorXd> E_);
    void refatorar();

    VectorXd FTRAN(VectorXd a);
    VectorXd BTRAN(VectorXd c);
    VectorXd solveInit(VectorXd b);
};

#endif