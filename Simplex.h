#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <ostream>
#include <utility>
#include <vector>
#include <algorithm>
#include <map>
#include <cstdlib>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

using Eigen::MatrixXd;
using Eigen::VectorXd;

using namespace std;

class Simplex
{
private:
    Eigen::SparseMatrix<double> B;
    VectorXd d;
    VectorXd y;
    int k;
    int n;

    vector<pair<int, VectorXd>> Ek;

    double *null;
    void *Symbolic, *Numeric;
    void FTRAN();
    void BTRAN();

public:
    Simplex(Eigen::SparseMatrix<double> B, int n, vector<pair<int, VectorXd>> Ek);
};
