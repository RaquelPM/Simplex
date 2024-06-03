#include <iostream>
#include <stdlib.h>
#include <stdio.h>
#include <numeric>
#include <vector>

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
    int n;
public:
    Simplex(int n);
};
