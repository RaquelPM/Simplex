#ifndef SCALING_H
#define SCALING_H

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

class Scaling
{
private:
    MatrixXd A_abs;

    double compute_min_vector(VectorXd v);
    double compute_min_aij(MatrixXd A);
    pair<double, double> compute_min_max_row_aij(MatrixXd A, int index);
    pair<double, double> compute_min_max_col_aij(MatrixXd A, int index);
    pair<double, double> compute_min_max_col_ratio(MatrixXd A);
    pair<double, double> compute_min_max_row_ratio(MatrixXd A);
    void geometric_scale(MatrixXd A, VectorXd b, VectorXd c, VectorXd u, VectorXd l, bool flag);
public:
    Scaling();
    void teste(MatrixXd A);
};

#endif