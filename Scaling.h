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

    double compute_min_vector(const VectorXd& v);
    double compute_min_aij();
    pair<double, double> compute_min_max_row_aij(int index);
    pair<double, double> compute_min_max_col_aij(int index);
    pair<double, double> compute_min_max_col_ratio();
    pair<double, double> compute_min_max_row_ratio();
    void geometric_scale(MatrixXd &A, VectorXd &b, VectorXd &c, VectorXd &l, VectorXd &u, int flag);

public:
    Scaling();
    void teste(MatrixXd A);
    void geometric_iterate(MatrixXd &A, VectorXd &b, VectorXd &c, VectorXd &l, VectorXd &u);
};

#endif