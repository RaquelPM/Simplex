#include <numeric>

#include "Data.h"

Eigen::MatrixXd Data::gen_random_non_singular_mat(int n, int it_max)
{

    Eigen::MatrixXd I = Eigen::MatrixXd::Identity(n, n);

    std::vector<int> idx_vec(n);
    std::iota(idx_vec.begin(), idx_vec.end(), 1);

    for (int k = 0; k < it_max; k++)
    {
        int r1 = idx_vec[rand() % n] - 1;
        std::swap(idx_vec.back(), idx_vec[r1]);
        int r2 = idx_vec[rand() % (n - 1)] - 1;

        I.row(r1) += ((double)rand() / RAND_MAX) * (rand() % 2 ? -1 : 1) * I.row(r2);
    }

    return I;
}

std::pair<int, Eigen::VectorXd> Data::gen_random_eta_mat(int n)
{
    int p = rand() % n;
    Eigen::VectorXd d = Eigen::VectorXd::Random(n);

    while (std::abs(d[p]) < 0.000001)
    {
        d = Eigen::VectorXd::Random(n);
    }

    return std::make_pair(p, d);
}