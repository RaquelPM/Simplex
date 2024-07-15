#include <numeric>

#include "Data.h"

Data::Data(int m, int n)
{
    this->m = m;
    this->n = n;

    this->A = Eigen::SparseMatrix<double>(m, n);
    this->b = VectorXd::Zero(m);
    this->c = VectorXd::Zero(n);
    this->u = VectorXd::Zero(n);
    this->l = VectorXd::Zero(n);
    this->B = std::vector<int>(m);
    this->N = std::vector<int>(n);
}

Data::Data(Eigen::SparseMatrix<double> &A, VectorXd &b, VectorXd &c, VectorXd &u, VectorXd &l, int m, int n)
{
    this->m = m;
    this->n = n;

    this->A = A;
    this->b = b;
    this->c = c;
    this->u = u;
    this->l = l;

    // preenchendo B e N
    for (int i = 0; i < n; i++)
    {
        if (i < n - m)
            this->N.push_back(i);
        else
            this->B.push_back(i);
    }
}

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