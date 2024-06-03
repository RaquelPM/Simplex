#include "Simplex.h"

#include <umfpack.h>

Simplex::Simplex(Eigen::SparseMatrix<double> B, int n, vector<pair<int, VectorXd>> Ek)
{
    this->B = B;
    this->k = 1;
    this->n = n;
    this->Ek = Ek;

    y = VectorXd::Zero(n);
    d = VectorXd::Zero(n);

    // LU decomposition
    null = (double *)NULL;

    (void)umfpack_di_symbolic(n, n, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, null, null);
    (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, null, null);

    this->FTRAN();
    this->BTRAN();
}

void Simplex::FTRAN()
{
    // solving Bd = a
    // algorithm:
    // B * E_1 * E_2 * ... * E_k * d = a

    // Replace (E_1 * ... * E_k * d) with v_1
    // Solve: B * v_1 = (L * U) * v_1 = a to find v_1
    // Now E_1 * E_2 * ... * E_k * d = v_1
    // Replace E_2 * ... * E_k * d with v_2
    // Solve: E_2 * v_2 = a to find v_2
    // Now E_2 * ... * E_k * d = v_2
    // repeat until Ek * d = v_k

    VectorXd a(n);
    a << 3, 1, 4;
    VectorXd v(n);

    (void)umfpack_di_solve(UMFPACK_A, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), v.data(), a.data(), Numeric, null, null);

    // verificando se de fato Bd = a
    if ((B * v - a).norm() <= 0.0000001)
    {
        std::cout << "Deu certo!" << std::endl;
    }

    // solving E_k * d = v_k

    for (int i = 0; i < k; i++)
    {
        int eta = Ek[i].first;
        VectorXd E = Ek[i].second;
        double p = v(eta) / E(eta);
        v(eta) = p;

        for (int j = 0; j < n; j++)
        {
            if (j == eta)
                continue;
            v(j) = v(j) - p * E(j);
        }
    }

    d = v;

    // verificando se deu certo
    MatrixXd B_linha = B;

    for (int i = 0; i < k; i++)
    {
        Eigen::MatrixXd E = Eigen::MatrixXd::Identity(n, n);
        E.col(Ek[i].first) = Ek[i].second;

        B_linha = B_linha * E;
    }

    if ((B_linha * d - a).norm() <= 0.0000001)
    {
        std::cout << "Deu certo!" << std::endl;
    }

    cout << d << endl;
}

void Simplex::BTRAN()
{
    // solving yB = c
    // algorithm:
    // y * (B * E_1 * E_2 * ... * E_k) = c

    // Replace (y * B * E_1 * E_2 * ... * E_k-1 ) with v_1
    // Solve v_1 * E_k = c
    // Replace y * B * E_1 * E_2 * ... * E_k-2 with v_2
    // Solve v_2 * E_k-1 = v_1
    // ..
    // Repeat until v_k * B = v_k * LU = v_k-1

    VectorXd c(n);
    c << 0, 12, 0;
    VectorXd v = c;

    // solving v_k * E_k = v_k-1

    for (int i = k - 1; i >= 0; i--)
    {
        int eta = Ek[i].first;
        VectorXd E = Ek[i].second;

        for (int j = 0; j < n; j++)
        {
            if (j == eta)
                continue;
            v(eta) -= v(j) * E(j);
        }

        v(eta) = v(eta) / E(eta);
    }

    // solving y * B = v_k

    (void)umfpack_di_solve(UMFPACK_At, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), y.data(), v.data(), Numeric, null, null);

    MatrixXd B_trans = B.transpose();

    // verificando se de fato yB = v_k
    if ((B_trans * y - v).norm() <= 0.0000001)
    {
        std::cout << "Deu certo!" << std::endl;
    }

    // verificando se deu certo
    MatrixXd B_linha = B;

    for (int i = 0; i < k; i++)
    {
        Eigen::MatrixXd E = Eigen::MatrixXd::Identity(n, n);
        E.col(Ek[i].first) = Ek[i].second;

        B_linha = B_linha * E;
    }

    MatrixXd B_linha_trans = B_linha.transpose();

    if ((B_linha_trans * y - c).norm() <= 0.0000001)
    {
        std::cout << "Deu certo!" << std::endl;
    }

    cout << y << endl;
}