#include "GS.h"

GS::GS(Eigen::SparseMatrix<double> B, int n)
{
    this->B = B;
    this->n = n;

    // LU decomposition
    null = (double *)NULL;

    (void)umfpack_di_symbolic(n, n, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, null, null);
    (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, null, null);
}

GS::~GS()
{
    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);
}

void GS::addEk(pair<int, VectorXd> E_)
{
    Ek.push_back(E_);
    //refatorar();
}

void GS::refatorar()
{
    int k = Ek.size();
    MatrixXd B_linha = B;
    for (int i = 0; i < k; i++)
    {
        Eigen::MatrixXd E = Eigen::MatrixXd::Identity(n, n);
        E.col(Ek[i].first) = Ek[i].second;

        B_linha = B_linha * E;
    }

    B = B_linha.sparseView();

    Ek.clear();

    umfpack_di_free_symbolic(&Symbolic);
    umfpack_di_free_numeric(&Numeric);

    (void)umfpack_di_symbolic(n, n, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), &Symbolic, null, null);
    (void)umfpack_di_numeric(B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), Symbolic, &Numeric, null, null);
}

VectorXd GS::FTRAN(VectorXd a)
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

    VectorXd v(n);

    (void)umfpack_di_solve(UMFPACK_A, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), v.data(), a.data(), Numeric, null, null);

    // solving E_k * d = v_k

    int k = Ek.size();

    for (int i = 0; i < k; i++)
    {
        int eta = Ek[i].first;
        VectorXd E = Ek[i].second;

        double p = v(eta) / E(eta);
        v(eta) = p;

        for (int j = 0; j < n; j++)
        {
            if (j != eta)
                v(j) = v(j) - p * E(j);
        }
    }

    return v;
}

VectorXd GS::BTRAN(VectorXd c)
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

    if (c.isZero())
        return c;

    VectorXd v = c;

    int k = Ek.size();

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

    VectorXd y = VectorXd::Zero(n);

    (void)umfpack_di_solve(UMFPACK_At, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), y.data(), v.data(), Numeric, null, null);

    return y;
}

VectorXd GS::solveInit(VectorXd b)
{
    // solving the initial solution
    // B*x_b = b - N*x_N

    VectorXd y = VectorXd::Zero(n);

    (void)umfpack_di_solve(UMFPACK_A, B.outerIndexPtr(), B.innerIndexPtr(), B.valuePtr(), y.data(), b.data(), Numeric, null, null);

    return y;
}