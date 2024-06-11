#include "Simplex.h"

extern double pInf;
extern double nInf;

void Simplex::findInitialSolution()
{
    Nb.resize(data.m, data.n - data.m);

    x.resize(data.n);
    c_B.resize(data.m);

    for (size_t i = 0; i < data.B.size(); i++)
    {
        c_B(i) = data.c[data.B[i]];
    }

    // VectorXd x_N(data.n - data.m);

    // for (size_t i = 0; i < data.N.size(); i++)
    // {
    //     int xi = data.N[i];
    //     Nb.col(i) = data.A.col(xi);
    // if (data.u[xi] == pInf && data.l[xi] == nInf)
    //     x_N[i] = 0;
    // else if (data.u[xi] == pInf)
    //     x_N[i] = data.l[xi];
    // else
    //     x_N[i] = data.u[xi];
    // }

    // // VectorXd effectXN = Nb * x_N;
    x << 1, 0, -2, 3, 2, 0, 0, 3, 0, 5, -1, 1;
}

pair<int, int> Simplex::chooseEnteringVariable()
{
    pair<int, int> variable(-1, 0);
    VectorXd y = gs.BTRAN(c_B);

    // calculando o custo reduzido

    VectorXd reduced_cost = data.c - data.A.transpose() * y;

    // escolhendo a variavel de entrada

    for (size_t i = 0; i < data.N.size(); i++)
    {
        int xi = data.N[i];
        if (reduced_cost(i) < 0 && x[xi] > data.l[xi])
        {
            variable.first = data.N[i];
            variable.second = -1;
            break;
        }
        else if (reduced_cost(i) > 0 && x[xi] < data.u[xi])
        {
            variable.first = data.N[i];
            variable.second = 1;
            break;
        }
    }

    return variable;
}