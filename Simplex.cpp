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

pair<int, double> Simplex::chooseLeavingVariable(pair<int, int> enteringVariable)
{
    pair<int, double> variable(-1, 0);

    int t_sign = enteringVariable.second;

    // calculando vetor direção
    VectorXd d = gs.FTRAN(data.A.col(enteringVariable.first));

    VectorXd step(data.m);

    // caso 1: l_b <= x_b - t
    // caso 2: x_b + t <= u_b

    for (size_t i = 0; i < data.B.size(); i++)
    {
        int xi = data.B[i];
        if (d(i) > 0)
        {
            step(i) = (x[xi] - data.l[xi]) / (d(i) * t_sign);
        }
        else if (d(i) < 0)
        {
            step(i) = (x[xi] - data.u[xi]) / (d(i) * t_sign);
        }
        else
        {
            step(i) = pInf;
        }
    }

    // // escolhendo a variavel de saida

    double min_step = pInf;

    for (size_t i = 0; i < data.B.size(); i++)
    {
        if (step(i) < min_step)
        {
            min_step = step(i);
            variable.first = data.B[i];
        }
    }

    min_step = min(min_step, data.u[enteringVariable.first] - data.l[enteringVariable.first]);

    if (min_step == data.u[enteringVariable.first] - data.l[enteringVariable.first])
        variable.first = enteringVariable.first;
    variable.second = min_step;

    return variable;
}

void Simplex::updateBasis(pair<int, int> enteringVariable, pair<int, double> leavingVariable)
{
    // atualizando x_B

    int t_signal = enteringVariable.second;

    for (size_t i = 0; i < data.B.size(); i++)
    {
        x[data.B[i]] += leavingVariable.second * t_signal;
    }

    x[enteringVariable.first] += leavingVariable.second * t_signal;

    // caso a variavel de entrada não seja básica
    if (enteringVariable.first == leavingVariable.first)
        return;

    // atualizando B
    for (size_t i = 0; i < data.B.size(); i++)
    {
        if (data.B[i] == leavingVariable.first)
        {
            data.B[i] = enteringVariable.first;
            break;
        }
    }

    // atualizando N
    for (size_t i = 0; i < data.N.size(); i++)
    {
        if (data.N[i] == enteringVariable.first)
        {
            data.N[i] = leavingVariable.first;
            break;
        }
    }

    // atualizando c_B
    for (size_t i = 0; i < data.B.size(); i++)
    {
        c_B(i) = data.c[data.B[i]];
    }
}