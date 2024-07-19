#include "Simplex.h"

extern double pInf;
extern double nInf;
extern double EPSILON_1;

#include <unistd.h>

void Simplex::findInitialSolution()
{
    // vetor solução
    x.resize(data.n);
    // vetor custos das variáveis básicas
    c_B.resize(data.m);
    // vetor direção
    d.resize(data.m);

    // preenchendo c_B
    for (int i = 0; i < data.m; i++)
    {
        c_B(i) = data.c[data.B[i]];
    }

    VectorXd x_N(data.n - data.m);
    MatrixXd Nb = MatrixXd::Zero(data.m, data.n - data.m);

    for (size_t i = 0; i < data.N.size(); i++)
    {
        int xi = data.N[i];
        Nb.col(i) = data.A.col(xi);
        if (data.u[xi] == pInf && data.l[xi] == nInf)
            x_N[i] = 0;
        else if (data.l[xi] == -pInf)
            x_N[i] = data.u[xi];
        else
            x_N[i] = data.l[xi];
    }

    VectorXd effectXn = Nb * x_N;

    // resolvendo sistema B*x_B = b - N*x_N
    VectorXd x_B = gs.solveInit(data.b - effectXn);

    x << x_N, x_B;
    cout << "x_inicial: " << x.transpose() << endl;

    //exit(0);

    // x << 1, 0, -2, 3, 2, 0, 0, 3, 0, 5, -1, 1;
    // x << 0, 0, 0, 0, 225, 117, 420;
    // x << 0, 0, 2, 3;
}

bool Simplex::computeInfeasibility()
{
    double inf = 0;
    P.clear();
    Q.clear();

    ub_phase = data.u;
    lb_phase = data.l;
    c_phase = VectorXd::Zero(data.n);

    for (int i = 0; i < data.m; i++)
    {
        int xi = data.B[i];
        if (x[xi] < data.l[xi])
        {
            inf += (data.l[xi] - x[xi]);
            P.push_back(xi);
            lb_phase(xi) = nInf;
            ub_phase(xi) = data.l[xi];
            c_phase[xi] = 1;
        }
        else if (x[xi] > data.u[xi])
        {
            inf += (x[xi] - data.u[xi]);
            Q.push_back(xi);
            ub_phase(xi) = pInf;
            lb_phase(xi) = data.u[xi];
            c_phase[xi] = -1;
        }
    }

    // cout << "ub_phase: " << ub_phase.transpose() << endl;
    // cout << "lb_phase: " << lb_phase.transpose() << endl;
    // cout << "c_phase: " << c_phase.transpose() << endl;

    return inf > EPSILON_1;
}

pair<int, int> Simplex::chooseEnteringVariable(bool phase)
{
    pair<int, int> variable(INT_MAX, 0);

    // altualizando c_B
    if(phase){
        for (int i = 0; i < data.m; i++)
        {
            c_B(i) = c_phase[data.B[i]];
        }
    } else {
        for (int i = 0; i < data.m; i++)
        {
            c_B(i) = data.c[data.B[i]];
        }
    }

    // calculando os duais
    VectorXd y = gs.BTRAN(c_B);
    //cout << "duais: " << y.transpose() << endl;

    // calculando o custo reduzido para todas as variáveis
    VectorXd reduced_cost;

    // caso a solução seja viavel
    if (!phase)
        reduced_cost = data.c - data.A.transpose() * y;
    // caso a solução seja inviavel
    else
        reduced_cost = c_phase - data.A.transpose() * y;

    //cout << "reduced_cost: " << reduced_cost.transpose() << endl;

    // escolhendo a variavel de entrada (problema de maximização)
    for (size_t i = 0; i < data.N.size(); i++)
    {
        int xi = data.N[i];
        if (reduced_cost(xi) < -EPSILON_1 && x[xi] - EPSILON_1 > data.l[xi])
        {
            if(xi < variable.first){
                variable.first = xi;
                variable.second = -1;
            }
        }
        else if (reduced_cost(xi) > EPSILON_1 && x[xi] + EPSILON_1 < data.u[xi])
        {
            if(xi < variable.first){
                variable.first = xi;
                variable.second = 1;
            }
        }
    }

    return variable;
}

pair<int, double> Simplex::chooseLeavingVariable(pair<int, int> enteringVariable, bool phase)
{
    pair<int, double> variable(INT_MAX, 0);

    int t_sign = enteringVariable.second;

    // calculando vetor direção
    d = gs.FTRAN(data.A.col(enteringVariable.first));
    //cout << "vetor d: " << d.transpose() << endl;

    VectorXd ub = phase ? ub_phase : data.u;
    VectorXd lb = phase ? lb_phase : data.l;

    VectorXd step(data.m);
    for (size_t i = 0; i < data.B.size(); i++)
    {
        int xi = data.B[i];
        double absD = abs(d(i));
        // caso d(i) seja zero
        if (absD < EPSILON_1)
            step(i) = pInf;
        // caso d(i) e -t_sign tenham sinais iguais
        else if (d(i) * -t_sign > EPSILON_1)
            step(i) = (ub[xi] - x[xi]) / absD;
        // caso d(i) e -t_sign tenham sinais opostos
        else if (d(i) * -t_sign < -EPSILON_1)
        {
            step(i) = (x[xi] - lb[xi]) / absD;
        }
    }

    //cout << "vetor t: " << step.transpose() << endl;

    // cout << endl;

    // // verificando a variavel limitante

    // for(size_t i =0; i < data.B.size(); i++){
    //     cout << data.B[i] << " ";
    // } cout << endl;

    // for(size_t i =0; i < data.m; i++){
    //     cout << step(i) << " ";
    // } cout << endl;

    // cout << endl;

    double min_step = pInf;
    for (size_t i = 0; i < data.B.size(); i++)
    {
        if (step(i) <= min_step)
        {           
            if(step(i) == min_step){
                // cout << "IGUAL!" << endl;
                // cout << step(i) << " " << min_step << endl;
                // cout << data.B[i] << " " << variable.first << endl;
                if(data.B[i] < variable.first){
                    variable.first = data.B[i];
                }
                //sleep(1);
            } else {
                // cout << "DI" << endl;
                variable.first = data.B[i]; 
                //sleep(1);
            }
            min_step = step(i);
        }
    }

    // caso a variável de saída seja a variável de entrada
    if (min_step > ub[enteringVariable.first] - lb[enteringVariable.first])
    {
        variable.first = enteringVariable.first;
        variable.second = ub[enteringVariable.first] - lb[enteringVariable.first];
        return variable;
    }

    variable.second = min_step;

    return variable;
}

void Simplex::updateBasis(pair<int, int> enteringVariable, pair<int, double> leavingVariable)
{
    // atualizando x_B
    int t_signal = enteringVariable.second;
    for (size_t i = 0; i < data.B.size(); i++)
    {
        x[data.B[i]] += leavingVariable.second * -t_signal * d(i);
    }

    // atualizando x_j
    x[enteringVariable.first] += leavingVariable.second * t_signal;

    // caso a variavel de entrada não seja básica
    if (enteringVariable.first == leavingVariable.first)
        return;

    // atualizando a base
    for (size_t i = 0; i < data.B.size(); i++)
    {
        if (data.B[i] == leavingVariable.first)
        {
            data.B[i] = enteringVariable.first;
            pair<int, VectorXd> E_(i, d);
            gs.addEk(E_);
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

    // refatorando
    if(gs.getEkSize() >= 20) gs.refatorar(data.B);
    //gs.refatorar(data.B);

    //cout << "soluçao: " << x.transpose() << endl;
}

double Simplex::objectiveFunction()
{
    return data.c.transpose() * x;
}