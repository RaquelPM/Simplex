#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "Data.h"
#include "GS.h"
class Simplex
{
private:
    Data &data;
    GS &gs;

    VectorXd x;
    VectorXd c_B;
    VectorXd d;

    vector<int> P;
    vector<int> Q;
    VectorXd c_phase;
    VectorXd ub_phase;
    VectorXd lb_phase;

public:
    Simplex(Data &data, GS &gs) : data(data), gs(gs) {}
    void findInitialSolution();
    double computeInfeasibility();
    pair<int, int> chooseEnteringVariable();
    pair<int, double> chooseLeavingVariable(pair<int, int> enteringVariable);
    void updateBasis(pair<int, int> enteringVariable, pair<int, double> leavingVariable);
    double objectiveFunction();
};

#endif