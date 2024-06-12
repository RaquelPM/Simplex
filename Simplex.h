#ifndef SIMPLEX_H
#define SIMPLEX_H

#include "Data.h"
#include "GS.h"
class Simplex
{
private:
    Data &data;
    GS &gs;
    MatrixXd Nb;
    VectorXd x;
    VectorXd c_B;

public:
    Simplex(Data &data, GS &gs) : data(data), gs(gs) {}
    void findInitialSolution();
    pair<int, int> chooseEnteringVariable();
    pair<int, double> chooseLeavingVariable(pair<int, int> enteringVariable);
    void updateBasis(pair<int, int> enteringVariable, pair<int, double> leavingVariable);
};

#endif