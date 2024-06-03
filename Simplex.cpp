#include "Simplex.h"
#include <umfpack.h>

Simplex::Simplex(int n)
{
    this->n = n;

    y = VectorXd::Zero(n);
    d = VectorXd::Zero(n);
}