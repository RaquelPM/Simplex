#include "Scaling.h"

#include <cmath>

extern double EPSILON_1;
extern double pInf;
extern double nInf;

Scaling::Scaling()
{
}

void Scaling::teste(MatrixXd A)
{
    cout << "min matriz: " << compute_min_aij(A) << endl;
    cout << "max matriz: " << A.maxCoeff() << endl;

    cout << "min max linhas: " << endl;

    pair<long double, long double> teste;
    for (long int i = 0; i < A.rows(); i++)
    {
        teste = compute_min_max_row_aij(A, i);
        cout << teste.first << " " << teste.second << endl;
    }

    cout << "min max colunas: " << endl;

    for (long int i = 0; i < A.cols(); i++)
    {
        teste = compute_min_max_col_aij(A, i);
        cout << teste.first << " " << teste.second << endl;
    }

    cout << "min max ratio row: min " << compute_min_max_row_ratio(A).first << " max: " << compute_min_max_row_ratio(A).second << endl;
    cout << "min max ratio col: min " << compute_min_max_col_ratio(A).first << " max: " << compute_min_max_col_ratio(A).second << endl;

    // geometric_scale(A, true);

    // geometric_iterate(A);

    cout << A << endl;
}

long double Scaling::compute_min_aij(MatrixXd A)
{
    // MatrixXd A_abs = A.cwiseAbs();
    long double min_temp, min = numeric_limits<long double>::max();

    VectorXd row = A_abs.row(0);

    for (long int i = 0; i < A.rows(); i++)
    {
        row = A_abs.row(i);
        min_temp = compute_min_vector(row);
        if (min > min_temp)
            min = min_temp;
    }

    return min;
}

long double Scaling::compute_min_vector(VectorXd v)
{
    long double min = numeric_limits<long double>::max();

    for (long int i = 0; i < v.size(); i++)
    {
        if (v(i) >= EPSILON_1 && min > v(i))
            min = v(i);
    }

    return min;
}

pair<long double, long double> Scaling::compute_min_max_row_aij(MatrixXd A, int index)
{
    // MatrixXd A_abs = A.cwiseAbs();
    VectorXd row = A_abs.row(index);

    pair<long double, long double> min_max;

    min_max.first = compute_min_vector(row);
    min_max.second = row.maxCoeff();

    return min_max;
}

pair<long double, long double> Scaling::compute_min_max_row_ratio(MatrixXd A)
{
    pair<long double, long double> min_max_ratio;
    min_max_ratio.first = pInf;
    min_max_ratio.second = 0;

    for (long int i = 0; i < A.rows(); i++)
    {
        pair<long double, long double> min_max = compute_min_max_row_aij(A, i);
        long double ratio = min_max.second / min_max.first;
        min_max_ratio.first = min(min_max_ratio.first, ratio);
        min_max_ratio.second = max(min_max_ratio.second, ratio);
    }

    return min_max_ratio;
}

pair<long double, long double> Scaling::compute_min_max_col_aij(MatrixXd A, int index)
{
    // MatrixXd A_abs = A.cwiseAbs();

    VectorXd col = A_abs.col(index);

    pair<long double, long double> min_max;

    min_max.first = compute_min_vector(col);
    min_max.second = col.maxCoeff();

    return min_max;
}

pair<long double, long double> Scaling::compute_min_max_col_ratio(MatrixXd A)
{
    pair<long double, long double> min_max_ratio;
    min_max_ratio.first = pInf;
    min_max_ratio.second = 0;

    for (long int i = 0; i < A.cols(); i++)
    {
        pair<long double, long double> min_max = compute_min_max_col_aij(A, i);
        long double ratio = min_max.second / min_max.first;
        min_max_ratio.first = min(min_max_ratio.first, ratio);
        min_max_ratio.second = max(min_max_ratio.second, ratio);
    }

    return min_max_ratio;
}

void Scaling::geometric_scale(MatrixXd &A, VectorXd &b, VectorXd &c, VectorXd &l, VectorXd &u, bool flag)
{
    int m = A.rows();
    int n = A.cols();

    pair<long double, long double> min_max;
    long double fac;

    for (int i = 0; i < 2; i++)
    {
        if (i == flag)
        {
            for (int j = 0; j < m; j++)
            {
                min_max = compute_min_max_row_aij(A, j);
                fac = 1 / sqrt(min_max.first * min_max.second);
                A.row(j) = A.row(j) * fac;
                b(j) = b(j) * fac;
            }
        }
        else
        {
            for (int j = 0; j < n; j++)
            {
                min_max = compute_min_max_col_aij(A, j);
                fac = 1 / sqrt(min_max.first * min_max.second);
                A.col(j) = A.col(j) * fac;
                c(j) = c(j) * fac;
                l(j) = l(j) * fac;
                u(j) = u(j) * fac;
            }
        }
    }
}

void Scaling::geometric_iterate(MatrixXd &A, VectorXd &b, VectorXd &c, VectorXd &l, VectorXd &u)
{
    A_abs = A.cwiseAbs();

    long double min_A = compute_min_aij(A);
    long double max_A = A_abs.maxCoeff();

    long double old_ratio, ratio = max_A / min_A;

    cout << "ratio antes do pré-processamento: " << ratio << endl;

    ratio = 0;

    pair<long double, long double> min_max_row_ratio = compute_min_max_row_ratio(A);
    pair<long double, long double> min_max_col_ratio = compute_min_max_col_ratio(A);

    bool flag = min_max_row_ratio.second > min_max_col_ratio.second;

    cout << "flag: " << flag << endl;

    for (int i = 1; i <= 15; i++)
    {
        cout << "iteração: " << i << endl;
        old_ratio = ratio;

        geometric_scale(A, b, c, l, u, flag);

        A_abs = A.cwiseAbs();

        min_A = compute_min_aij(A);
        max_A = A_abs.maxCoeff();

        ratio = max_A / min_A;

        if (i > 1 && ratio > 0.9 * old_ratio)
        {
            break;
        }
    }

    cout << "min max A: " << min_A << " " << max_A << endl;
    cout << "ratio após o pré-processamento: " << ratio << endl;
}