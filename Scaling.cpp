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
    A_abs = A.cwiseAbs();
    cout << "min matriz: " << compute_min_aij() << endl;
    cout << "max matriz: " << A_abs.maxCoeff() << endl;

    cout << "min max linhas: " << endl;

    pair<double, double> teste;
    for (long int i = 0; i < A_abs.rows(); i++)
    {
        teste = compute_min_max_row_aij(i);
        cout << teste.first << " " << teste.second << endl;

        if (teste.first == 1 and teste.second == 1)
        {
            cout << "i: " << i << endl;
            cout << A_abs.row(i) << endl;
            exit(0);
        }
    }

    cout << "min max colunas: " << endl;

    for (long int i = 0; i < A.cols(); i++)
    {
        teste = compute_min_max_col_aij(i);
        cout << teste.first << " " << teste.second << endl;
    }

    cout << "min max ratio row: min " << compute_min_max_row_ratio().first << " max: " << compute_min_max_row_ratio().second << endl;
    cout << "min max ratio col: min " << compute_min_max_col_ratio().first << " max: " << compute_min_max_col_ratio().second << endl;

    // geometric_scale(A, true);

    // geometric_iterate(A);

    // cout << A << endl;
}

double Scaling::compute_min_aij()
{
    double min_temp, min = numeric_limits<double>::max();

    VectorXd row = A_abs.row(0);

    for (long int i = 0; i < A_abs.rows(); i++)
    {
        row = A_abs.row(i);
        min_temp = compute_min_vector(row);
        if (min > min_temp)
            min = min_temp;
    }

    return min;
}

double Scaling::compute_min_vector(const VectorXd &v)
{
    double min = numeric_limits<double>::max();

    for (long int i = 0; i < v.size(); i++)
    {
        if (v(i) > EPSILON_1 && min > v(i))
        {
            min = v(i);
        }
    }

    return min;
}

pair<double, double> Scaling::compute_min_max_row_aij(int index)
{
    VectorXd row = A_abs.row(index);

    pair<double, double> min_max;

    min_max.first = compute_min_vector(row);
    min_max.second = row.maxCoeff();

    return min_max;
}

pair<double, double> Scaling::compute_min_max_row_ratio()
{
    pair<double, double> min_max_ratio;
    min_max_ratio.first = pInf;
    min_max_ratio.second = 0;

    for (long int i = 0; i < A_abs.rows(); i++)
    {
        pair<double, double> min_max = compute_min_max_row_aij(i);
        double ratio = min_max.second / min_max.first;
        min_max_ratio.first = min(min_max_ratio.first, ratio);
        min_max_ratio.second = max(min_max_ratio.second, ratio);
    }

    return min_max_ratio;
}

pair<double, double> Scaling::compute_min_max_col_aij(int index)
{
    VectorXd col = A_abs.col(index);

    pair<double, double> min_max;

    min_max.first = compute_min_vector(col);
    min_max.second = col.maxCoeff();

    return min_max;
}

pair<double, double> Scaling::compute_min_max_col_ratio()
{
    pair<double, double> min_max_ratio;
    min_max_ratio.first = pInf;
    min_max_ratio.second = 0;

    for (long int i = 0; i < A_abs.cols(); i++)
    {
        pair<double, double> min_max = compute_min_max_col_aij(i);
        double ratio = min_max.second / min_max.first;
        min_max_ratio.first = min(min_max_ratio.first, ratio);
        min_max_ratio.second = max(min_max_ratio.second, ratio);
    }

    return min_max_ratio;
}

void Scaling::geometric_scale(MatrixXd &A, VectorXd &b, VectorXd &c, VectorXd &l, VectorXd &u, int flag)
{
    int m = A_abs.rows();
    int n = A.cols();

    pair<double, double> min_max;
    double fac;

    for (int i = 0; i < 2; i++)
    {
        if (i == flag)
        {
            for (int j = 0; j < m; j++)
            {
                min_max = compute_min_max_row_aij(j);
                if (min_max.second == 0)
                    continue;
                fac = 1 / sqrt(min_max.first * min_max.second);
                A.row(j) = A.row(j) * fac;
                b(j) = b(j) * fac;
            }
        }
        else
        {
            for (int j = 0; j < n; j++)
            {
                min_max = compute_min_max_col_aij(j);
                if (min_max.second == 0)
                    continue;
                double r = sqrt(min_max.first * min_max.second);
                fac = 1 / r;
                A.col(j) = A.col(j) * fac;
                c(j) = c(j) * fac;
                l(j) = l(j) * r;
                u(j) = u(j) * r;
            }
        }
    }
}

void Scaling::geometric_iterate(MatrixXd &A, VectorXd &b, VectorXd &c, VectorXd &l, VectorXd &u)
{
    A_abs = A.cwiseAbs();

    double min_A = compute_min_aij();
    double max_A = A_abs.maxCoeff();

    double old_ratio, ratio = max_A / min_A;

    // cout << "ratio antes do pré-processamento: " << ratio << endl;

    ratio = 0;

    pair<double, double> min_max_row_ratio = compute_min_max_row_ratio();
    pair<double, double> min_max_col_ratio = compute_min_max_col_ratio();

    int flag = min_max_row_ratio.second > min_max_col_ratio.second;
    // int flag = 0;

    cout << "flag: " << flag << endl;

    for (int i = 1; i <= 15; i++)
    {
        cout << "iteração: " << i << endl;
        old_ratio = ratio;

        geometric_scale(A, b, c, l, u, flag);

        A_abs = A.cwiseAbs();

        min_A = compute_min_aij();
        max_A = A_abs.maxCoeff();

        ratio = max_A / min_A;

        if (i > 1 && ratio > 0.9 * old_ratio)
        {
            break;
        }
    }

    // cout << "min max A: " << min_A << " " << max_A << endl;
    // cout << "ratio após o pré-processamento: " << ratio << endl;
}