#include "Scaling.h"

#include <cmath>

extern double EPSILON_1;
extern double pInf;
extern double nInf;

Scaling::Scaling(){
    
}

void Scaling::teste(MatrixXd A){
    cout << "min matriz: " << compute_min_aij(A) << endl;
    cout << "max matriz: " << A.maxCoeff() << endl;

    cout << "min max linhas: " << endl;

    pair<double, double> teste;
    for(long int i = 0; i < A.rows(); i++){
        teste = compute_min_max_row_aij(A, i);
        cout << teste.first << " " << teste.second << endl;
    }

    cout << "min max colunas: " << endl;

    for(long int i = 0; i < A.cols(); i++){
        teste = compute_min_max_col_aij(A, i);
        cout << teste.first << " " << teste.second << endl;
    }

    cout << "min max ratio row: min " << compute_min_max_row_ratio(A).first << " max: " << compute_min_max_row_ratio(A).second << endl;
    cout << "min max ratio col: min " << compute_min_max_col_ratio(A).first << " max: " << compute_min_max_col_ratio(A).second << endl;
}


double Scaling::compute_min_aij(MatrixXd A){
    MatrixXd A_abs = A.cwiseAbs();
    double min_temp, min = numeric_limits<double>::max();

    VectorXd row = A_abs.row(0);

    for(long int i = 0; i < A.rows(); i++){
        row = A_abs.row(i);
        min_temp = compute_min_vector(row);
        if(min > min_temp) 
            min = min_temp;
    }

    return min;
}

double Scaling::compute_min_vector(VectorXd v){
    double min = numeric_limits<double>::max();

    for(long int i =0; i < v.size(); i++){
        if(v(i) >= EPSILON_1 && min > v(i)) 
            min = v(i);
    }

    return min;
}

pair<double, double> Scaling::compute_min_max_row_aij(MatrixXd A, int index){
    MatrixXd A_abs = A.cwiseAbs();

    VectorXd row = A_abs.row(index);

    pair<double, double> min_max;

    min_max.first = compute_min_vector(row);
    min_max.second = row.maxCoeff();

    return min_max;
}

pair<double, double> Scaling::compute_min_max_row_ratio(MatrixXd A){
    pair<double, double> min_max_ratio;
    min_max_ratio.first = pInf;
    min_max_ratio.second = 0;

    for(long int i = 0; i < A.rows(); i++){
        pair<double, double> min_max = compute_min_max_row_aij(A, i);
        double ratio = min_max.second/min_max.first;
        min_max_ratio.first = min(min_max_ratio.first, ratio);
        min_max_ratio.second = max(min_max_ratio.second, ratio);
    }

    return min_max_ratio;
}

pair<double, double> Scaling::compute_min_max_col_aij(MatrixXd A, int index){
    MatrixXd A_abs = A.cwiseAbs();

    VectorXd col = A_abs.col(index);

    pair<double, double> min_max;

    min_max.first = compute_min_vector(col);
    min_max.second = col.maxCoeff();

    return min_max;
}


pair<double, double> Scaling::compute_min_max_col_ratio(MatrixXd A){
    pair<double, double> min_max_ratio;
    min_max_ratio.first = pInf;
    min_max_ratio.second = 0;

    for(long int i = 0; i < A.cols(); i++){
        pair<double, double> min_max = compute_min_max_col_aij(A, i);
        double ratio = min_max.second/min_max.first;
        min_max_ratio.first = min(min_max_ratio.first, ratio);
        min_max_ratio.second = max(min_max_ratio.second, ratio);
    }

    return min_max_ratio;
}


void Scaling::geometric_scaling(MatrixXd A, VectorXd b, VectorXd c, VectorXd u, VectorXd l, bool flag){
    int m = A.rows();
    int n = A.cols();

    pair<double, double> min_max;
    double fac;

    for(int i =0; i < 2; i++){
        if(flag){
            for(int j =0; j < m; j++){
                min_max = compute_min_max_row_aij(A, j);
                fac = 1 / sqrt(min_max.first * min_max.second);
                A.row(j) = A.row(j) * fac;
                b(i) = b(i)* fac;
            }
        }
        else {
            for(int j =0; j < n; j++){
                min_max = compute_min_max_col_aij(A, j);
                fac = 1 / sqrt(min_max.first * min_max.second);
                A.col(j) = A.col(j) * fac;
                c(i) = c(i)*fac;
                l(i) = l(i)*fac;
                u(i) = u(i)*fac;
            }
        }
    }
}