#ifndef mpsReader_h
#define mpsReader_h

#include <iostream>
#include <vector>
#include <string>
#include <fstream>

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/src/Core/Matrix.h"

#include "Data.h"

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;

class mpsReader
{
public:
    mpsReader(string fileName);

    string Name;
    int n_rows;
    int n_rows_eq;
    int n_rows_inq;
    int n_cols;

    MatrixXd A;
    VectorXd b;
    VectorXd lb;
    VectorXd ub;
    VectorXd c;

    vector<string> row_labels;
    vector<string> col_labels;
    vector<string> row_list;
    vector<string> col_list;

private:
    long col_pos;
    long rhs_pos;
    long bnd_pos;

    bool bnd_exist;

    void _findPos2Start(ifstream &readFile);
    void _preprocScan(ifstream &readFile);
    int _checkSectionName(string checkWord) const;
    void _nextLine(ifstream &readFile);
    int _getIndex(vector<string> &list, string item) const;
    void _extractData(ifstream &readFile);
    void _getAraw(ifstream &readFile, MatrixXd &Araw);
    void _getbraw(std::ifstream &readFile, VectorXd &braw);
    void _splitRaw(MatrixXd &Araw, VectorXd &braw, VectorXd &c, MatrixXd &A, VectorXd &b);
    void _getBnds(std::ifstream &readFile);
};

#endif