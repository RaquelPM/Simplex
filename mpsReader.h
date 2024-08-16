/*
 *  mpsReader.h
 *  cppipm
 *
 *  Created by Yiming Yan on 11/07/2014.
 *  Copyright (c) 2014 Yiming Yan. All rights reserved.
 *
 * =================================================================
 *
 * Problem format
 * min 1/2 x'Qx + c'x
 * s.t.
 *      Ax = b
 *      lb <= x <= ub
 *
 * Note that in order to have the above format, slack variables will
 * be added if needed.
 *
 * After call trans2standardForm() function, we get
 * min 1/2 x'Qx + c'x
 * s.t.
 *      Ax = b
 *      x >= 0
 *
 * =================================================================
 * Accepted format: mps, qps, free fromatted mps,
 * free formatted qps
 *
 *
 * In the ROWS section, each row of the constraint matrix must have a
 * row type and a row name specified. The code for indicating row type
 * is as follows:
 *
 *      type        meaning
 * ---------------------------
 *      E           equality
 *      L           less than or equal
 *      G           greater than or equal
 *      N           objective
 *
 * *** N will only be recognised as objective function.
 *
 * RANGES and SOS are not accepted currently.
 *
 * For BOUNDS, we accept only
 *      type            meaning
 *  ---------------------------------------------------
 *      LO              lower bound        lb <= x (< +inf)
 *      UP              upper bound        (0 <=) x <= ub
 *
 * *** Thus lb is always finite.
 *
 * For details about MPS format, see
 *      http://lpsolve.sourceforge.net/5.5/mps-format.htm
 *
 */

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
    mpsReader();
    void read(string fileName, int pre);

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
    vector<int> restricoes;

    int preprocess;

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
    void _splitC(MatrixXd &Araw, VectorXd &c);
};


#endif