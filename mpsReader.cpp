#include "mpsReader.h"
#include "Scaling.h"

mpsReader::mpsReader(string fileName)
{
    ifstream readFile(fileName);

    if (readFile.is_open())
    {
        // wall_clock timer;
        // timer.tic();

        // get rid of comments or blank line
        _findPos2Start(readFile);

        // get problem dimention
        _preprocScan(readFile);

        // extract data
        _extractData(readFile);

        // readFile.close();
        // time = timer.toc();

        // output
        //_printData();
    }
    else
    {
        cout << "Error: MPSREADER - File not found" << endl;
    }
    readFile.close();
}

mpsReader::mpsReader()
{
}

void mpsReader::read(string fileName, int pre)
{
    ifstream readFile(fileName);
    preprocess = pre;

    if (readFile.is_open())
    {
        // wall_clock timer;
        // timer.tic();

        // get rid of comments or blank line
        _findPos2Start(readFile);

        // get problem dimention
        _preprocScan(readFile);

        // extract data
        _extractData(readFile);

        // readFile.close();
        // time = timer.toc();

        // output
        //_printData();
    }
    else
    {
        cout << "Error: MPSREADER - File not found" << endl;
    }
    readFile.close();
}

void mpsReader::_findPos2Start(ifstream &readFile)
{
    long pos;
    string line;
    string firdWord;
    while (true)
    {
        pos = readFile.tellg();
        getline(readFile, line);
        if (line.empty())
            continue;
        else if (line.find("*") == 0)
            continue;
        else
        {
            readFile.seekg(pos, ios::beg); // go back one line
            break;
        }
    }
}

void mpsReader::_preprocScan(ifstream &readFile)
{
    n_rows = 0;
    n_rows_eq = 0;
    n_rows_inq = 0;
    n_cols = 0;

    bnd_exist = false;

    string tmp = "_";
    string tmpItem = "";

    string firstWord;
    readFile >> firstWord;

    while (!readFile.eof())
    {
        if (!firstWord.empty() && firstWord.find("*") != 0)
        {
            // ======= check termination ======
            if (_checkSectionName(firstWord) == 10)
                break;

            if (_checkSectionName(firstWord) == 8)
            {
                cout << "Error: MPSREASER - Currently cannot handle SOS" << endl;
                break;
            }

            if (_checkSectionName(firstWord) == 9)
            {
                cout << "Error: MPSREASER - Currently cannot handle RANGES" << endl;
                break;
            }
            if (_checkSectionName(firstWord) == 7)
            {
                cout << "Error: MPSREASER - Currently cannot handle OBJSENSE" << endl;
                break;
            }

            if (_checkSectionName(firstWord) == 6)
            {
                cout << "Error: MPSREASER - Currently cannot handle QUADOBJ" << endl;
                break;
            }
            // ======= end check termination ======

            // ======= preprocess sections ======
            if (_checkSectionName(firstWord) == 1) // name
            {
                readFile >> Name;
                _nextLine(readFile);
                readFile >> firstWord;
            }
            else if (_checkSectionName(firstWord) == 2) // rows
            {
                // update firstWord and row_list
                readFile >> firstWord;

                while (_checkSectionName(firstWord) == -1)
                {

                    // count n_rows, n_rows_eq, n_rows_inq
                    /*if (firstWord.compare("E") == 0)
                    {
                        n_rows++;
                        n_rows_eq++;
                    }*/
                    if (firstWord.compare("L") == 0 || firstWord.compare("G") == 0 || firstWord.compare("E") == 0)
                    {
                        n_rows++;
                        n_rows_inq++;
                    }
                    else if (firstWord.compare("N") == 0)
                    {
                        n_rows++;
                    }
                    // store row labels
                    row_labels.push_back(firstWord);

                    // get row_list
                    readFile >> tmpItem;
                    row_list.push_back(tmpItem);

                    _nextLine(readFile);
                    readFile >> firstWord;
                }
            }
            else if (_checkSectionName(firstWord) == 3) // cols
            {
                // get postion
                col_pos = readFile.tellg();

                readFile >> firstWord;

                while (_checkSectionName(firstWord) == -1) // continue if the keyword is not a feild name
                {
                    if (firstWord.compare(tmp) != 0)
                    {
                        // update column list
                        col_list.push_back(firstWord);

                        // count columns
                        n_cols++;
                        tmp = firstWord;
                    }
                    _nextLine(readFile);
                    readFile >> firstWord;
                }
            }
            else if (_checkSectionName(firstWord) == 4) // rhs
            {
                rhs_pos = readFile.tellg();
                _nextLine(readFile);
                readFile >> firstWord;
            }
            else if (_checkSectionName(firstWord) == 5) // bounds
            {
                bnd_exist = true;

                bnd_pos = readFile.tellg();
                _nextLine(readFile);
                readFile >> firstWord;
            }
            else // not a keywod of sections
            {
                _nextLine(readFile);
                readFile >> firstWord;
            }
            // ======= end preprocess sections termination ======
        }
    }
}

void mpsReader::_extractData(ifstream &readFile)
{
    MatrixXd Araw = MatrixXd::Zero(n_rows, n_cols);
    VectorXd braw = VectorXd::Zero(n_rows);

    // get Araw
    _getAraw(readFile, Araw);
    // cout << Araw << endl;

    // get rhs
    _getbraw(readFile, braw);
    // braw.print("braw: ");

    // get c
    c = VectorXd::Zero(n_cols + n_rows_inq);
    _splitC(Araw, c);

    // cout << "braw: " << braw << endl;

    // get bounds

    lb = VectorXd::Zero(n_cols + n_rows_inq);
    ub = VectorXd::Zero(n_cols + n_rows_inq);
    ub.fill(numeric_limits<double>::infinity());

    if (bnd_exist)
        _getBnds(readFile);

    // split Araw to A, Aeq, c
    // and splict braw to b and beq
    A = MatrixXd::Zero(n_rows_inq + n_rows_eq, n_cols + n_rows_inq);
    b = VectorXd::Zero(n_rows_inq + n_rows_eq);
    Scaling sc;
    if(preprocess) sc.geometric_iterate(Araw, braw, c, lb, ub);
    _splitRaw(Araw, braw, c, A, b);
}

void mpsReader::_getAraw(ifstream &readFile, MatrixXd &Araw)
{
    string line, colName, rowName;
    double value;
    int colIdx = 0, rowIdx = 0;

    // go to the the postion of cols
    readFile.seekg(col_pos, ios::beg);
    _nextLine(readFile);
    do
    {
        // clear
        colName = "";
        rowName = "";
        value = 0.0;

        // read one line
        getline(readFile, line);
        istringstream thisLine(line);

        // cout << thisLine.str() << endl;

        thisLine >> colName;

        // break if get to next section
        if (_checkSectionName(colName) != -1)
            break;

        // get col index
        colIdx = _getIndex(col_list, colName);

        while (thisLine >> rowName >> value)
        {
            // get row index
            rowIdx = _getIndex(row_list, rowName);

            Araw(rowIdx, colIdx) = value;
        }

    } while (true);
}

void mpsReader::_getbraw(std::ifstream &readFile, VectorXd &braw)
{
    std::string line, colName, rowName;
    double value;
    int rowIdx = 0;
    // go to the position of rhs
    readFile.seekg(rhs_pos, std::ios::beg);
    _nextLine(readFile);

    do
    {
        // read one line
        getline(readFile, line);
        std::istringstream thisLine(line);

        thisLine >> colName;

        // cout << thisLine.str() << endl;

        if (_checkSectionName(colName) != -1)
            break;

        while (thisLine >> rowName >> value)
        {
            // get row index
            rowIdx = _getIndex(row_list, rowName);

            braw(rowIdx) = value;
        }

    } while (true);
}

void mpsReader::_getBnds(std::ifstream &readFile)
{
    std::string label, colName, rowName;

    double value;

    int colIdx = 0;

    readFile.seekg(bnd_pos, std::ios::beg);
    _nextLine(readFile);

    do
    {
        string line;
        getline(readFile, line);
        istringstream iss(line);
        iss >> label >> rowName >> colName >> value;
        colIdx = _getIndex(col_list, colName);
        // cout << "l: " << label << " " << colName << endl;
        if (label == "LO")
            lb(colIdx) = value;
        else if (label == "UP")
            ub(colIdx) = value;
        else if (label == "FR")
        {
            ub(colIdx) = numeric_limits<double>::infinity();
            lb(colIdx) = -numeric_limits<double>::infinity();
        }
        else if (label == "FX")
        {
            ub(colIdx) = value;
            lb(colIdx) = value;
        }
        {
            if (_checkSectionName(label) == 10)
            {
                std::cout << /*"Error: MPSREADER only accept LO and UP for Bounds"*/ std::endl;
                break;
            }
        }

    } while (true);
}

void mpsReader::_splitC(MatrixXd &Araw, VectorXd &c)
{
    for (int i = n_rows - 1; i >= 0; i--)
    {
        if (row_labels[i] == "N"){
            c.head(n_cols) = Araw.row(i).transpose();
            row_labels[i] = "X";
            Araw.row(i).setZero();
        }
            
    }
}

void mpsReader::_splitRaw(MatrixXd &Araw, VectorXd &braw, VectorXd &c, MatrixXd &A, VectorXd &b)
{
    int counter = 0, counter_inq = 0;

    for (int i = 0; i < n_rows; i++)
    {
        if (row_labels[i] != "X"){
            A.block(counter, 0, 1, n_cols) = Araw.row(i);
            //b.row(counter) = braw.row(i);
            if (row_labels[i] == "L"){
                restricoes.push_back(-1);
                A(counter, n_cols + counter_inq) = -1;
                lb(n_cols + counter_inq) = -numeric_limits<double>::infinity();
                ub(n_cols + counter_inq) = braw(i);
                // cout << "braw(i): " << braw(i) << " " << n_cols + counter_inq << endl;
                counter_inq++;
            }
            else if (row_labels[i] == "G"){
                restricoes.push_back(1);
                A(counter, n_cols + counter_inq) = -1;
                ub(n_cols + counter_inq) = numeric_limits<double>::infinity();
                lb(n_cols + counter_inq) = braw(i);
                // cout << "braw(i): " << braw(i) << endl;
                counter_inq++;
            }
            else if (row_labels[i] == "E"){
                // restricoes.push_back(0);
                // adicionado:
                restricoes.push_back(1);
                A(counter, n_cols + counter_inq) = -1;
                ub(n_cols + counter_inq) = braw(i);
                lb(n_cols + counter_inq) = braw(i);
                // cout << "braw(i): " << braw(i) << endl;
                counter_inq++;
            }
            counter++;
        }
    }
}

int mpsReader::_checkSectionName(string checkWord) const
{

    if (checkWord.compare("NAME") == 0)
        return 1;
    else if (checkWord.compare("ROWS") == 0)
        return 2;
    else if (checkWord.compare("COLUMNS") == 0)
        return 3;
    else if (checkWord.compare("RHS") == 0)
        return 4;
    else if (checkWord.compare("BOUNDS") == 0)
        return 5;
    else if (checkWord.compare("QUADOBJ") == 0)
        return 6;
    else if (checkWord.compare("OBJSENSE") == 0)
        return 7;
    else if (checkWord.compare("SOS") == 0)
        return 8;
    else if (checkWord.compare("RANGES") == 0)
        return 9;
    else if (checkWord.compare("ENDATA") == 0)
        return 10;
    else if (checkWord.compare("FR") == 0)
        return 11;
    else
        return -1;
}

void mpsReader::_nextLine(ifstream &readFile)
{
    readFile.ignore(numeric_limits<streamsize>::max(), '\n');
}
int mpsReader::_getIndex(vector<string> &list, string item) const
{
    int idx = (int)(find(list.begin(), list.end(), item) - list.begin());

    if (idx >= (int)list.size())
        idx = -1;

    return idx;
}