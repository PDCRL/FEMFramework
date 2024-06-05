#ifndef CONSTRUCTOR_H
#define CONSTRUCTOR_H

#include "constraints.h"

#include <vector>
#include <string>

class Constructor {
public:
    Constructor();
    pair<vector<double>, vector<vector<double>>> model(vector<double> K1, vector<double> K2, vector<double> matpar);
    pair<vector<vector<double>>, vector<vector<double>>> Constitutrel(vector<vector<double>>& L, vector<vector<double>>& A, vector<vector<string>>& funcname, vector<vector<double>>& parameters, int nmemb, vector<vector<double>>& forcememb);
};

#endif