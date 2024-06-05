#ifndef OBJECTFUNC_H
#define OBJECTFUNC_H

#include "constraints.h"
#include "constructor.h"

#include <vector>
#include <string>

vector<vector<double>> dot(vector<vector<double>> a, vector<vector<double>> b);
double constraint_eqn(const vector<double> &x, vector<double> &grad, void *data1);
double objectfunc(const vector<double> &d, vector<double> &grad, void *data1);   


// class Objectfunc {
// public:
//     Objectfunc();
//     vector<vector<double>> dot(vector<vector<double>> a, vector<vector<double>> b);
//     double constraint_eqn(const vector<double> &x, vector<double> &grad, void *data1);
//     double objectfunc(const vector<double> &d, vector<double> &grad, void *data1);   
// };

#endif