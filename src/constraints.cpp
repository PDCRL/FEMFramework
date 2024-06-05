#include "constraints.h"

#include<iostream>
#include<vector>

vector<Oceane::Nodeptr> m_nodes;
vector<Oceane::Elementptr> m_elements;
vector<Oceane::Nodeptr> m_boundaryElements;
vector<vector<double>> m_loads;
vector<vector<double>> m_fixity;
vector<vector<double>> L;
vector<vector<double>>  A;
vector<vector<double>> F;
vector<vector<double>> Breq;
vector<vector<double>> BC;
vector<vector<double>> S;
vector<vector<double>> deltasupp;
vector<vector<double>> parameters;
vector<vector<string>> funcname;
int nnode;
int nmemb;
double Fm;
double avgforce;
vector<int> nodenum;
vector<int> xyz;
int ns;
vector<vector<double>> Br;
vector<int> qs;
vector<int> nodepos,reactpos;
vector<vector<double>> f,deltacheck;
double multiplier;
vector<double> d,Fr;
vector<vector<double>> d0;
double fval;
vector<int> sb(2);
int sz;
// struct ConstraintData {
//     vector<vector<double>> BC;
//     vector<vector<double>> F;
// };

// struct ObjectiveData {
//     vector<vector<double>> L;
//     vector<vector<double>> A;
//     vector<vector<double>> parameters;
//     vector<vector<string>> funcname;
//     int nmemb;
// };