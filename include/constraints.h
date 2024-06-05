#ifndef CONSTRAINTS_H
#define CONSTRAINTS_H

#include<iostream>
#include<vector>
#include<map>

#include "element.h"

using namespace std;

extern vector<Oceane::Nodeptr> m_nodes;
extern vector<Oceane::Elementptr> m_elements;
extern vector<Oceane::Nodeptr> m_boundaryElements;
extern vector<vector<double>> m_loads;
extern vector<vector<double>> m_fixity;
extern vector<vector<double>> L;
extern vector<vector<double>>  A;
extern vector<vector<double>> F;
extern vector<vector<double>> Breq;
extern vector<vector<double>> BC;
extern vector<vector<double>> S;
extern vector<vector<double>> deltasupp;
extern vector<vector<double>> parameters;
extern vector<vector<string>> funcname;
extern int nnode;
extern int nmemb;
extern double Fm;
extern double avgforce;
extern vector<int> nodenum;
extern vector<int> xyz;
extern int ns;
extern vector<vector<double>> Br;
extern vector<int> qs;
extern vector<int> nodepos,reactpos;
extern vector<vector<double>> f,deltacheck;
extern double multiplier;
extern vector<double> d,Fr;
extern vector<vector<double>> d0;
extern double fval;
extern vector<int> sb;
extern int sz;
struct ConstraintData {
    vector<vector<double>> BC;
    vector<vector<double>> F;
};

struct ObjectiveData {
    vector<vector<double>> L;
    vector<vector<double>> A;
    vector<vector<double>> parameters;
    vector<vector<string>> funcname;
    int nmemb;
};
#endif