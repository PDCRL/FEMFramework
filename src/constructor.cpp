#include "constructor.h"
#include "constraints.h"
#include <fstream>
#include <sstream>
#include <iostream>

Constructor::Constructor() {};

pair<vector<double>, vector<vector<double>>> Constructor::model(vector<double> K1, vector<double> K2, vector<double> matpar){
    vector<double> beta;
    vector<vector<double>> dbeta(3, vector<double>(3, 0.0));
    beta.push_back(matpar[0]*K1[0]);
    beta.push_back(matpar[1]*pow(K2[0],matpar[2]) + matpar[3]);
    beta.push_back(0);

    dbeta[0][0] = matpar[0];
    dbeta[1][1] = matpar[1] * pow(K2[0], matpar[2]-1) * matpar[2];
    pair<vector<double>, vector<vector<double>>> result;
    result.first = beta;
    result.second = dbeta;

    return result;
}

pair<vector<vector<double>>, vector<vector<double>>> Constructor::Constitutrel(vector<vector<double>>& L, vector<vector<double>>& A, vector<vector<string>>& funcname, vector<vector<double>>& parameters, int nmemb, vector<vector<double>>& forcememb)
{
    //cout<<"Hello this is Constructor.cpp"<<endl;
    vector<vector<double>> fdiff, strain, delta;
    for(int i = 0; i < nmemb; i++)
    {
        vector<double> temp;
        temp.push_back(0);
        fdiff.push_back(temp);
        strain.push_back(temp);
        delta.push_back(temp);
    }

    for(int i = 0; i < nmemb; i++)
    {
        pair<vector<double>, vector<vector<double>>> abcd;
        vector<double> sigma;
        for(int j = 0; j < forcememb[0].size(); j++){
            sigma.push_back(forcememb[i][j]/A[i][j]);
            
        }
        vector<double> K1, K2;
        for(int j = 0; j < sigma.size(); j++){
            K1.push_back(sigma[j]);
            K2.push_back(sigma[j]*sigma[j]);
        }
        // K2.push_back(sigma[0]*sigma[0]);
        string namestr = funcname[i][0];
        
        if(namestr == "model"){
            vector<double> temp;
            for(int j = 1; j < parameters[0].size(); j++){
                temp.push_back(parameters[i][j]);
                // cout << parameters[i][j] << " ";
            }
            // cout << endl;
            Constructor con1;
            abcd = con1.model(K1, K2, temp);
        }
        vector<double> beta = abcd.first;
        vector<vector<double>> dbeta = abcd.second;
        // for(int j = 0; j < beta.size(); j++){
        //     cout << beta[j] << endl;
        // }
        
        strain[i][0] = beta[0]+beta[1]*K1[0]+beta[2]*K2[0];
        delta[i][0]=strain[i][0] * L[i][0];
        fdiff[i][0]=(dbeta[0][0]+dbeta[0][1]*2*K1[0]+dbeta[1][0]*K1[0]+dbeta[1][1]*2*K2[0]+beta[1]+dbeta[2][0]*K2[0]+dbeta[2][1]*2*K1[0]*K2[0]+2*beta[2]*K1[0])*L[i][0]/A[i][0];
        
    }
    pair<vector<vector<double>>, vector<vector<double>>> v;
    v.first=delta;
    v.second=fdiff;
    
    return v;
}