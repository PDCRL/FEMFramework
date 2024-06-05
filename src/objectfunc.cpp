#include "objectfunc.h"
#include <bits/stdc++.h> 
//#include <nlopt.hpp>

// Objectfunc::Objectfunc() {};

// Objectfunc obj2;

vector<vector<double>> dot(vector<vector<double>> a, vector<vector<double>> b){
    int rows = a.size();
    int columns = b[0].size();
    vector<vector<double>> ans(rows, vector<double>(columns, 0.0));
    for(int i = 0; i < rows; i++){
        for(int j = 0; j < columns; j++){
            for(int k = 0; k < a[0].size(); k++){
                ans[i][j] += a[i][k] * b[k][j];
            }
        }
    }

    return ans;
}

double constraint_eqn(const vector<double> &x, vector<double> &grad, void *data1){
        //cout<<"Hi"<<endl;
        vector<vector<double>> tempx;
        vector<double> temp; 
        double sum=0.0;
        for(int i=0;i<x.size();i++)
        {
            temp.push_back(x[i]);
            tempx.push_back(temp);
            cout<<x[i]<<" ";
        }
        cout << "X" << endl;

        for(int i = 0; i < tempx.size(); i++){
            cout << tempx[i][0] << " ";
        }
        cout << endl << "Tempx" << endl;

        vector<vector<double>> tempy = dot(BC,tempx);
        for(int i=0;i<tempy.size();i++)
        {
            tempy[i][0]=tempy[i][0]-F[i][0];
            cout<<tempy[i][0]<<" ";
            sum+=(tempy[i][0]*tempy[i][0]);
        }
        cout << endl << "Tempy" << endl;
        sum=sqrt(sum);
        // cout<<sum<<endl;
        if(grad.size()>0)
        {
            for(int i=0;i<tempx.size();i++)
            {
                vector<vector<double>> bcd;
                vector<double> bcde;
                for(int j=0;j<BC.size();j++)
                {
                    bcde.push_back(2*BC[j][i]);
                    // cout<<2*BC[j][i]<<" ";
                }
                bcd.push_back(bcde);
                // cout<<endl;
                grad[i]=(dot(bcd,tempy)[0][0])/(2*sum);
                // cout<<"HI  " << grad[i]<<endl;
            }
        }
        cout << "___________Grad__________" << endl;
        for(int i = 0; i < grad.size(); i++){
            cout << grad[i] << " ";
        }
        cout << endl;
        //cout << "___________Sum__________" << endl;
        //cout << sum << endl;
        
        
        // vector<double> tempg, dotg;
        // vector<vector<double>> tempx;
        // double result = 0.0;
        // for(int i = 0; i < x.size(); i++){
        //     vector<double> temp;
        //     temp.push_back(x[i]);
        //     tempx.push_back(temp);
        // }
        // struct ConstraintData* data = (struct ConstraintData*)data1;
        // vector<vector<double>> main = dot(data->BC, tempx);
        // for(int i = 0; i < main.size(); i++){
        //     main[i][0] -= data->F[i][0];
        //     result += main[i][0]*main[i][0];
        // }
        // result = sqrt(result);
        // if(grad.size() > 0){
        //     for(int i = 0; i < x.size(); i++){
        //         double temp = 0.0;
        //         for(int j = 0; j < data->BC.size(); j++){
        //             temp += data->BC[j][i] * main[j][0];
        //         }
        //         grad[i] = temp;
        //         cout << grad[i] << " ";
        //     }
        //     cout << endl;
        // }
        // cout<<result<<endl;
        //cout<<"value"<<" "<<sum<<endl;
        return sum;
    }

double objectfunc(const vector<double> &d, vector<double> &grad, void *data1){
    //cout<<"Hello this is Objectfunc.cpp"<<endl;
    double result = 0.0;
    vector<double> tempg, dotg;
    vector<vector<double>> tempd, cond, temperr;
    // double result = 0.0;
    //cout << "Hi" << endl;
    for(int i = 0; i < d.size(); i++){
        vector<double> temp;
        temp.push_back(d[i]);
        tempd.push_back(temp);
    }
    for(int i = 0; i < nmemb; i++){
        cond.push_back(tempd[i]);
    }
    struct ObjectiveData* data = (struct ObjectiveData*)data1;
    Constructor con1;
    pair<vector<vector<double>>, vector<vector<double>>> con = con1.Constitutrel(data->L, data->A, data->funcname, data->parameters, data->nmemb, cond);
    vector<vector<double>> delta = con.first; 
    vector<vector<double>> fdiff = con.second;
    vector<vector<double>> errvec = dot(S, tempd);
    for(int i = 0; i < errvec.size(); i++){
        vector<double> temp;
        double ins = errvec[i][0] - multiplier*delta[i][0] - deltasupp[i][0];
        // Line 292 should be checked.
        temp.push_back(ins);
        temperr.push_back(temp);
        result += ins*ins;
    }
    result = sqrt(result);
    if(grad.size() > 0){
        vector<vector<double>> diag(fdiff.size(), vector<double>(fdiff.size(), 0.0));
        for(int i = 0; i < diag.size(); i++){
            diag[i][i] = -multiplier*fdiff[i][0];
        }
        vector<vector<double>> arr1 = dot(diag, temperr);
        vector<vector<double>> arr2 = dot(Breq, temperr);
        vector<vector<double>> arr3;
        for(int i = 0; i < ns; i++){
            vector<double> temp;
            temp.push_back(1);
            arr3.push_back(temp);
        }
        vector<double> squeeze;
        for(int i = 0; i < arr1.size(); i++){
            squeeze.push_back(arr1[i][0]);
        }

        for(int i = 0; i < arr2.size(); i++){
            squeeze.push_back(arr2[i][0]);
        }

        for(int i = 0; i < arr3.size(); i++){
            squeeze.push_back(arr3[i][0]);
        }

        for(int i = 0; i < d.size(); i++){
            grad[i] = squeeze[i]/result;
        }

    }

    return result; 
}