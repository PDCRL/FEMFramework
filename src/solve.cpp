#include "solve.h"
#include "objectfunc.h"
#include <bits/stdc++.h> 
#include <nlopt.hpp>

Solve::Solve() {};

void Solve::Solver()
{
    
    nlopt::opt opt(nlopt::LD_SLSQP, d0.size());
    // cout << "Guru " << d0.size() << endl;
    // cout << "guru " << A.size() << "\n";
    struct ObjectiveData data;
    data.L = L;
    data.A = A;
    data.parameters = parameters;
    data.nmemb = nmemb;
    data.funcname = funcname;
    //Objectfunc obj1; 
    opt.set_min_objective(objectfunc, &data);
    
    struct ConstraintData data1;
    data1.BC = BC;
    data1.F = F;
    opt.add_equality_constraint(constraint_eqn, &data1);
    //cout<<"Hello this is Solve3.cpp"<<endl;
    opt.set_xtol_rel(1e-10);
    //cout<<"Hello this is Solve2.cpp"<<endl;
    opt.set_maxeval(2);
    //cout<<"Hello this is Solve1.cpp"<<endl;
    nlopt::result res = opt.optimize(d, fval);
    //cout<<"Hello this is Solve4.cpp"<<endl;
    cout << opt.last_optimize_result() << endl;
    

    // // nlopt_opt opt;
    // // opt = nlopt_create(NLOPT_LD_SLSQP, d0.size());
    // // nlopt_set_min_objective(opt, myfunc, NULL);
    // // nlopt_add_equality_constraint(opt, myconstraint, NULL, 1e-8);
    // // nlopt_set_xtol_rel(opt, 1e-4);
    // // nlopt_set_maxeval(opt, 2000);
    // // cout << "check1\n";
    // // if (nlopt_optimize(opt, d, &fval) < 0) {
    // //     cout << "nlopt failed" << endl;
    // // } else {
    // //     cout << fval << endl;
    // // }
    // // cout << "check2\n";

    for(int i = nmemb+2*nnode-ns; i < d0.size(); i++){
        Fr.push_back(d[i]);
    }
    sz = Fr.size();
    cout<<nmemb<<" "<<nnode<<" "<<ns<<endl;
    // Post Processing
    cout << "Member Forces" << endl << " " << endl;
    for(int i = 0; i < nmemb; i++){
        //elementptr datatype?
        cout << "F" << m_elements[i]->getIndex() << " = " << round(d[i] * 1e6) / 1e6 << endl;
    }
    vector<int> xyzs, nodenumber;
    for(int i = 0; i < BC.size(); i++){
        for(int j = 1; j < 3; j++){
            if(BC[i][j] == 1){
                xyzs.push_back(i);
                nodenumber.push_back(j);
            }
        }
    }
    int k = 0;
    string c;
    cout << "______________________________________________" << endl << "Node Displacements" << endl << " " << endl;
    for(int i = nmemb; i < nmemb + sb[0]; i++){
        if (nodenumber[i] == 0){
            c = 'x';
        }
        else if (nodenumber[i] != 0){
            c = 'y';
        }

        cout << 'u' << c << xyzs[k] + 1 << " = " << d[i]/multiplier << endl;
        k += 1;
    }
    k = 0;
    cout << "______________________________________________" << endl << "Support Reactions" << endl << " " << endl;
    for(int i = 0; i < sz; i++){
        if (nodenumber[i] == 0){
            c = 'x';
        }
        else if (nodenumber[i] == 0){
            c = 'y';
        }

        k = round(Fr[i] * 1e6)/1e6;

        cout << 'F' << c << xyzs[i] + 1 << " = " << k << endl;
        // k += 1;
    }
    
    /*
    -> d0.size()
    -> sb[0]
    -> sz
    */
}

