#include <iostream>
#include <vector>
#include <cmath>
#include <math.h>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <iomanip>
#include <utility>
#include <algorithm>
#include <bits/stdc++.h> 
#include <nlopt.h>
#include <domain.h>

using namespace std;

void reader(vector<vector<double>>& x,string name)
{
    const std::string csvFilePath = name;
    std::ifstream csvFile(csvFilePath);

    if (!csvFile.is_open()) {
        std::cerr << "Error: Unable to open the CSV file." << std::endl;
        return;
    }

    std::vector<std::vector<double>> csvData;
    std::string line;
    bool skipFirstLine = true;
    bool skipFirstField = false;



    while (std::getline(csvFile, line)) {
        //std::cout << line << "\n";
        if (skipFirstLine) {
            skipFirstLine = false;
            continue;
        }

        if(name=="Mem Data.csv")
        {
            //cout<<"Hello2"<<endl;
            skipFirstField = true;
        }

        std::stringstream ss(line);
        std::vector<double> fields;
        std::string field;
        while (std::getline(ss, field, ',')) {
            // cout<<"Hello"<<endl;
            if(field=="")
                break;
            if(skipFirstField){
                skipFirstField = false;
                //cout<<"123"<<endl;
                continue;   
            }
            // cout<<"Hello1"<<endl;
            //cout<<field<<endl;
            fields.push_back(std::stod(field));
        }
        // cout<<"Hello"<<endl;
        x.push_back(fields);
    }
    

    csvFile.close();    
}

void funcname_reader(vector<vector<string>>& x, string name){
    const std::string csvFilePath = name;
    std::ifstream csvFile(csvFilePath);

    if (!csvFile.is_open()) {
        std::cerr << "Error: Unable to open the CSV file." << std::endl;
        return;
    }

    std::vector<std::vector<string>> csvData;
    std::string line;
    
    bool skipFirstLine = true;

    while (std::getline(csvFile, line)) {
        //std::cout << line << "\n";
        if (skipFirstLine) {
            skipFirstLine = false;
            continue;
        }
        //     if(field=="")
        // break;

        std::stringstream ss(line);
        std::vector<string> fields;
        std::string field;
        
        
        if(std::getline(ss, field, ',')){
            // cout << field << endl;
            fields.push_back(field);
            x.push_back(fields);
        }
        if(field == ""){
            break;
        }
    }
    

    csvFile.close();
}

pair<vector<double>, vector<vector<double>>> model(vector<double> K1, vector<double> K2, vector<double> matpar){
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

pair<vector<vector<double>>, vector<vector<double>>> Constitutrel(vector<vector<double>> L, vector<vector<double>> A, vector<vector<string>> funcname, vector<vector<double>> member, int nmemb, vector<vector<double>> forcememb)
{
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
            K1.push_back(sigma[0]);
            K2.push_back(sigma[0]*sigma[0]);
        }
        // K2.push_back(sigma[0]*sigma[0]);
        string namestr = funcname[i][0];
        
        if(namestr == "model"){
            vector<double> temp;
            for(int j = 4; j < member.size(); j++){
                temp.push_back(member[i][j]);
            }
            abcd = model(K1, K2, temp);
        }
        vector<double> beta = abcd.first;
        vector<vector<double>> dbeta = abcd.second;
        
        strain[i][0] = beta[0]+beta[1]*K1[0]+beta[2]*K2[0];
        delta[i][0]=strain[i][0] * L[i][0];
        fdiff[i][0]=(dbeta[0][0]+dbeta[0][1]*2*K1[0]+dbeta[1][0]*K1[0]+dbeta[1][1]*2*K2[0]+beta[1]+dbeta[2][0]*K2[0]+dbeta[2][1]*2*K1[0]*K2[0]+2*beta[2]*K1[0])*L[i][0]/A[i][0];
        
    }
    pair<vector<vector<double>>, vector<vector<double>>> v;
    v.first=delta;
    v.second=fdiff;
    
    return v;
}

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

double objectfunc(vector<double> d, vector<double> grad){
    double result = 0.0;
    vector<double> tempg, dotg;
    vector<vector<double>> tempd, cond, temperr;
    double result = 0.0;
    for(int i = 0; i < d.size(), i++){
        vector<double> temp;
        temp.push_back(d[i]);
        tempd.pushback(temp);
    }
    for(int i = 0; i < nmemb; i++){
        cond.push_back(tempd[i]);
    }
    vector<vector<vector<double>>> con = Constitutrel(L, A, funcname, m_elements, nmemb, cond);
    vector<vector<double>> delta = con[0]; 
    vector<vector<double>> fdiff = con[1];
    vector<vector<double>> errvec = dot(S, tempd);
    for(int i = 0; i < errvec.size(), i++){
        vector<double> temp;
        double ins = errvec[i][0] - multiplier*delta[i][0] - deltasupp[i][0];
        temp.push_back(ins);
        temperr.pushback(temp);
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

int main(){
    ifstream inpuflname("compound_planar.lsx");
    //vector<vector<double>> node; //std::
    vector<Nodeptr> m_nodes;
    //vector<vector<double>> member; //std::
    vector<Elementptr> m_elements;
    vector<vector<string>> funcname; //
    //vector<vector<double>> bc; //std::
    vector<Elementptr> m_boundaryElements;
    //vector<vector<double>> Force; //std::
    vector<Oceane::Load> m_loads;
    //vector<vector<double>> Suppsettle; //std::
    vector<Oceane::Fixity> m_fixity;
    vector<vector<double>> L;
    vector<vector<double>>  A;
    vector<vector<double>> F;
    double Fm;
    double avgforce;
    int nnode;
    int nmemb;

    // Read Node Data
    reader(m_nodes,"Node Data.csv");
    reader(m_elements,"Mem Data.csv");
    reader(m_boundaryElements,"Bound Con.csv");
    reader(m_loads,"Force Data.csv");
    reader(m_fixity,"Supp set.csv");
    funcname_reader(funcname, "Mem Data.csv");

    nnode = m_nodes.size();
    nmemb = m_elements.size();

    for (int i = 0; i < nmemb; i++){
        vector<double> temp;
        temp.push_back(0);
        L.push_back(temp);
    }

    // A.resize(nmemb, vector<double>(1, 0.0));
    for (int i = 0; i < nmemb; i++){
        vector<double> temp;
        temp.push_back(m_elements[i][3]);
        A.push_back(temp);
    }
    // F.resize(Force.size(), vector<double>(1, 0.0));
    for (int i = 0; i < m_loads.size(); i++){
        // vector<int> temp;
        for(int j = 1; j < 3; j++){
            vector<double> temp;
            temp.push_back(m_loads[i][j]);
            F.push_back(temp);
        }
    }

    double Fm = 0;
    for (int i = 0; i < F.size(); i++) {
        double absF = abs(F[i][0]);
        if (absF > Fm) {
            Fm = absF;
        }
    }
    avgforce = Fm;
    
    // double avgforce = Fm;

    // Calculating the length for each member
    for (int i = 0; i < nmemb; i++) {
        double lx = m_nodes[m_elements[i][2]-1][1] - m_nodes[m_elements[i][1]-1][1];
        double mx = m_nodes[m_elements[i][2]-1][2] - m_nodes[m_elements[i][1]-1][2];
        L[i][0] = ceil(sqrt(pow(lx, 2) + pow(mx, 2))*1e8)/1e8;
    }

    // Setting up the equilibrium equations by arranging the corresponding equations for each node by tension coefficient method and assembling them in final master matrix(B)
    vector<vector<double>> B(2 * nnode, vector<double>(nmemb, 0.0));
    for (int i = 0; i < nnode; i++) {
        vector<int> mem;
        vector<int> nodesn;
        for (int j = 0; j < nmemb; j++) {
            for(int k = 1; k < 3; k++){
                if(m_elements[j][k] == i+1){
                    mem.push_back(j);
                    nodesn.push_back(k);
                }
            }
        }
        int s = mem.size();
        // for(int i = 0; i < s; i++){
        //     cout << mem[i] << '\t';
        // }
        // cout << endl;
        // cout << s << endl;
        for (int k : mem) {
            B[2 * i][k] = (-m_nodes[m_elements[k][1] - 1][1] + m_nodes[m_elements[k][2] - 1][1]) / L[k][0];
            B[2 * i + 1][k] = (-m_nodes[m_elements[k][1] - 1][2] + m_nodes[m_elements[k][2] - 1][2]) / L[k][0];
        }
        mem.clear();
        nodesn.clear();
        for (int j = 0; j < nmemb; j++) {
            if (m_elements[j][2] == i + 1) {
                mem.push_back(j);
                nodesn.push_back(0);
            }
        }
        // for(int i = 0; i < mem.size(); i++){
        //     cout << mem[i] << '\t';
        // }
        // cout << endl;
        for (int k : mem) {
            B[2 * i][k] = (m_nodes[m_elements[k][1] - 1][1] - m_nodes[m_elements[k][2] - 1][1]) / L[k][0];
            B[2 * i + 1][k] = (m_nodes[m_elements[k][1] - 1][2] - m_nodes[m_elements[k][2] - 1][2]) / L[k][0];
        }
        // B = -B
    }
    for (int i = 0; i < B.size(); i++) {
        for (int j = 0; j < B[i].size(); j++) {
            B[i][j] = -B[i][j];
        }
    }

    // print_matrix(B);
    vector<int> nodenum;
    vector<int> xyz;
    for (int i = 1; i < 3; i++) {
        for(int j = 0; j < m_boundaryElements.size(); j++){
            if (m_boundaryElements[j][i] == 0){
                nodenum.push_back(i-1);
                xyz.push_back(j);
            }
        }
    }

    int ns = nodenum.size();
    // for(int i=0;i<ns;i++)
    // {
    //     cout<<nodenum[i]<<" "<<endl;
    //     cout<<xyz[i]<<endl;
    // }
    vector<vector<double>> Br(ns, vector<double>(nmemb, 0.0));
    for(int i = 0; i < ns; i++){
        int temp = 2*xyz[i] + nodenum[i];
        for(int j = 0; j < B[0].size(); j++){
            Br[i][j] = B[temp][j];
        }
    }

    // cout << ns << endl;
    // print_matrix(Br);
    vector<vector<double>> Breq;
//     vector<int> temp2;
//     vector<int> temp3;
//     for(int i = 0; i < ns; i++){
//         temp2.push_back(2*xyz[i] + nodenum[i]);
//     }
//     for(int i = 0; i < 2*nnode; i++){
//         vector<int> ::iterator it = find(temp2.begin(), temp2.end(), i);
//         if(it != temp2.end()){
//             temp3.push_back(i)
//         }
//     }
//     for(int i = 0; i < ns, i++){
//         x = temp3[i];
//         Breq.push_back(B[x]);
//     }
    vector<int> temp2;
    for(int i = 0; i < ns; i++){
        temp2.push_back(2*xyz[i] + nodenum[i]);
    }
    sort(temp2.begin(), temp2.end());
    int index = 0;
    for(int i = 0; i < 2*nnode-1; i++){
        if(i != temp2[index]){
            vector<double> temp;
            for(int j = 0; j < B[i].size(); j++){
                temp.push_back(B[i][j]);
            }
            Breq.push_back(temp);
        }
        index += 1;
    }
    // print_matrix(Breq);
    vector<int> qs;
    vector<vector<int>> Suppsettle(m_fixity[0].size(), vector<int>());
    for(int i = 0; i < m_fixity.size(); i++){
        for(int j = 0; j < m_fixity[0].size(); j++){
            Suppsettle[j].push_back(m_fixity[i][j]);
        }
    }
    vector<int> temp4;
    vector<int> temp5;
    for(int i = 0; i < Suppsettle.size(); i++){
        for(int j = 0; j < Suppsettle[0].size(); j++){
            temp4.push_back(Suppsettle[i][j]);
        }
    }
    for(int i = 0; i < ns; i++){
        temp5.push_back(nnode*(nodenum[i] + 1) + xyz[i]);
    }
    for(int i = 0; i < ns; i++){
        qs.push_back(temp4[temp5[i]]);
    }
    // for(int i = 0; i < qs.size(); i++){
    //     cout << qs[i] << '\t';
    // }
    vector<vector<double>> deltasupp;
    for (int i = 0; i < nmemb; i++) {
        int x = 0;
        for (int j = 0; j < ns; j++) {
            x += Br[j][i] * qs[j];
        }
        vector<double> y;
        y.push_back(x);
        deltasupp.push_back(y);
    }
    // print_matrix(deltasupp);
    vector<vector<double>> BC(2 * nnode, std::vector<double>(2 * nnode+nmemb, 0));
    for (int i = 0; i < 2 * nnode; i++) {
        for (int j = 0; j < 2 * nnode; j++) {
            BC[i][j] = B[i][j];
        }
    }
    // print_matrix(BC);
    vector<int> nodepos;
    for (int i = 0; i < ns; i++) {
        nodepos.push_back(2 * xyz[i] + nodenum[i]);
    }
    vector<int> reactpos;
    for (int i = nmemb + 2 * nnode - ns; i < nmemb + 2 * nnode; i++) {
        reactpos.push_back(i);
    }
    for (int i = 0; i < reactpos.size(); i++) {
        BC[(2 * nnode * reactpos[i] + nodepos[i]) % BC.size()][(2 * nnode * reactpos[i] + nodepos[i]) / B.size()] = -1;
    }
    // print_matrix(BC);
    vector<int> sb(2) ;
    sb[0]=Breq.size();
    sb[1]=Breq[0].size();
    std::vector<vector<double>> S(sb[1],(std::vector<double> (sb[1] + sb[0] + ns,0)));
    for(int i = 0; i < sb[1]; i++){
        for(int j = sb[1]; j < sb[1] + sb[0]; j++){
            S[i][j] = Breq[j-sb[1]][i];
        }
    }
    // print_matrix(S);
    vector<vector<double>> d0;
    for(int i = 0; i < sb[0] + sb[1]+ns; i++){
        vector<double> temp;
        temp.push_back(avgforce);
        d0.push_back(temp);
    }
    int max_for_Fm = 0;
    for(int i = 0; i < nmemb; i++){
        if (abs(F[i][0]) > max_for_Fm){
            max_for_Fm = abs(F[i][0]);
        }
    }
    Fm = max_for_Fm;
    avgforce = Fm;

    double numsigdigf = 1 + floor(log10(avgforce));

    vector<vector<double>> f;
    for(int i = 0; i < nmemb; i++){ 
        vector<double> temp;
        temp.push_back(avgforce);
        f.push_back(temp);
    }

    // print_matrix(d0);
    // cout << endl << Fm << endl << avgforce << endl << numsigdigf << endl;
    // print_matrix(f);

    vector<vector<double>> deltacheck = Constitutrel(L, A, funcname, m_elements, nmemb, f).first;

    double de1 = 0;
    for (int i = 0; i < deltacheck.size(); i++) {
        for(int j = 0; j < deltacheck[i].size(); j++){
            de1 = max(de1, abs(deltacheck[i][j]));
        }
        
    }

    double nosigdigd = 1 + floor(log10(de1));
    double multiplier = pow(10, numsigdigf - nosigdigd);
    
    return 0;

    //constraint_eqn
    // the input can be vector or vector<vector>. see virtual methods on how to solve that.

    // double f(const vector<double> &x, vector<double> &grad, void* f_data) {
    //     vector<int> tempx, tempg;
    //     for(int i = 0; i < x.size(), i++){
    //         vector<int> temp;
    //         temp.push_back(x[i]);
    //         tempx.pushback(temp);
    //     }
    //     if(grad.size() > 0){
    //         for(int i = 0; i < x.size(); i++){
    //             vector<vector<int>> main = dot(BC, x)
    //             for(int i = 0; i < main.size(); i++){
    //                 main[i] = main[i] - F[i];
    //             }
    //         }
    //     }
    // }


    double constraint_eqn(vector<double> x, vector<double> grad){
        vector<double> tempg, dotg;
        vector<vector<double>> tempx;
        double result = 0.0;
        for(int i = 0; i < x.size(), i++){
            vector<double> temp;
            temp.push_back(x[i]);
            tempx.pushback(temp);
        }
        vector<vector<double>> main = dot(BC, tempx);
        for(int i = 0; i < main.size(); i++){
            main[i][0] -= F[i][0];
            result += main[i][0]*main[i][0];
        }
        result = sqrt(result);
        if(grad.size() > 0){
            for(int i = 0; i < x.size(); i++){
                double temp = 0.0;
                for(int j = 0; j < BC.size(); j++){
                    temp += BC[j][i] + main[j][0];
                }
                grad[i] = temp;
            }
        }
        return result;
    }
    //nlopt and shit

    // nlopt::opt(NLOPT_LD_SLSQP, d0.size());
    /*
    set_min_objective
    add_equality_constraint
    */
    // void nlopt::opt::set_min_objective(nlopt::vfunc objectfunc, NULL);
    // void nlopt::opt::add_equality_constraint(nlopt::vfunc constraint_eqn, NULL);
    // void nlopt::opt::set_ftol_rel(1e-10);
    // void nlopt::opt::set_maxeval(2000);

    vector<double> d, Fr;
    double fval;
    for(int i = 0; i < d0.size(); i++){
        d.push_back(d0[i][0]);
    }

    nlopt_opt opt;
    opt = nlopt_create(NLOPT_LD_SLSQP, d0.size());
    nlopt_set_min_objective(opt, objectfunc, NULL);
    nlopt_add_equality_constraint(opt, constraint_eqn, &constraint_data, 1e-10);
    nlopt_set_ftol_rel(opt, 1e-10);
    nlopt_set_maxeval(opt, 2000);
    nlopt_optimize(opt, d, fval);

    for(int i = nmemb+2*nnode-ns; i < d.size(), i++){
        Fr.push_back(d[i]);
    }
    int sz = Fr.size();
    // Post Processing
    cout << "Member Forces" << endl << '' << endl;
    for(int i = 0; i < nmemb; i++){
        //elementptr datatype?
    }
    vector<int> xyzs, nodenumber;
    for(int i = 0; i < bc.size(); i++){
        for(int j = 1; j < 3; j++){
            if(bc[i][j] == 1){
                xyz.push_back(i);
                nodenumber.push_back(j);
            }
        }
    }
    double k = 0.0;
    string c;
    cout << '______________________________________________' << endl << 'Node Displacements' << endl << ' ' << endl;
    for(int i = nmemb; i < nmemb + sb[0]; i++){
        if (nodenumber[i] == 0){
            c = 'x';
        }
        elif (nodenumber[i] == 0){
            c = 'y';
        }

        cout << 'u' << c << xyzs[k] + 1 << ' = ' << d[i]/multiplier << endl;
        k += 1
    }
    cout << '______________________________________________' << endl << 'Support Reactions' << endl << ' ' << endl;
    for(int i = 0; i < sz[0]; i++){
        if (nodenumber[i] == 0){
            c = 'x';
        }
        elif (nodenumber[i] == 0){
            c = 'y';
        }

        k = round(Fr[i] * 1e6)/1e6;

        cout << 'F' << c << xyzs[k] + 1 << ' = ' << k << endl;
        k += 1
    }


}


/*
-> Constiturel: Almost done, ending needs to be added..........................................DONE
-> Model, SteelModel...........................................................................NOT NEEDED
-> Check the functioning of dot(MatMul)........................................................DONE
-> Objectfunc: grad[i] needs to be done .......................................................DONE
-> getMatrixA and getMatrixF: Needs to be visited again........................................DONE
-> NLOpt: Need to check the functioning and documentation(Vamsi is working on it)..............DONE
-> NLOpt: opt.optimize is remaining............................................................DONE
-> Printing of the output .....................................................................DONE


-> Testing of the truss2d
-> Incorporate PlaneStress into the framework
-> Change the variables in truss according to the framework
*/