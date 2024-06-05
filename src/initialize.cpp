#include "initialize.h"
#include "constraints.h"
#include "constructor.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>

Initialize::Initialize() {};

void Initialize::startInitializing()
{
    //cout<<"Hello this is Initialize.cpp"<<endl;
    nnode = m_nodes.size();
    nmemb = m_elements.size();
    //cout<<nnode<<" "<<nmemb<<endl;
    for (int i = 0; i < nmemb; i++){
        vector<double> temp;
        temp.push_back(0);
        L.push_back(temp);
        //cout<<temp[0]<<endl;
    }

    // A.resize(nmemb, vector<double>(1, 0.0));
    for (int i = 0; i < nmemb; i++){
        vector<double> temp;
        temp.push_back(parameters[i][0]);
        A.push_back(temp);
        //cout<<temp[0]<<endl;
    }
    // F.resize(Force.size(), vector<double>(1, 0.0));
    for (int i = 0; i < m_loads.size(); i++){
        // vector<int> temp;
        for(int j = 1; j < 3; j++){
            vector<double> temp;
            temp.push_back(m_loads[i][j]);
            F.push_back(temp);
            //cout<<temp[0]<<endl;
        }
    }
    //cout<<L.size()<<" "<<A.size()<<" "<<F.size()<<endl;
    Fm = 0;
    for (int i = 0; i < F.size(); i++) {
        double absF = abs(F[i][0]);
        if (absF > Fm) {
            Fm = absF;
        }
    }
    avgforce = Fm;
    //cout<<Fm<<" "<<avgforce<<endl;
    // double avgforce = Fm;

    // Calculating the length for each member
    vector<vector<double>> elements, nodes1, bc;

    for (int i = 0; i < nmemb; i++) {
        vector<double> temp;
        // double ind, x, y;
        Oceane::Index index = m_elements[i]->getIndex();
        temp.push_back(static_cast<int>(index));
        vector<Oceane::Nodeptr> temp1 = m_elements[i]->getNodes();
        Oceane::Nodeptr node1, node2;
        node1 = temp1[0];
        node2 = temp1[1];
        Oceane::Index x, y;
        x = node1->getIndex();
        y = node2->getIndex();
        temp.push_back(static_cast<int>(x));
        temp.push_back(static_cast<int>(y));
        elements.push_back(temp);
    }

    for (int i = 0; i < nnode; i++) {
        vector<double> temp;
        // double ind, x, y;
        Oceane::Index index = m_nodes[i]->getIndex();
        temp.push_back(static_cast<int>(index));
        vector<double> nodes = m_nodes[i]->getCoord();
        temp.push_back(nodes[0]);
        temp.push_back(nodes[1]);
        nodes1.push_back(temp);
    }

    for (int i = 0; i < nnode; i++) {
        vector<double> temp;
        // double ind, x, y;
        Oceane::Index index = m_nodes[i]->getIndex();
        temp.push_back(static_cast<int>(index));
        vector<double> nodes = m_nodes[i]->getCoord();
        temp.push_back(nodes[0]);
        temp.push_back(nodes[1]);
        bc.push_back(temp);
    }



    for (int i = 0; i < nmemb; i++) {
        double lx = nodes1[elements[i][2]-1][1] - nodes1[elements[i][1]-1][1];
        double mx = nodes1[elements[i][2]-1][2] - nodes1[elements[i][1]-1][2];
        L[i][0] = ceil(sqrt(pow(lx, 2) + pow(mx, 2))*1e8)/1e8;
        //cout<<L[i][0]<<endl;
    }

    //cout<<bc.size()<<" "<<elements.size()<<" "<<nodes1.size()<<endl;

        // Setting up the equilibrium equations by arranging the corresponding equations for each node by tension coefficient method and assembling them in final master matrix(B)
    vector<vector<double>> B(2 * nnode, vector<double>(nmemb, 0.0));
    for (int i = 0; i < nnode; i++) {
        vector<int> mem;
        vector<int> nodesn;
        for (int j = 0; j < nmemb; j++) {
            for(int k = 1; k < 3; k++){
                if(elements[j][k] == i+1){
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
            B[2 * i][k] = (-nodes1[elements[k][1] - 1][1] + nodes1[elements[k][2] - 1][1]) / L[k][0];
            B[2 * i + 1][k] = (-nodes1[elements[k][1] - 1][2] + nodes1[elements[k][2] - 1][2]) / L[k][0];
        }
        mem.clear();
        nodesn.clear();
        for (int j = 0; j < nmemb; j++) {
            if (elements[j][2] == i + 1) {
                mem.push_back(j);
                nodesn.push_back(0);
            }
        }
        // for(int i = 0; i < mem.size(); i++){
        //     cout << mem[i] << '\t';
        // }
        // cout << endl;
        for (int k : mem) {
            B[2 * i][k] = (nodes1[elements[k][1] - 1][1] - nodes1[elements[k][2] - 1][1]) / L[k][0];
            B[2 * i + 1][k] = (nodes1[elements[k][1] - 1][2] - nodes1[elements[k][2] - 1][2]) / L[k][0];
        }
        // B = -B
    }
    for (int i = 0; i < B.size(); i++) {
        for (int j = 0; j < B[i].size(); j++) {
            B[i][j] = -B[i][j];
            //cout<<B[i][j]<<endl;
        }
    }

    // print_matrix(B);
    for (int i = 1; i < 3; i++) {
        for(int j = 0; j < m_boundaryElements.size(); j++){
            if (bc[j][i] == 0){
                nodenum.push_back(i-1);
                xyz.push_back(j);
                //cout<<i-1<<" "<<j<<endl;
            }
        }
    }

    ns = nodenum.size();
    //cout<<ns<<endl;
    // for(int i=0;i<ns;i++)
    // {
    //     cout<<nodenum[i]<<" "<<endl;
    //     cout<<xyz[i]<<endl;
    // }
    for(int i=0;i<ns;i++)
    {
        vector<double> temp;
        for(int j=0;j<nmemb;j++)
        {
            temp.push_back(0);
        }
        Br.push_back(temp);
    }
    for(int i = 0; i < ns; i++){
        int temp = 2*xyz[i] + nodenum[i];
        for(int j = 0; j < B[0].size(); j++){
            Br[i][j] = B[temp][j];
            //cout<<Br[i][j]<<endl;
        }
    }

    // cout << ns << endl;
    // print_matrix(Br);
    // vector<vector<double>> Breq;
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
                //cout<<B[i][j]<<endl;
            }
            Breq.push_back(temp);
        }
        index += 1;
    }
    // print_matrix(Breq);
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
        //cout<<temp4[temp5[i]]<<endl;
        qs.push_back(temp4[temp5[i]]);
    }
    // for(int i = 0; i < qs.size(); i++){
    //     cout << qs[i] << '\t';
    // }
    // vector<vector<double>> deltasupp;
    for (int i = 0; i < nmemb; i++) {
        int x = 0;
        for (int j = 0; j < ns; j++) {
            x += Br[j][i] * qs[j];
        }
        vector<double> y;
        y.push_back(x);
        deltasupp.push_back(y);
        //cout<<x<<endl;
    }
    // print_matrix(deltasupp);
    // vector<vector<double>> BC(2 * nnode, std::vector<double>(2 * nnode+nmemb, 0));
    for (int i = 0; i < 2*nnode; i++){
        vector<double> temp;
        for(int j = 0; j < 2*nnode + nmemb; j++){
            temp.push_back(0);
        }
        BC.push_back(temp);
    }
    for (int i = 0; i < 2 * nnode; i++) {
        for (int j = 0; j < nmemb; j++) {
            BC[i][j] = B[i][j];
            //cout<<BC[i][j]<<endl;
        }
    }
    // print_matrix(BC);
    
    for (int i = 0; i < ns; i++) {
        //cout<<2 * xyz[i] + nodenum[i]<<endl;
        nodepos.push_back(2 * xyz[i] + nodenum[i]);
    }
    for (int i = nmemb + 2 * nnode - ns; i < nmemb + 2 * nnode; i++) {
        //cout<<i<<endl;
        reactpos.push_back(i);
    }
    for (int i = 0; i < reactpos.size(); i++) {
        BC[(2 * nnode * reactpos[i] + nodepos[i]) % BC.size()][(2 * nnode * reactpos[i] + nodepos[i]) / B.size()] = -1;
    }
        for (int i = 0; i < 2 * nnode; i++) {
        for (int j = 0; j < 2 * nnode; j++) {
            //BC[i][j] = B[i][j];
            //cout<<BC[i][j]<<endl;
        }
    }
    // print_matrix(BC);
    //vector<int> sb(2) ;
    sb[0]=Breq.size();
    sb[1]=Breq[0].size();
    //cout<<sb[0]<<" "<<sb[1]<<endl;
    // vector<vector<double>> S(sb[1],(vector<double> (sb[1] + sb[0] + ns,0)));
    for (int i = 0; i < sb[1]; i++){
        vector<double> temp;
        for(int j = 0; j < sb[1] + sb[0] + ns; j++){
            temp.push_back(0);
        }
        S.push_back(temp);
    }
    for(int i = 0; i < sb[1]; i++){
        for(int j = sb[1]; j < sb[1] + sb[0]; j++){
            S[i][j] = Breq[j-sb[1]][i];
            
        }
    }
    for (int i = 0; i < sb[1]; i++){
        for(int j = 0; j < sb[1] + sb[0] + ns; j++){
            //cout<<S[i][j]<<endl;
        }
    }
    // print_matrix(S);
    //vector<vector<double>> d0;
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
    //cout<<avgforce<<" "<<numsigdigf<<endl;
    
    for(int i = 0; i < nmemb; i++){ 
        vector<double> temp;
        temp.push_back(avgforce);
        f.push_back(temp);
        // cout<<avgforce<<endl;
    }

    // print_matrix(d0);
    // cout << endl << Fm << endl << avgforce << endl << numsigdigf << endl;
    // print_matrix(f);
    Constructor con;
    deltacheck = con.Constitutrel(L, A, funcname, parameters, nmemb, f).first;

    double de1 = 0;
    for (int i = 0; i < deltacheck.size(); i++) {
        for(int j = 0; j < deltacheck[i].size(); j++){
            // cout<<deltacheck[i][j]<<endl;
            de1 = max(de1, abs(deltacheck[i][j]));
        }
    }
    cout<<de1<<endl;

    double nosigdigd = 1 + floor(log10(de1));
    multiplier = pow(10, numsigdigf - nosigdigd);
    cout<<nosigdigd<<" "<<multiplier<<endl;
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


    double fval;
    for(int i = 0; i < d0.size(); i++){
        d.push_back(d0[i][0]);
    }

    cout << d.size() <<" "<<"Hi"<<d[1]<< endl;
    
}