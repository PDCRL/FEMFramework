#include <iostream>
#include "domain.h"
#include "abaqusio.h"
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
using namespace Oceane;

class Truss: public Domain {
    public:
    vector<vector<string>> funcname;
    vector<vector<double>> bc;
    vector<vector<double>> Loads;
    vector<vector<double>> Fixity;
    vector<vector<double>> BoundaryConds;
    vector<vector<double>> Suppsettle;
    vector<double> L;
    vector<double> A;
    vector<double> F;
    double Fm;
    double avgforce;
    int nnode;
    int nmemb;
    double area_cross_section;
    string inputflname;


    Truss(string filename): Domain() {
        inputflname = filename;
        Fm = 0.0;
        avgforce = 0.0;
        nnode = 0;
        nmemb = 0;
        area_cross_section = 0.0;
    }

    void read() {
        // Read Nodes
        ifstream inpFile(inputflname + "_node_data.csv");
        if (inpFile.is_open()) {    
            cout << "File opened successfully" << endl;
        }

        Index nodeId;
        double x, y;

        string nodeId_s, x_s, y_s;

        string first_line;
        getline(inpFile, first_line); // skip the first line
        while (getline(inpFile, first_line)) {
            getline(inpFile, nodeId_s, ',');
            getline(inpFile, x_s, ',');
            getline(inpFile, y_s);

            nodeId = (Index)stoi(nodeId_s);
            x = stod(x_s);
            y = stod(y_s);

            Nodeptr m_nodes = make_shared<m_nodes>(nodeId, x, y);
            this->add_Nodeptr(m_nodes);
        }

        // Read Elements
        ifstream inpFile(inputflname + "_member_data.csv");
        if (inpFile.is_open()) {
            cout << "File opened successfully" << endl;
        }

        Index nodeId, m1, m2;
        double matparam1, matparam2, matparam3, matparam4;
        map<Index,vector<Index> > elements;
        string material, nodeId_s, m1_s, m2_s, area_s, matparam1_s, matparam2_s, matparam3_s, matparam4_s;
        
        string first_line;
        getline(inpFile, first_line); // skip the first line
        while (getline(inpFile, first_line)) {
            getline(inpFile, material, ',');
            getline(inpFile, nodeId_s, ',');
            getline(inpFile, m1_s, ',');
            getline(inpFile, m2_s, ',');
            getline(inpFile, area_s, ',');
            getline(inpFile, matparam1_s, ',');
            getline(inpFile, matparam2_s, ',');
            getline(inpFile, matparam3_s, ',');
            getline(inpFile, matparam4_s);

            nodeId = (Index)stoi(nodeId_s);
            m1 = (Index) stoi(m1_s);
            m2 = (Index) stoi(m2_s);

            Elementptr element= make_shared<Element>(nodeId, getNodeptr(m1), getNodeptr(m2));
            this->add_Elementptr(element);
            area_cross_section = stod(area_s);
        }

        // material = "model" as per the example
        matparam1 = stod(matparam1_s);
        matparam2 = stod(matparam2_s);
        matparam3 = stod(matparam3_s);
        matparam4 = stod(matparam4_s);

        // Read Loads
        ifstream inpFile(inputflname + "_force_data.csv");
        if (inpFile.is_open()) {
            cout << "File opened successfully" << endl;
        }

        Index forceId;
        double x, y;

        string forceId_s, Fx, Fy;

        string first_line;
        getline(inpFile, first_line); // skip the first line
        while (getline(inpFile, first_line)) {
            getline(inpFile, forceId_s, ',');
            getline(inpFile, Fx, ',');
            getline(inpFile, Fy);

            forceId = (Index) stoi(forceId_s);
            x = stod(Fx);
            y = stod(Fy);

            Loads.push_back(vector<double>{x,y});
        }

        // Read Fixity
        ifstream inpFile(inputflname + "_support.csv");
        if (inpFile.is_open()) {
            cout << "File opened successfully" << endl;
        }

        Index fixityId;
        double x, y;

        string fixityId_s, ux, uy;

        string first_line;
        getline(inpFile, first_line); // skip the first line
        while (getline(inpFile, first_line)) {
            getline(inpFile, fixityId_s, ',');
            getline(inpFile, ux, ',');
            getline(inpFile, uy);

            fixityId = (Index) stoi(fixityId_s);
            x = stod(ux);
            y = stod(uy);

            Fixity.push_back(vector<double>{x,y});
        }

        // Read Boundary Conditions
        ifstream inpFile(inputflname + "_boundary_conditions.csv");
        if (inpFile.is_open()) {
            cout << "File opened successfully" << endl;
        }

        Index bcId;
        double x, y;

        string bcId_s, ux, uy;

        string first_line;
        getline(inpFile, first_line); // skip the first line
        while (getline(inpFile, first_line)) {
            getline(inpFile, bcId_s, ',');
            getline(inpFile, ux, ',');
            getline(inpFile, uy);

            bcId = (Index) stoi(bcId_s);
            x = stod(ux);
            y = stod(uy);

            BoundaryConds.push_back(vector<double>{x,y});
        }
    }

    vector<double> getMatrixA() {
        vector<double> A;
        map<Index, double>::iterator it = area_cross_section.begin(); 
    
        // Iterating over the map using Iterator till map end. 
        while (it != m_area.end()) { 
            A.push_back(it->second);
        }

        return A;
    }

    vector<double> getMatrixF() {
        vector<double> F;
        vector<vector<double>>::iterator it = Loads.begin(); 
    
        // Iterating over Loads. 
        for (const auto& row : Loads) {
            for (const double& value : row) {
                F.push_back(value);
            }
        }
        return F;
    }

    vector<vector<vector<int>>> Constiturel(vector<double> L, vector<double> A, vector<vector<string>> funcname, vector<vector<double>> member, int nmemb, vector<vector<int>> forcememb){
        vector<vector<int>> fdiff, strain, delta;
        for(int i = 0; i < nmemb; i++){
            vector<int> temp;
            temp.push_back(1);
            fdiff.push_back(temp);
            strain.push_back(temp);
            delta.push_back(temp);
        }

        for(int i = 0; i < nmemb; i++){
            int sigma = (forcememb[i][0]*1)/A[i];
            vector<int> K1, K2;
            K1.push_back(sigma);
            K2.push_back(pow(sigma, 2));
            string temp = funcname[i+1][0];
            vector<string> namestr;
            int space;
            for(int i = 0; i < temp.size(), i++){
                if(int(temp[i]) == 32){
                    space = i;
                    break;
                }
            }
            namestr.push_back(temp.substr(0, space));
            namestr.push_back(temp.substr(space + 1, temp.size() - space));
        }
    }

    void init_variables() {
        nnode = getNodes().size();
        nmemb = getElements().size();

        L.resize(nmemb, 0.0);
        A = getMatrixA();
        F = getMatrixF();

        double Fm = 0;
        for (int i = 0; i < sizeof(F) / sizeof(F[0]); i++) {
            double absF = abs(F[i]);
            if (absF > Fm) {
                Fm = absF;
            }
        }
        
        double avgforce = Fm;

        // Calculating the length for each m_elements
        for (int i = 0; i < nmemb; i++) {
            double lx = m_nodes[m_elements[i][2]-1][1] - m_elements[m_elements[i][1]-1][1];
            double mx = m_nodes[m_elements[i][2]-1][2] - m_elements[m_elements[i][1]-1][2];
            L[i] = sqrt(pow(lx, 2) + pow(mx, 2));
        }

        // Setting up the equilibrium equations by arranging the corresponding equations for each m_nodes by tension coefficient method and assembling them in final master matrix(B)
        vector<vector<double>> B(2 * nnode, vector<double>(nmemb, 0.0));
        for (int i = 0; i < nnode; i++) {
            vector<int> mem;
            vector<int> nodesn;
            for (int j = 0; j < m_elements.size(); j++) {
                if (m_elements[j][1] == i + 1 || m_elements[j][2] == i + 1) {
                    mem.push_back(j);
                    if (m_elements[j][1] == i + 1) {
                        nodesn.push_back(1);
                    } else {
                        nodesn.push_back(2);
                    }
                }
            }
            int s = mem.size();
            for (int k : mem) {
                B[2 * i][k] = (-m_nodes[m_elements[k][1] - 1][1] + m_nodes[m_elements[k][2] - 1][1]) / L[k];
                B[2 * i + 1][k] = (-m_nodes[m_elements[k][1] - 1][2] + m_nodes[m_elements[k][2] - 1][2]) / L[k];
            }
            mem.clear();
            nodesn.clear();
            for (int j = 0; j < m_elements.size(); j++) {
                if (m_elements[j][2] == i + 1) {
                    mem.push_back(j);
                    nodesn.push_back(0);
                }
            }
            for (int k : mem) {
                B[2 * i][k] = (m_nodes[m_elements[k][1] - 1][1] - m_nodes[m_elements[k][2] - 1][1]) / L[k];
                B[2 * i + 1][k] = (m_nodes[m_elements[k][1] - 1][2] - m_nodes[m_elements[k][2] - 1][2]) / L[k];
            }
        }
        // B = -B
        for (int i = 0; i < B.size(); i++) {
            for (int j = 0; j < B[i].size(); j++) {
                B[i][j] = -B[i][j];
            }
        }

        // Finding positions of restraints in movement in different direction for all nodes, storing respective equations in Br and removing them from the main system of equations to be solved.
        vector<int> nodenum;
        vector<int> xyz;
        for (int j = 0; j < bc.size(); j++) {
            if (bc[1][j] == 0|| bc[2][j] == 0) {
                xyz.push_back(j);
                if (bc[1][j] == 0) {
                    nodenum.push_back(1);
                } else {
                    nodenum.push_back(2);
                }
            }
        }

        int ns = nodenum.size();
        vector<vector<int>> Br;//(ns, vector<int>(nmemb, 0));
        vector<int> temp1;
        for (int i = 0; i < ns; i++) {
            int x = 2*xyz[i] + nodenum[i];
            temp1.push_back(x);
        }
        for(int i = 0; i < ns, i++){
            int x = temp1[i];
            Br.push_back(B[x]);
        }

        // Breq
        vector<vector<int>> Breq;
        vector<int> temp2;
        vector<int> temp3;
        for(int i = 0; i < ns; i++){
            temp2.push_back(2*xyz[i] + nodenum[i]);
        }
        for(int i = 0; i < 2*nnode; i++){
            vector<int> ::iterator it = find(temp2.begin(), temp2.end(), i);
            if(it != temp2.end()){
                temp3.push_back(i);
            }
        }
        for(int i = 0; i < ns, i++){
            int x = temp3[i];
            Breq.push_back(B[x]);
        }
        // for (int i = 0; i < ns; i++) {
        //     for (int j = 0; j < nmemb; j++) {
        //         Br[i][j] = B[2 * xyz[i] + nodenum[i] - 1][j];
        //     }
        // }
        
        vector<int> qs;
        vector<vector<int>> Suppsettle(m_fixity[0].size(), vector<int>());
        for(int i = 0; i < m_fixity.size(); i++){
            for(int i = 0; i < m_fixity[0].size(); i++){
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

        //deltasupp
        vector<vector<int>> deltasupp;
        for (int i = 0; i < nmemb; i++) {
            int x = 0;
            for (int j = 0; j < ns; j++) {
                x += Br[j][i] * qs[j];
            }
            vector<int> y;
            y.push_back(x);
            deltasupp.push_back(y);
        }
        // for (int i = 0; i < ns; i++) {
        //     qs.push_back(Suppsettle[xyz[i]][nnode * (nodenum[i] + 1) + xyz[i]]);
        // }

        // vector<vector<int>> deltasupp(ns, vector<int>(1,0));
        // for (int i = 0; i < ns; i++) {
        //     for (int j = 0; j < nmemb; j++) {
        //         deltasupp[i][0] += Br[j][i] * qs[i];
        //     }
        // }

        // BC
        vector<vector<int>> BC(2 * nnode, vector<int>(2 * nnode+nmemb, 0));
        for (int i = 0; i < 2 * nnode; i++) {
            for (int j = 0; j < 2 * nnode; j++) {
                BC[i][j] = B[i][j];
            }
        }

        vector<int> nodepos;
        for (int i = 0; i < ns; i++) {
            nodepos.push_back(2 * xyz[i] + nodenum[i]);
        }
        vector<int> reactpos;
        for (int i = nmemb + 2 * nnode - ns; i < nmemb + 2 * nnode; i++) {
            reactpos.push_back(i);
        }
        for (int i = 0; i < reactpos.size(); i++) {
            for (int j = 0; j < nodepos.size(); j++) {
                BC[(2 * nnode * reactpos[i] + nodepos[j]) % BC.size()][(2 * nnode * reactpos[i] + nodepos[j]) / B.size()] = -1;
            }
        }

        // Assembling the master matrix(S), the system of non-linear equations to be solved
        vector<int> sb(2) ;
        sb[0]=Breq.size();
        sb[1]=Breq[0].size();
        vector<vector<int>> S(sb[1],(vector<int> (2*sb[1]+ns,0)));

        for(int i = 0; i < sb[1]; i++){
            for(int j = sb[1]; j < sb[1] + sb[0]; j++){
                S[i][j] = B[j-sb[1]][i];
            }
        }

        // For initial guess d0
        vector<vector<int>> d0;//(sb[0]+sb[1]+ns,(vector<int> (1,avgforce)));
        for(int i = 0; i < sb[0] + sb[1]+ns; i++){
            vector<int> temp;
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

        for(int i = 0; i < nmemb; i++){ 
            vector<int> temp;
            temp.push_back(avgforce);
            d0.push_back(temp);
        }

        double numsigdigf = 1 + floor(log10(avgforce));

        vector<vector<double>> f(nmemb, vector<double>(1, avgforce));
        
        vector<double> deltacheck = Constitutrel(L, A, funcname, m_elements, nmemb, f)[0];

        double de1 = 0;
        for (int i = 0; i < deltacheck.size(); i++) {
            de1 = max(de1, abs(deltacheck[i]));
        }

        double nosigdigd = 1 + floor(log10(de1));
        double multiplier = pow(10, numsigdigf - nosigdigd);
    }

    void solve() {
        
    }
};