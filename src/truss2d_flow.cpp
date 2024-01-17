#include <iostream>
#include "domain.h"
#include "abaqusio.h"
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;
using namespace Oceane;

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

int main() {
    ifstream inpuflname("compound_planar.lsx");
    vector<Oceane::Nodeptr> m_nodes;
    vector<Oceane::Elementptr> m_elements;
    vector<vector<string>> funcname;
    vector<Oceane::Elementptr> m_boundaryElements;
    vector<Oceane::Load> m_loads;
    vector<Oceane::Fixity> m_fixity;
    vector<double> L;
    vector<vector<double>> A;
    vector<vector<double>> F;
    double Fm;
    double avgforce;
    int nnode;
    int nmemb;

    // Read Input data
    Oceane::Domain domain;
    Oceane::AbaqusIO("compound_planar.csv",&domain);

    nnode = m_nodes.size();
    nmemb = m_elements.size();

    L.resize(nmemb, 0.0);
    A = Oceane::AbaqusIO::getMatrixA();
    F = Oceane::AbaqusIO::getMatrixF();

    double Fm = 0;
    for (int i = 0; i < sizeof(F) / sizeof(F[0]); i++) {
        double absF = abs(F[i]);
        if (absF > Fm) {
            Fm = absF;
        }
    }
    
    double avgforce = Fm;

    // Calculating the length for each member
    for (int i = 0; i < nmemb; i++) {
        double lx = node[member[i][2]-1][1] - node[member[i][1]-1][1];
        double mx = node[member[i][2]-1][2] - node[member[i][1]-1][2];
        L[i] = sqrt(pow(lx, 2) + pow(mx, 2));
    }

    // Setting up the equilibrium equations by arranging the corresponding equations for each node by tension coefficient method and assembling them in final master matrix(B)
    vector<vector<double>> B(2 * nnode, vector<double>(nmemb, 0.0));
    for (int i = 0; i < nnode; i++) {
        vector<int> mem;
        vector<int> nodesn;
        for (int j = 0; j < member.size(); j++) {
            if (member[j][1] == i + 1 || member[j][2] == i + 1) {
                mem.push_back(j);
                if (member[j][1] == i + 1) {
                    nodesn.push_back(1);
                } else {
                    nodesn.push_back(2);
                }
            }
        }
        int s = mem.size();
        for (int k : mem) {
            B[2 * i][k] = (-node[member[k][1] - 1][1] + node[member[k][2] - 1][1]) / L[k];
            B[2 * i + 1][k] = (-node[member[k][1] - 1][2] + node[member[k][2] - 1][2]) / L[k];
        }
        mem.clear();
        nodesn.clear();
        for (int j = 0; j < member.size(); j++) {
            if (member[j][2] == i + 1) {
                mem.push_back(j);
                nodesn.push_back(0);
            }
        }
        for (int k : mem) {
            B[2 * i][k] = (node[member[k][1] - 1][1] - node[member[k][2] - 1][1]) / L[k];
            B[2 * i + 1][k] = (node[member[k][1] - 1][2] - node[member[k][2] - 1][2]) / L[k];
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
    std::vector<std::vector<int>> Br(ns, std::vector<int>(nmemb, 0));
    for (int i = 0; i < ns; i++) {
        for (int j = 0; j < nmemb; j++) {
            Br[i][j] = B[2 * xyz[i] + nodenum[i] - 1][j];
        }
    }
    // Breq????????
    // qs needs some correction
    
    std::vector<int> qs;
    for (int i = 0; i < ns; i++) {
        qs.push_back(Suppsettle[xyz[i]][nnode * (nodenum[i] + 1) + xyz[i]]);
    }

    std::vector<std::vector<int>> deltasupp(ns, std::vector<int>(1,0));
    for (int i = 0; i < ns; i++) {
        for (int j = 0; j < nmemb; j++) {
            deltasupp[i][0] += Br[j][i] * qs[i];
        }
    }

    std::vector<std::vector<int>> BC(2 * nnode, std::vector<int>(2 * nnode+nmemb, 0));
    for (int i = 0; i < 2 * nnode; i++) {
        for (int j = 0; j < 2 * nnode; j++) {
            BC[i][j] = B[i][j];
        }
    }

    std::vector<int> nodepos;
    for (int i = 0; i < ns; i++) {
        nodepos.push_back(2 * xyz[i] + nodenum[i]);
    }
    std::vector<int> reactpos;
    for (int i = nmemb + 2 * nnode - ns; i < nmemb + 2 * nnode; i++) {
        reactpos.push_back(i);
    }
    for (int i = 0; i < reactpos.size(); i++) {
        for (int j = 0; j < nodepos.size(); j++) {
            BC[(2 * nnode * reactpos[i] + nodepos[j]) % BC.size()][(2 * nnode * reactpos[i] + nodepos[j]) / B.size()] = -1;
        }
    }
    std::vector<std::vector<int>> Breq,BreqT;

    // Assembling the master matrix(S), the system of non-linear equations to be solved
    std::vector<int> sb(2) ;
    sb[0]=Breq.size();
    sb[1]=Breq[0].size();
    std::vector<vector<int>> S(sb[1],(std::vector<int> (2*sb[1]+ns,0)));

    for (int i = 0; i < sb[1]; i++) {
        for (int j = 0; j < sb[0]; j++) {
            BreqT[i][j] = Breq[j][i];
        }
    }

    for (int i = 0; i < sb[1]; i++) 
    {
        for (int j = sb[1]; j < sb[1]+sb[0]; j++)
        {
            S[i][j]=BreqT[i][j-sb[1]];
        }
    }

    // For initial guess d0
    std::vector<std::vector<int>> d0(sb[0]+sb[1]+ns,(std::vector<int> (1,avgforce)));

    double numsigdigf = 1 + floor(log10(avgforce));

    vector<vector<double>> f(nmemb, vector<double>(1, avgforce));

    double deltacheck = Constitutrel(L, A, funcname, member, nmemb, f)[0];

    double de1 = 0;
    for (int i = 0; i < deltacheck.size(); i++) {
        de1 = max(de1, abs(deltacheck[i]));
    }


    double nosigdigd = 1 + floor(log10(de1));
    double multiplier = pow(10, numsigdigf - nosigdigd);

    return 0;
}
