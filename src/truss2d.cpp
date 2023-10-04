#include <iostream>
#include <vector>
#include <cmath>
#include <fstream>

using namespace std;

int main() {
    ifstream inpuflname("compound_planar.lsx");
    vector<vector<double>> node;
    vector<vector<double>> member;
    vector<vector<string>> funcname;
    vector<vector<double>> bc;
    vector<vector<double>> Force;
    vector<vector<double>> Suppsettle;
    vector<vector<double>> L;
    vector<vector<double>> A;
    vector<vector<double>> F;
    double Fm;
    double avgforce;
    int nnode;
    int nmemb;

    // Read Node Data
    while (!inpuflname.eof()) {
        vector<double> row;
        double value;
        inpuflname >> value;
        row.push_back(value);
        node.push_back(row);
    }

    // Read Member Data
    while (!inpuflname.eof()) {
        vector<double> row;
        double value;
        inpuflname >> value;
        row.push_back(value);
        member.push_back(row);
    }

    // Read Function Name
    while (!inpuflname.eof()) {
        vector<string> row;
        string value;
        inpuflname >> value;
        row.push_back(value);
        funcname.push_back(row);
    }

    // Read Boundary Conditions
    while (!inpuflname.eof()) {
        vector<double> row;
        double value;
        inpuflname >> value;
        row.push_back(value);
        bc.push_back(row);
    }

    // Read Force Data
    while (!inpuflname.eof()) {
        vector<double> row;
        double value;
        inpuflname >> value;
        row.push_back(value);
        Force.push_back(row);
    }

    // Read Support Settlements
    while (!inpuflname.eof()) {
        vector<double> row;
        double value;
        inpuflname >> value;
        row.push_back(value);
        Suppsettle.push_back(row);
    }

    nnode = node.size();
    nmemb = member.size();

    L.resize(nmemb, vector<double>(1, 0.0));
    A.resize(nmemb, vector<double>(1, 0.0));
    F.resize(Force.size(), vector<double>(1, 0.0));

    for (int i = 0; i < nmemb; i++) {
        double lx = node[member[i][2]-1][1] - node[member[i][1]-1][1];
        double mx = node[member[i][2]-1][2] - node[member[i][1]-1][2];
        L[i][0] = sqrt(pow(lx, 2) + pow(mx, 2));
    }

    vector<vector<double>> B(2 * nnode, vector<double>(nmemb, 0.0));
    for (int i = 0; i < nnode; i++) {
        vector<int> mem;
        vector<int> nodesn;
        for (int j = 0; j < member.size(); j++) {
            if (member[j][1] == i + 1 || member[j][2] == i + 1) {
                mem.push_back(j);
                if (member[j][1] == i + 1) {
                    nodesn.push_back(member[j][2]);
                } else {
                    nodesn.push_back(member[j][1]);
                }
            }
        }
        int s = mem.size();
        for (int k = 0; k < s; k++) {
            B[2 * i][mem[k]] = 1;
            B[2 * i + 1][mem[k]] = 1;
        }
    }

    return 0;
}