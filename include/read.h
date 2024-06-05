#ifndef READ_H
#define READ_H

#include "constraints.h"

#include <vector>
#include <string>

class Read {
public:
    Read();

    void node_reader(vector<Oceane::Nodeptr>& x, string name);
    void element_reader(vector<vector<string>>& x, vector<Oceane::Elementptr>& y, vector<vector<double>>& z, vector<Oceane::Nodeptr> nodes, string name);
    void read(vector<vector<double>>& x,string name);
    void startreading();

private:
    std::string filename;
};

#endif // READ_H
