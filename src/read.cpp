#include "read.h"
#include "constraints.h"
#include <fstream>
#include <sstream>
#include <iostream>

Read::Read() {};

void Read::node_reader(vector<Oceane::Nodeptr>& x, string name)
{

    const std::string csvFilePath = name;
    std::ifstream csvFile(csvFilePath);
    bool skipFirstLine = true;

    std::string line;
    while (std::getline(csvFile, line)) {
        if (skipFirstLine) {
            skipFirstLine = false;
            continue;
        }
        std::stringstream ss(line);
        // Oceane::Nodeptr fields;

        std::string index, x2, y;
        std::getline(ss, index, ',');
        std::getline(ss, x2, ',');
        std::getline(ss, y, ',');

        Oceane::Index id = std::stod(index);
        vector<double> vec;
        vec.push_back(stod(x2));
        vec.push_back(stod(y));
        Oceane::Nodeptr fields = std::make_shared<Oceane::Node>(id, vec);

        // cout<<"Hello"<<endl;
        x.push_back(fields);
    }

    csvFile.close();
}

void Read::element_reader(vector<vector<string>>& x, vector<Oceane::Elementptr>& y, vector<vector<double>>& z, vector<Oceane::Nodeptr> nodes, string name){
    const std::string csvFilePath = name;
    std::ifstream csvFile(csvFilePath);

    if (!csvFile.is_open()) {
        std::cerr << "Error: Unable to open the CSV file." << std::endl;
        return;
    }

    bool skipFirstLine = true;

    // std::vector<std::vector<string>> csvData;
    std::string line;
    

    while (std::getline(csvFile, line)) {
        if (skipFirstLine) {
            skipFirstLine = false;
            continue;
        }
        std::stringstream ss(line);
        std::vector<string> fieldx;
        
        std::vector<double> fieldz;
        std::string field;
        std::string a1, a2, a3, a4;

        std::getline(ss, a1, ',');
        if(a1 == ""){
            break;
        }
        fieldx.push_back(a1);
        x.push_back(fieldx);
        std::getline(ss, a2, ',');
        Oceane::Index id = std::stod(a2);
        std::getline(ss, a3, ',');
        std::getline(ss, a4, ',');
        Oceane::Nodeptr n1, n2;
        for(int i = 0; i < nodes.size(); i++){
            Oceane::Index index = nodes[i]->getIndex();
            if(static_cast<int>(index) == stoi(a3)){
                n1 = nodes[i];
            }
            else if(static_cast<int>(index) == stoi(a4)){
                n2 = nodes[i];
            }
        }
        vector<Oceane::Nodeptr> vec;
        vec.push_back(n1);
        vec.push_back(n2);

        Oceane::Elementptr fieldy = std::make_shared<Oceane::Element>(id, vec);
        y.push_back(fieldy);

        while(std::getline(ss, field, ',')){
            if(field == ""){
                break;
            }
            fieldz.push_back(std::stod(field));
        }
        z.push_back(fieldz);
        
        // if(std::getline(ss, field, ',')){
        //     // cout << field << endl;
        //     fieldx.push_back(field);
        //     x.push_back(fieldx);
        // }

        
    }
    
    csvFile.close();
}

void Read::read(vector<vector<double>>& x,string name)
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

void Read::startreading()
{
    //cout<<"Hello this is Read.cpp"<<endl;
    Read reader1;
    reader1.node_reader(m_nodes,"../newfem/files/Node Data.csv");
    reader1.element_reader(funcname, m_elements, parameters, m_nodes, "../newfem/files/Mem Data.csv");
    reader1.node_reader(m_boundaryElements,"../newfem/files/Bound Con.csv");
    reader1.read(m_loads,"../newfem/files/Force Data.csv");
    reader1.read(m_fixity,"../newfem/files/Supp set.csv");
    // cout << m_loads.size()<<endl;
    // for(int i=0;i<m_nodes.size();i++)
    // {
    //     for(int j=0;j<3;j++)
    //     {
    //         cout<<m_nodes[i][j]<<endl;
    //     }
    // }
}