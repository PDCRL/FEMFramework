#ifndef ABAQUSIO_H
#define ABAQUSIO_H

#include <fstream>
#include <map>
#include <vector>

#include "domain.h"
#include "meshio.h"
namespace Oceane {


using NodeMap = std::map<Index,std::vector<double> > ;
using ElementMap = std::map<std::string,std::map<Oceane::Index,std::vector<Index> > >;
using NsetMap = std::map<std::string,std::vector<Index> >;
using ElsetMap = std::map<std::string,std::vector<Index> >;
using LoadVec = std::vector<Load>;
using FixityVec = std::vector<Fixity>;
using LoadCsv = std::map<Index, std::vector<double>>;
using FixityCsv = std::map<Index, std::vector<double>>;
using BoundCondn = std::map<Index, std::vector<double>>;

class AbaqusIO:public MeshIO
{
public:
    AbaqusIO(std::string szFile,Domain* domain);
    void read();
    void read_csv();

    //TEST FUNCTIONS
    std::map<Index,std::vector<double> > Nodes() {return m_nodemap;}
    ElementMap Elements() {return m_elementmap;}


private:
    void readNodes();
    void readElements( std::string);
    void readNset(std::string);
    void readElset(std::string);
    void readFixity();
    void readLoad();
    void readBoundary();
    void addToDomain();


    std::fstream    m_file;
    std::string prefix_file_name;
    Domain*         m_domain;
    NodeMap m_nodemap;
    ElementMap m_elementmap;
    ElsetMap  m_elsetmap;
    NsetMap   m_nsetmap;
    LoadVec   m_loads;
    FixityVec m_fixity;
    LoadCsv m_loads_csv;
    FixityCsv m_fixity_csv;
    BoundCondn m_bcmap;
};

} // namespace Oceane

#endif // ABAQUSIO_H
