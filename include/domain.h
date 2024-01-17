#ifndef DOMAIN_H
#define DOMAIN_H

#include <vector>
#include <string>
#include <iostream>
#include <map>
#include "fwd.h"
#include "constitutivefunction.h"
#include "element.h"
#include "node.h"
#include "edge.h"


namespace Oceane {

struct Load {
    std::string _elset;
    Index _face;
    double _load;
};

struct Fixity {
    std::string _elset;
    Index _face;
};



class AbaqusIO;
class Domain
{

public:
    Domain();
    void add_Nodeptr(Nodeptr node);
    void add_Elementptr(Elementptr elem);
    void addElset(std::string name,std::vector<Elementptr> elset );
    void addNset(std::string name, std::vector<Nodeptr> nset);
    void setNsetmap(std::map<std::string,std::vector<Nodeptr> > nsetmap);
    void setElsetmap(std::map<std::string,std::vector<Elementptr> > elsetmap );
    void setLoads(std::vector<Load> loads);
    void setFixity(std::vector<Fixity> fix);
    const std::vector<Elementptr>& getElements();
    const std::vector<Nodeptr>& getNodes();
    size_t getTotalUnknowns();
    size_t getTotalVariables();
    const Nodeptr getNodeptr(Index id) const ;
    const Elementptr getElementptr(Index id) const ;
    const std::vector<Nodeptr> getNodes_at(std::vector<Index> ids) const;
    const constitutiveFnPtr getConstitutiveFunction() const;
    void process();


    void Init();
    void printElements();
    void printNodes();
    void Print();

    virtual void read() = 0;
    virtual void init_variables() = 0;
    virtual void solve() = 0;
    virtual void display() = 0;
    void process(Domain * obj);
private:
    std::vector<Nodeptr> m_nodes;
    std::vector<Elementptr> m_elements;
    std::vector<Elementptr> m_boundaryElements;
    std::vector<Edgeptr> m_edges;
    constitutiveFnPtr m_constitutiveFn;
    std::map<std::string,std::vector<Nodeptr>>  m_Nsetmap;
    std::map<std::string,std::vector<Elementptr>>  m_Elsetmap;
    std::vector<Oceane::Load> m_loads;
    std::vector<Oceane::Fixity> m_fixity;
    std::map <Index,std::vector<Index> > m_nodeConnectivity;


};

} //namespace Oceane
#endif // DOMAIN_H













