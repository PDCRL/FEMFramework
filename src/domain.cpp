#include "domain.h"
#include<iostream>

namespace Oceane {

Domain::Domain()
{
    m_constitutiveFn = std::make_shared<Oceane::PlaneStress>();
}

void Domain::add_Nodeptr(Nodeptr node)
{
    m_nodes.push_back(node);
}

void Domain:: add_Elementptr(Elementptr elem)
{
    m_elements.push_back(elem);
}

void Domain::addElset(std::string name, std::vector<Elementptr> elset)
{
    m_Elsetmap[name]=elset ;

}

void Domain::addNset(std::string name, std::vector<Nodeptr> nset)
{
    m_Nsetmap[name] = nset;
}

void Domain::setNsetmap(std::map<std::string,std::vector<Nodeptr> > nsetmap)
{
    m_Nsetmap=nsetmap;
}

void Domain::setElsetmap(std::map<std::string,std::vector<Elementptr> > elsetmap )
{m_Elsetmap=elsetmap;}

void Domain::setLoads(std::vector<Load> loads)
{
    m_loads=loads;
}

void Domain::setFixity(std::vector<Fixity> fix)
{
    m_fixity = fix;
}

const std::vector<Elementptr>& Domain::getElements()
{
    return m_elements;
}

const std::vector<Nodeptr>& Domain::getNodes()
{
    return m_nodes;
}



void Oceane::Domain::printNodes()
{
    if(m_nodes.empty()) std::cout<<"EMPTY NODES\n";
    std::cout<<"******NODE DETAILS******\n\n";
    for(auto& iter:m_nodes)
    {
               auto coords = iter->getCoord();
         std::cout<<iter->getIndex()<<"\t"<<coords[0]<<"\t"<<coords[1]<<"\n";
    }
}

void Oceane::Domain::printElements()
{
    if(m_elements.empty()) std::cout<<"EMPTY ELEMENTS\n";
    std::cout<<"*****ELEMENT DETAILS******\n";
    for(auto& iter:m_elements)
    {
        auto nodes =iter->getNodes();

        std::cout<<iter->getIndex()<<"\t";
        for(auto& node:nodes)
        {
            std::cout<<node->getIndex()<<"\t";
        }
        std::cout<<"\n";
    }
}

const Nodeptr Domain::getNodeptr(Index id) const
{
    return m_nodes.at(id-1);
}
const Elementptr Domain::getElementptr(Index id) const
{
    return m_elements.at(id-1);
}

const std::vector<Nodeptr> Domain::getNodes_at(std::vector<Index> ids) const
{
    std::vector<Nodeptr> nodes;
    nodes.reserve(ids.size());
    for(auto iter:ids)
    {
        nodes.push_back(m_nodes.at(iter-1));
    }
    return nodes;
}

const constitutiveFnPtr Domain::getConstitutiveFunction() const {
    return m_constitutiveFn;
}

void Domain::Init()
{


//    Fill the incidents elements data in nodes.
    for(auto& elm:m_elements)
    {
        auto Id = elm->getIndex();
        for(auto& nodes : elm->getNodes())
        {
            nodes->addSharedElementIndex(Id);
            nodes->addSharedElement(elm);
        }
    }

//    collect boundary elements  and create and collect Edge
    for(auto& elm:m_elements)
    {
        if(elm->isBoundary())
        {
            m_boundaryElements.push_back(elm);
            elm->setEdgeptrs();
            auto edges=elm->getEdgeptrs();
//            adding elemenptr in edge
            for(auto&& iter:edges) iter->setElementptr(elm);

            m_edges.insert(m_edges.end(),edges.begin(),edges.end());
        }
    }
//    add the loads in corresponding edges
    for(auto& load:m_loads)
    {
        for(auto& elem:m_Elsetmap[load._elset])
        {
            auto edge=elem->getEdge(load._face);
            edge->setPressure(load._load);
        }
    }

//    add fixity in corresponding nodes
    {
        for(auto& iter:m_fixity)
        {
            auto face =iter._face;
            for(auto& elem:m_Elsetmap[iter._elset])
            {
                elem->getEdge(face)->getNode(0)->setFix();
                elem->getEdge(face)->getNode(1)->setFix();
            }
        }
    }
}

void Domain::Print()
{
    std::cout<<"NODES\n";
    for(auto& node:m_nodes) {
        std::cout<<node->getIndex()<<"\t"<<node->getCoord()[0]<<"\t"
                <<node->getCoord()[1]<<"\n";
    }

    std::cout<<"ELEMENT NODE CONNECTIVITY\n";
    for(auto& element:m_elements) {
        for(auto node:element->getNodes()) {
            std::cout<<node->getIndex()<<"\t";
        }
        std::cout<<"\n";
    }

    std::cout<<"FIXITY DETAILS\n";
    for(auto fix:m_fixity) {
        std::cout<<fix._elset<<"\t"<<fix._face<<"\n";
    }

    std::cout<<"SET OF ELEMENTS\n";
    for(auto set:m_Elsetmap) {
        std::cout<<set.first<<"\t";
        for(auto elem:set.second)
            std::cout<<elem->getIndex()<<"\t";
        std::cout<<"\n";
    }

    std::cout<<"\nLOAD DETAILS\n";
    std::cout<<"ELSET \t FACE \t LOAD \n";
    for(auto load:m_loads) {
        std::cout<<load._elset<<"\t"<<load._face<<"\t"<<load._load<<"\n";
    }
}



} //namespace Oceane
