#include "element.h"

#include <algorithm>
namespace Oceane {

Element::Element(Index id,std::vector<Nodeptr> nodes):m_id(id),m_nodes(nodes),
    m_edges{},m_isBoundary(false),m_init(false) {}

Element::~Element()
{

}

void Element::getBmats_at(double x, double y, Oceane::Matrix &bsmat, Oceane::Matrix &bdmat, double &detJ)
{
    this->getBdmat_at(x,y,bdmat,detJ);
    this->getBsmat_at(x,y,bsmat,detJ);
}

void Element::addNode(Nodeptr node) {
    m_nodes.push_back(node);
}

std::vector<Nodeptr> Element::getNodes()
{
    return  m_nodes;
}

Index Element:: getIndex()
{
    return m_id;
}

void Element::printNodes()
{
    std::cout<<"NODES\n";
    for(auto& node:m_nodes)
    {
        node->printNode();
    }
}

Edgeptr Element::getEdge(Index i)    {
    auto fun=[this](Index i){
        for(auto& iter:this->m_edges) {
            if(iter->getIndex()==i)
                return true;
        }
        return false;
    };
    assert(fun(i)==true);
        for(auto& edge:m_edges)        {
            if(edge->getIndex() == i)
                return  edge;
        }

}

Index Element::elements_shared_by_face(Index i) const
{
    std::vector<Nodeptr> face = getEdgeNodes(i);
    std::vector<Index> v1 = face.at(0)->getSharedElementsIndex() ; std::sort(v1.begin(),v1.end());
    std::vector<Index> v2 = face.at(1)->getSharedElementsIndex() ; std::sort(v2.begin(),v2.end());
    std::vector<Index> intersection ;
    set_intersection(v1.begin(),v1.end(),v2.begin(),v2.end(),std::inserter(intersection,intersection.begin()));
    return intersection.size();
}

std::vector<Nodeptr> Element::getEdgeNodes(Index i) const
{
    auto  endNode = m_nodes.size();
    std::vector<Nodeptr> temp;
    if(i !=endNode-1)
    {
        temp.push_back(m_nodes.at(i));
        temp.push_back(m_nodes.at(i+1));
    }
    else
    {
        temp.push_back(m_nodes.at(endNode-1));
        temp.push_back(m_nodes.at(0));
    }

    return temp;
}


bool Element:: isBoundary(){

    if(m_init==false)
    {
        for(size_t i=0; i<m_nodes.size(); i++)
        {
            if(elements_shared_by_face(i) == 1)
            {
                m_isBoundary= true;
                m_init=true;
                break;
            }
        }
    }

    return m_isBoundary;
}

void Element::setEdgeptrs()
{
//    Assuming no. ef edges are equal to no of nodes
    for(size_t i=0; i<m_nodes.size(); ++i)
    {
        if(this->elements_shared_by_face(i)==1)
        {
            std::vector<Nodeptr> facenodes= this->getEdgeNodes(i);
            Edgeptr edge = std::make_shared<Edge>(i,this->getIndex(),facenodes);
            m_edges.push_back(edge);
        }
    }

}

std::vector<Index> Element::get_StressDof_indices() {
    std::vector<Index> stsdof;
    for(auto& iter:m_nodes) {
        auto nodedof = iter->get_StressDof_Indices();
        stsdof.insert(stsdof.end(),nodedof.begin(),nodedof.end());
    }
    return stsdof;
}

std::vector<Index> Element::get_displacementDof_indices() {
    std::vector<Index> dispdof;
    for(auto& iter:m_nodes) {
        auto nodedof = iter->get_DisplacementDof_Indices();
        dispdof.insert(dispdof.end(),nodedof.begin(),nodedof.end());
    }
    return dispdof;
}











} // namespace Oceane
