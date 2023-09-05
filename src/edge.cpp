#include "edge.h"

namespace Oceane {

Edge::Edge(Index face, Index elementId,std::vector<Nodeptr> nodes)
    :m_id(face),m_elementId(elementId),m_pressure(0.0),
      m_isFix(false)
{
    m_first = nodes[0];
    m_second=nodes[1];

    double x1 = m_first->getCoord()[0];
    double y1 = m_first->getCoord()[1];

    double x2 = m_second->getCoord()[0];
    double y2 = m_second->getCoord()[1];


    //Direction cosine of the line of face
    double Nx= (x2-x1)/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    double Ny=(y2-y1)/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

    //Direction cosine of the normal
    m_nx = Ny;  m_ny=-Nx;

    m_length = sqrt((x1-x2)*(x1-x2)+(y2-y1)*(y2-y1));

}


Edge::Edge(Index edgeId, Index id,Nodeptr first, Nodeptr second)
    :m_id(edgeId),m_elementId(id),
    m_first(first),m_second(second),m_pressure(0.0),m_isFix(false)
{

    double x1 = m_first->getCoord()[0];
    double y1 = m_first->getCoord()[1];

    double x2 = m_second->getCoord()[0];
    double y2 = m_second->getCoord()[1];


    //Direction cosine of the line of face
    double Nx= (x2-x1)/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));
    double Ny=(y2-y1)/sqrt((x2-x1)*(x2-x1)+(y2-y1)*(y2-y1));

    //Direction cosine of the normal
    m_nx = Ny;  m_ny=-Nx;

    m_length = sqrt((x1-x2)*(x1-x2)+(y2-y1)*(y2-y1));


}

void Edge::printEdge()
{
       std::cout<<"\n Element ID"<<"\t"<<m_elementId<<"\n";
       std::cout<<"Nodes\t"<<m_first->getIndex()<<"\t"<<m_second->getIndex()<<"\n";
}

Index Edge::getElementId()
{
    return m_elementId;
}

void Edge::setElementptr(Elementptr elem)
{
    m_element=elem;
}

Nodeptr Edge::getNode(Index id)
{
    assert(id==0 or id==1);
    if(id==0) return m_first;
    else if (id==1) {
       return m_second;
    }
    else
    {
       std::cerr<<"Invalid Node number";
    }
}

} // namespace Oceane
