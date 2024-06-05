#include "../include/node.h"

namespace Oceane {

Node::Node()
{

}

Node::Node(Index id, std::vector<double> coord):m_id(id),m_coord(coord),
    m_sharedElements{},
    m_sharedElementsIndex{},m_fix(false)
{

}

Node::Node(Index id, double  x, double  y):m_id(id),
    m_coord{},
    m_sharedElements{},
    m_sharedElementsIndex{},
    m_fix(false)
  {
      m_coord.push_back(x);m_coord.push_back(y);
  }

std::vector<Index> Node::getSharedElementsIndex() const
{
    return  m_sharedElementsIndex;
}

Index Node::getIndex()
{
    return m_id;
}

void Node::printNode()
{
    std::cout<<m_id;
    for(auto iter:m_coord)
        std::cout<<"\t"<<iter;
    std::cout<<"\n";
}
double Node::getDof(const double* x,DoFType p)
{
    return x[m_dofMap[p]];
}

Index Node::getDofPos(DoFType p)
{
    return m_dofMap[p];
}
std::vector<Index> Node::getDofPos()
{
    std::vector<Index> pos;
    for(auto iter:m_dofMap)
    {
        pos.push_back(iter.second);
    }
    return pos;
}
std::vector<double> Node::getCoord()
{
    return m_coord;
}
void Node::addDof(DoFType p,Index pos)
{
    m_dofMap[p] =pos;
}

void Node::setFix()
{
    m_fix =true;
}

bool Node::isFixed()
{
    return m_fix;
}

void Node::addSharedElementIndex(Index id) {
    m_sharedElementsIndex.push_back(id);
}

void Node::addSharedElement(Elementptr elem){
    m_sharedElements.push_back(elem);
}


std::vector<Index> Node::get_StressDof_Indices() {
    std::vector<Index> stsdof;
    stsdof.push_back(m_dofMap[P]);
    stsdof.push_back(m_dofMap[PX]);
    stsdof.push_back(m_dofMap[PY]);
    stsdof.push_back(m_dofMap[PXX]);
    stsdof.push_back(m_dofMap[PYY]);
    stsdof.push_back(m_dofMap[PXY]);

    return  stsdof;
}

std::vector<Index> Node::get_DisplacementDof_Indices() {
    std::vector<Index> dispdof;
    dispdof.push_back(m_dofMap[UX]);
    dispdof.push_back(m_dofMap[UY]);

    return dispdof;
}

} // namespace Oceane
