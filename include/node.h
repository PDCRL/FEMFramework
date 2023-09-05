#ifndef NODE_H
#define NODE_H

#include<memory>
#include<vector>
#include "fwd.h"
#include <iostream>
#include <map>
namespace Oceane {

enum DoFType
{
    P=0,
    PX,
    PY,
    PXX,
    PYY,
    PXY,
    UX,
    UY,
    RX,
    RY
};


class Edge;
class Element;
using Elementptr = std::shared_ptr<Oceane::Element>;
class Domain;
class Node
{

    friend Edge;
    friend Domain;
public:
    Node();
    Node(Index id, std::vector<double> coord);
    Node(Index id, double  x, double  y);

    std::vector<Index> getSharedElementsIndex() const;
    const std::vector<Elementptr>& getSharedElements() const;
    Index getIndex();

    std::vector<double> getCoord();

    void printNode();
    double getDof(const double* x,DoFType p);
    Index getDofPos(DoFType p);


    std::vector<Index> getDofPos();
    void addDof(DoFType p,Index pos);
    void setFix();
    bool isFixed();

    std::vector<Index> get_StressDof_Indices();
    std::vector<Index> get_DisplacementDof_Indices();


private:


    Index m_id;
    std::vector <double> m_coord;
    std::vector<Elementptr> m_sharedElements;
    std::vector<Index> m_sharedElementsIndex;
    bool m_fix;
    //dof is mapped against its location in array.
    std::map<DoFType,Index> m_dofMap;


//    following functions are only accessed throgh Domain class, Node is declared as friend class
    void addSharedElementIndex(Index id);
    void addSharedElement(Elementptr elem);
//   following function is only accessed through dofmanager class.

};

using Nodeptr = std::shared_ptr<Node>;

} // namespace Oceane

#endif // NODE_H
