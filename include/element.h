#ifndef ELEMENT_H
#define ELEMENT_H
#include<iostream>
#include <memory>
#include <vector>
#include "matrix.h"
#include "node.h"
#include "edge.h"
namespace Oceane {
class Element
{
public:

    Element(Index id,std::vector<Nodeptr> nodes);
    virtual ~Element();


//    virtual functions
    virtual void getBsmat_at(double x, double y, Oceane::Matrix& bsmat, double& detj)=0;
    virtual void getBdmat_at(double x, double y, Oceane::Matrix& bdmat, double& detj)=0;
    virtual void getBmats_at(double x, double y, Oceane::Matrix& bdmat,
                             Oceane::Matrix& bsmat, double& detJ);
    //Add node in element
    void addNode(Nodeptr node);

    std::vector<Nodeptr> getNodes();
    Index getIndex();
    bool isBoundary();  
    void printNodes();

    /* return the Edgeptr indicated by the number */
    Edgeptr getEdge(Index i);




    /* Edge is boundary which contain two nodes. if no of nodes are equal to no of edges following
       function will work otherwise implement the virtual function. Domain::Init has to be called
       before invoking this function.
   */
    virtual void setEdgeptrs();

    /* Domain call this function to collect Edge. Element::setEdgeptrs has to be called
    before calling this function*/
    std::vector<Edgeptr> getEdgeptrs(){  return  m_edges;  }

    std::vector<Index> get_StressDof_indices();
    std::vector<Index> get_displacementDof_indices();
protected:
    Index m_id;
    std::vector<Nodeptr> m_nodes;
    std::vector<Edgeptr> m_edges;
    bool m_isBoundary;
    bool m_init;


private:
    Index elements_shared_by_face(Index) const ;
    std::vector<Nodeptr> getEdgeNodes(Index i) const;
};

using Elementptr =std::shared_ptr<Element>;

} // namespace Oceane

#endif // ELEMENT_H
