#ifndef EDGE_H
#define EDGE_H
#include<memory>
#include "node.h"
#include "function.h"

namespace Oceane {




class Edge
{
    friend Node;

public:
//    constructors
    Edge(Index edgeId, Index elementId,std::vector<Nodeptr> nodes);
    Edge(Index edgeId,Index id,Nodeptr first, Nodeptr second);
    void setElementptr(Elementptr elem);
    void printEdge();
    Nodeptr getNode(Index id);

    Index getIndex(){      return m_id;    }
    void setPressure(double pressure){ m_pressure=pressure; }
    double getPressure()    {        return  m_pressure;    }
    double getNx()  {   return m_nx;    }
    double getNy()  {   return m_ny;    }
    double getTxx() {   return m_pressure*m_nx;    }
    double getTyy() {   return m_pressure*m_ny;    }
    double getLength()  {   return m_length;   }
    bool isFix(){return m_isFix;}
    Index getElementId();


protected:
    Index m_id;
    double m_nx,m_ny;  //directin cosine of normal.
    Index m_elementId;
    Elementptr m_element;
    Nodeptr m_first;
    Nodeptr m_second;
    double m_pressure;
    bool m_isFix;
    double m_length;
private:
//    Following functions are accessed by Domain class only to set properties. Edge is declared as
//    friend in Domain class.

    void  setFix(){ m_isFix=true;}
};


using Edgeptr = std::shared_ptr<Edge>;
} // namespace Oceane

#endif // EDGE_H
