#ifndef CPS4_H
#define CPS4_H
#include "fwd.h"
#include "matrix.h"
#include "node.h"
#include "element.h"

namespace Oceane {

template<typename BEngine>
class Quad4 :public Element
{
public:

    virtual ~Quad4(){}
    Quad4(Index id,std::vector<Nodeptr> nodes):Element(id,nodes)
    {
        Oceane::Vector xvec(4);
        Oceane::Vector yvec(4);
        int xi=-1,yi=-1;
        for(auto node:m_nodes)
        {
            auto coord =node->getCoord();
            xvec[++xi]=coord[0];
            yvec[++yi]=coord[1];
        }

        m_B.addCoordinates(xvec,yvec);
    }

    void getBsmat_at(double x, double y, Oceane::Matrix& matrix, double& detJ)
    {
        m_B.getBsmat_at(x,y,matrix,detJ);
    }
    void getBdmat_at(double x, double y, Oceane::Matrix& matrix, double& detJ)
    {
        m_B.getBdmat_at(x,y,matrix,detJ);
    }

private:

    BEngine m_B;   
};







/*...........................................................
 * Global element is typedef-ed as CPS4......................
  ...........................................................*/
class Global
{
public:
    Global()=default;
    void addCoordinates(Oceane::Vector xvec,Oceane::Vector yvec)
    {
        m_xvec =xvec; m_yvec = yvec;
    }
    void getBsmat_at(double x, double y, Oceane::Matrix& bsmat, double& detJ);
    void getBdmat_at(double x, double y, Oceane::Matrix& bdmat, double& detJ);
private:
   Oceane::Vector m_xvec;
   Oceane::Vector m_yvec;
};








/*............................................................
 * ...Parametric class is typedef-ed as CPS4P.................
 * .........................................................*/
class Shape {
public:
    Shape(){
        m_dNx.resize(4);
        m_dNy.resize(4);
        m_JacobMatrix.resize(2,2);
    }

   void addCoordinate(Oceane::Vector xcoord, Oceane::Vector ycoord);
   void update(double x, double y);
   Oceane::Matrix getJacobMatrix();
   Oceane::Vector getdNx();
   Oceane::Vector getdNy();

private:
   Oceane::Vector m_dNx;
   Oceane::Vector m_dNy;
   Oceane::Vector m_xcoord;
   Oceane::Vector m_ycoord;
   Oceane::Matrix m_JacobMatrix;

};
class Parametric
{
public:
    Parametric()=default;
    /* addCoordintes() is called by Quad4 constructor*/
    void addCoordinates(Oceane::Vector xvec,Oceane::Vector yvec);
    void getBsmat_at(double x, double y, Oceane::Matrix& bsmat, double& detJ);
    void getBdmat_at(double x, double y, Oceane::Matrix& bdmat, double& detJ);

private:
    void updateShape(double x, double y);
    Shape m_shape;
};





typedef Quad4<Global> CPS4;
typedef Quad4<Parametric> CPS4P;




} // namespace Oceane

#endif // CPS4_H
