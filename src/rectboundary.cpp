#include "rectboundary.h"

namespace Oceane {

namespace Boundary  {
namespace Rectangle {
namespace Traction  {

Txx::Txx(Edgeptr edge, Index node):LinearConstraint(edge)
{
    auto N = m_edge->getNode(node);
    m_pos.push_back(N->getDofPos(PYY)); // 0
    m_pos.push_back(N->getDofPos(PXY)); // 1
    m_val =m_edge->getTxx();
    m_gradvals.push_back(Nx);
    m_gradvals.push_back(-Ny);
    m_name ="TXX";
}

void Txx:: Eval(const double* x, double& result)
{
    double Pyy=x[m_pos[0]];
    double Pxy =x[m_pos[1]];
    result = Pyy*Nx-Pxy*Ny;
}

Tyy::Tyy(Edgeptr edge,Index node):LinearConstraint(edge)
{
   auto N = m_edge->getNode(node);
   m_pos.push_back(N->getDofPos(PXX)); // 0
   m_pos.push_back(N->getDofPos(PXY)); // 1
   m_val = m_edge->getTyy();
   m_gradvals.push_back(Ny);
   m_gradvals.push_back(-Nx);
   m_name = "TYY";
}

void Tyy::Eval(const double* x,double& result)
{
   auto Pxx = x[m_pos[0]];
   double Pxy = x[m_pos[1]];
   result = Pxx* Ny-Pxy*Nx;
}

Fxx::Fxx(Edgeptr edge):LinearConstraint(edge)
{
   m_pos.push_back(edge-> getNode(0)->getDofPos(PY));
   m_pos.push_back(edge-> getNode(1)->getDofPos(PY));

   if(Nx==0)
   {
       m_gradvals.push_back(Ny);
       m_gradvals.push_back(-Ny);
   }

   if(Ny==0)
   {
       m_gradvals.push_back(-Nx);
       m_gradvals.push_back(Nx);
   }

   auto id= edge->getIndex();
   int m=1;
   if(id==2 or id==3 ) m=-1;

   m_val = m_edge->getTxx() * m_edge->getLength()*m;   //Force along X direction.
   m_name ="FXX";
}

void Fxx:: Eval(const double* x, double& result)
{
    double Py1 = x[m_pos[0]];
    double Py2 = x[m_pos[1]];

    if(Nx==0)
        result = (Py1-Py2)*Ny;
    if(Ny==0)
        result =(Py2-Py1)*Nx;
}

Fyy::Fyy(Edgeptr edge):LinearConstraint(edge)
{
    m_pos.push_back(m_edge-> getNode(0)->getDofPos(PX));
    m_pos.push_back(m_edge-> getNode(1)->getDofPos(PX));

    if(Nx==0)
    {
        m_gradvals.push_back(-Ny);
        m_gradvals.push_back(Ny);
    }

    if(Ny==0)
    {
        m_gradvals.push_back(Nx);
        m_gradvals.push_back(-Nx);
    }

    auto id= edge->getIndex();
    int m=1;
    if(id==2 or id==3 ) m=-1;

    m_val = edge->getTyy()* edge->getLength()*m;
    m_name = "FYY";
}

void Fyy::Eval(const double* x, double& result)
{
    double Px1 = x[m_pos[0]];
    double Px2 = x[m_pos[1]];


    if(Nx==0)
    {
        result = (Px2-Px1)*Ny;
    }
    if(Ny==0)
    {
        result= (Px1-Px2)*Nx;
    }

}

Moment::Moment(Edgeptr edge):LinearConstraint(edge)
{
    x2 =edge-> getNode(1)->getCoord()[0];
    x1 =edge-> getNode(0)->getCoord()[0];
    y1 =edge-> getNode(0)->getCoord()[1];
    y2= edge-> getNode(1)->getCoord()[1];

    m_pos.push_back(edge-> getNode(0)->getDofPos(P));
    m_pos.push_back(edge-> getNode(0)->getDofPos(PX));
    m_pos.push_back(edge->getNode(0)->getDofPos(PY));
    m_pos.push_back(edge-> getNode(1)->getDofPos(P));
    m_pos.push_back(edge-> getNode(1)->getDofPos(PX));
    m_pos.push_back(edge->getNode(1)->getDofPos(PY));

    if(Nx==0.0)
    {
        m_val = (x2*x2-x1*x1)/2 * edge->getTyy();  //Integral of xTy

        m_gradvals.push_back(Ny);
        m_gradvals.push_back(-x1*Ny);
        m_gradvals.push_back(-y1*Ny);
        m_gradvals.push_back(-Ny);
        m_gradvals.push_back(x2*Ny);
        m_gradvals.push_back(y2*Ny);
    }

    if(Ny==0.0)
    {    
        m_val = -(y2*y2-y1*y1)/2 * edge->getTxx();  //Integral of xTy

        m_gradvals.push_back(-Nx);
        m_gradvals.push_back(x1*Nx);
        m_gradvals.push_back(y1*Nx);
        m_gradvals.push_back(Nx);
        m_gradvals.push_back(-x2*Nx);
        m_gradvals.push_back(-y2*Nx);
    }
    m_name = "MOM";
}

void Moment::Eval(double const* x, double& result)
{
    if(Nx ==0.0)
    {
        double P1 =x[m_pos[0]];
        double PX1 = x[m_pos[1]];
        double PY1 = x[m_pos[2]];
        double P2=x[m_pos[3]];
        double PX2=x[m_pos[4]];
        double PY2  = x[m_pos[5]];

        result = P1*Ny - PX1*x1*Ny - PY1*y1*Ny -P2*Ny +PX2*x2*Ny +PY2*y2*Ny;
    }
    if(Ny == 0.0)
    {
        double P1 =x[m_pos[0]];
        double PX1 = x[m_pos[1]];
        double PY1 = x[m_pos[2]];
        double P2=x[m_pos[3]];
        double PX2=x[m_pos[4]];
        double PY2  = x[m_pos[5]];
        result = -P1*Nx +PX1*x1*Nx + PY1*y1*Nx +P2*Nx - PX2*x2*Nx -PY2*y2*Nx;
    }
}

}   //traction namespace

namespace Displacement {

Txx::Txx(Edgeptr edge, Index node):LinearConstraint(edge)
{

    auto N = m_edge->getNode(node);
    m_pos.push_back(N->getDofPos(PYY));
    m_pos.push_back(N->getDofPos(RX));

    m_val =0.0;

    m_gradvals.push_back(Nx);
    m_gradvals.push_back(-1);
}

void Txx::Eval(const double* x, double& result)
{
    double Pyy=x[m_pos[0]];
    double Rxx = x[m_pos[1]];
    result = Pyy*Nx - Rxx;
}

Tyy::Tyy(Edgeptr edge,Index node):LinearConstraint(edge)
{
   auto N = m_edge->getNode(node);
   m_pos.push_back(N->getDofPos(PXX));
   m_pos.push_back(N->getDofPos(RY));
   m_val = 0.0;
   m_gradvals.push_back(Ny);
   m_gradvals.push_back(-1);
}

void Tyy::Eval(const double* x,double& result)
{

   auto Pxx = x[m_pos[0]];
   double Ryy = x[m_pos[1]];
   result = Pxx* Ny - Ryy;
}

Fxx::Fxx(Edgeptr edge):LinearConstraint(edge)
{
   m_pos.push_back(edge-> getNode(0)->getDofPos(PXX));
   m_pos.push_back(edge-> getNode(1)->getDofPos(PXX));
   m_pos.push_back(edge-> getNode(0)->getDofPos(RX));
   m_pos.push_back(edge-> getNode(1)->getDofPos(RY));

   m_gradvals.push_back(-Nx);
   m_gradvals.push_back(Nx);
   m_gradvals.push_back(-(m_edge->getLength())/2);
   m_gradvals.push_back(-(m_edge->getLength())/2);
   m_val = 0.0;   //Force along X direction.
}

void Fxx:: Eval(const double* x, double& result)
{
    double Pyy1 = x[m_pos[0]];
    double Pyy2 = x[m_pos[1]];
    double Rx1 =x[m_pos[2]];
    double Rx2 =x[m_pos[3]];
    result = (Pyy2-Pyy1)*Nx - (Rx2+Rx1)/2*m_edge->getLength();
}

Fyy::Fyy(Edgeptr edge):LinearConstraint(edge)
{
    m_pos.push_back(m_edge-> getNode(0)->getDofPos(PXX));
    m_pos.push_back(m_edge-> getNode(1)->getDofPos(PXX));
    m_pos.push_back(m_edge-> getNode(0)->getDofPos(RY));
    m_pos.push_back(m_edge-> getNode(1)->getDofPos(RY));
    m_gradvals.push_back(Nx);
    m_gradvals.push_back(-Nx);
    m_gradvals.push_back(-(m_edge->getLength())/2);
    m_gradvals.push_back(-(m_edge->getLength())/2);
    m_val = 0.0;
}

void Fyy::Eval(const double* x, double& result)
{
    double Pxx1 = x[m_pos[0]];
    double Pxx2 = x[m_pos[1]];
    double Ry1 = x[m_pos[2]];
    double Ry2 = x[m_pos[3]];
    result = (Pxx2-Pxx1)*Ny - (Ry1+Ry2)/2*m_edge->getLength();
}

Moment::Moment(Edgeptr edge):LinearConstraint(edge)
{
    m_pos.push_back(edge-> getNode(0)->getDofPos(P));
    m_pos.push_back(edge-> getNode(0)->getDofPos(PX));
    m_pos.push_back(edge->getNode(0)->getDofPos(PY));
    m_pos.push_back(edge-> getNode(1)->getDofPos(P));
    m_pos.push_back(edge-> getNode(1)->getDofPos(PX));
    m_pos.push_back(edge->getNode(1)->getDofPos(PY));
    m_pos.push_back(m_edge-> getNode(0)->getDofPos(RY));
    m_pos.push_back(m_edge-> getNode(1)->getDofPos(RY));


    x2 =edge-> getNode(1)->getCoord()[0];
    x1 =edge-> getNode(0)->getCoord()[0];
    y1 =edge-> getNode(0)->getCoord()[1];
    y2= edge-> getNode(1)->getCoord()[1];
    length = m_edge->getLength();
    if(Nx==0.0)
    {

        m_val = (x2*x2-x1*x1)/2 * edge->getTxx();  //Integral of xTy

        m_gradvals.push_back(Ny);
        m_gradvals.push_back(-x1*Ny);
        m_gradvals.push_back(-Ny);
        m_gradvals.push_back(x2*Ny);
    }

    if(Ny==0.0)
    {
        m_pos.push_back(edge-> getNode(0)->getDofPos(P));
        m_pos.push_back(edge-> getNode(0)->getDofPos(PY));
        m_pos.push_back(edge-> getNode(1)->getDofPos(P));
        m_pos.push_back(edge-> getNode(1)->getDofPos(PY));
        m_pos.push_back(m_edge-> getNode(0)->getDofPos(RX));
        m_pos.push_back(m_edge-> getNode(1)->getDofPos(RX));

        m_val = -(y2*y2-y1*y1)/2 * edge->getTxx();  //Integral of xTy

        m_gradvals.push_back(-Nx);
        m_gradvals.push_back(y1*Nx);
        m_gradvals.push_back(Nx);
        m_gradvals.push_back(-y2*Nx);
    }
}

void Moment:: Eval(double const* x, double& result)
{

    if(Nx ==0.0)  {

        double P1 =x[m_pos[0]];
        double PX1 = x[m_pos[1]];
        double P2=x[m_pos[2]];
        double PX2=x[m_pos[3]];

        result = (x2*PX2-x1*PX1+P1-P2)*Ny - (PX1+PX2)/2*length/2;
    }
    if(Ny == 0.0)
    {
        double P1 =x[m_pos[0]];
        double PY1 = x[m_pos[1]];
        double P2=x[m_pos[2]];
        double PY2=x[m_pos[3]];

        result =( P2-P1 + y1*PY1 -y2*PY2)*Ny - (PY1+PY2)/2*length/2;
    }

}

}



}   //Rectangle namespace
}   //Boundary namespace
}   //Oceane namespace
