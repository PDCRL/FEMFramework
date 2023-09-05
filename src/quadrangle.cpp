#include "quadrangle.h"

namespace Oceane {
namespace Boundary {
namespace Quad {

Txx::    Txx(Edgeptr edge, Index node):LinearConstraint(edge)
{

    m_val =m_edge->getTxx();
    auto N = m_edge->getNode(node);
    m_pos.push_back(N->getDofPos(PYY));
    m_pos.push_back(N->getDofPos(PXY));
}



void Txx:: Eval(const double* x, double& result)
{
    double Nx = m_edge->getNx();
    double Ny = m_edge->getNy();
    double Pyy=x[m_pos[0]];
    double Pxy= -x[m_pos[1]];
    result = Pyy*Nx - Pxy*Ny;
}

Tyy::Tyy(Edgeptr edge,Index node):LinearConstraint(edge)
{
   m_val = m_edge->getPressure();
   auto N = m_edge->getNode(node);
   m_pos.push_back(N->getDofPos(PXX));
   m_pos.push_back(N->getDofPos(PXY));
}

}   //namespace Quad

}   //namespace Bounadry

}   //namespace Oceane

