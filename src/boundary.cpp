#include "boundary.h"
namespace  Oceane {

namespace Boundary {

std::string LinearConstraint::getName()
{
    return m_name;
}

Edgeptr LinearConstraint::getEdge()
{
    return m_edge;
}

LinearConstraint::LinearConstraint(Edgeptr edge):
    Oceane::GenericFunction ()
{
    m_edge=edge;
    Nx = edge->getNx();
    Ny = edge->getNy();
    m_isSparse=true;
}


std::vector<double> LinearConstraint::SparseGradientValues(const double *)
{
    return m_gradvals;
}


} // namespace Boundary


} // namespace Oceane
