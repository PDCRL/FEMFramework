#ifndef CONSTRAINT_H
#define CONSTRAINT_H

#include<string>
#include "function.h"
#include "edge.h"
#include "dofmanager.h"
namespace Oceane {
namespace Boundary {

//constrain function is generally parse. so position of variables required
class LinearConstraint: public Oceane::GenericFunction
{
public:
    LinearConstraint(Edgeptr edge);
    ~LinearConstraint()=default;
    std::vector<double> SparseGradientValues(const double*);
    Edgeptr getEdge();
    std::string getName();

protected:

    Edgeptr m_edge;
    std::vector<double> m_gradvals; //gradients stored in order
    double Nx, Ny;
    std::string m_name;
};

}//namespace Boundary
} //namespace Oceane

#endif // CONSTRAINT_H
