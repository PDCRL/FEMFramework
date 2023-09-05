#ifndef OBJECTIVEFUNCTION_H
#define OBJECTIVEFUNCTION_H
#include "function.h"
#include "dofmanager.h"
#include "fevalues.h"
#include "domain.h"
namespace Oceane {
class ObjectiveFunction:public Oceane::GenericFunction
{
public:

    ObjectiveFunction(Domain* domain, Index varCount);
    ObjectiveFunction(Domain* domain,Oceane::DofManager* dofmanager);
    ~ObjectiveFunction()=default;
    void Eval(const double* x,double& res);
    bool Gradient(const double* x, double* grad);
    static int FCOUNT;
private:
    inline Oceane::Vector getDofVector(std::vector<Oceane::Index> indices, const double*x);
    int m_quadpoints{10};
    Oceane::Index m_Nvar;
    std::vector<Oceane::FEvalues>   m_vals;
};

} //namespace Ocean
#endif // OBJECTIVEFUNCTION_H
