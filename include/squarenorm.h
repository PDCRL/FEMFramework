#ifndef SQUARENORM_H
#define SQUARENORM_H
#include "function.h"
#include "dofmanager.h"
#include "fevalues.h"
#include "domain.h"
namespace Oceane {
class SquareNorm:public Oceane::GenericFunction
{
public:
    SquareNorm(Domain* domain,Oceane::DofManager* dofmanager);
    ~SquareNorm()=default;
    void Eval(const double* x,double& res);
    bool Gradient(const double* x, double* grad);
private:

    inline Oceane::Vector getDofVector(std::vector<Oceane::Index> indices, const double*x);
    int m_quadpoints;
    Oceane::Index m_Nvar;
    std::vector<Oceane::FEvalues>   m_vals;
};

} //namespace Ocean
#endif // SQUARENORM_H
