#include "objectivefunction.h"
namespace Oceane {
int ObjectiveFunction::FCOUNT =0;
ObjectiveFunction::ObjectiveFunction(Domain* domain,
                                     size_t varCount)
{
    m_Nvar = varCount;
    for(auto&& iter:domain->getElements())
    {
        m_vals.emplace_back(m_quadpoints,iter,domain);
    }
}


ObjectiveFunction::ObjectiveFunction(Domain *domain,
                                     Oceane::DofManager *dofmanager)

{
    m_Nvar = dofmanager->getVariableCount();
    for(auto&& iter:domain->getElements())
    {
        m_vals.emplace_back(m_quadpoints,iter,domain);
    }

}

Oceane::Vector ObjectiveFunction::e
getDofVector(std::vector<Oceane::Index> indices, const double *x)
{
    Oceane::Vector dofs;
    dofs.resize(indices.size());
    for(size_t i=0; i<indices.size(); ++i)
    {
       dofs[i]= x[indices[i]];
    }
    return dofs;
}

void ObjectiveFunction::Eval(const double *x, double &result)
{
    FCOUNT++;
    static int Fcount=0;
    ++Fcount;
    result =0.0;
    for(auto&& iter:m_vals)
    {
        auto stressDof= iter.getStressDofIndices();
        auto strainDof = iter.getStrainDofIndices();
        Oceane::Vector phi = this->getDofVector(stressDof,x);
        Oceane::Vector uvec = this->getDofVector(strainDof,x);
        auto elemvals=iter.getQuadPointVals();
        for(auto&& vals:elemvals)
        {

            Oceane::Matrix Bs = vals.getBsmat();
            Oceane::Matrix Bd = vals.getBdmat();
            double dV = vals.get_dV();
            Oceane::Vector cnfn = Bd*uvec - Bs*phi;
            double cnfnnorm = cnfn.norm();
            double delval = cnfn.norm()*dV;
            result+=delval;

        }
    }
}


bool ObjectiveFunction::Gradient(const double *x, double *grad)
{
    std::vector<double> gradient(m_Nvar,0.0);
    Oceane::Vector phi;
    Oceane::Vector uvec;
    Oceane::Vector stsgrad(24),dispgrad(8);
    for (auto& Elmiter: m_vals)
    {
        auto stressIndices = Elmiter.getStressDofIndices();
        auto strainIndices = Elmiter.getStrainDofIndices();
        phi=getDofVector(stressIndices,x);
        uvec = getDofVector(strainIndices,x);
        stsgrad.setZero();
        dispgrad.setZero();

        auto elemvals=Elmiter.getQuadPointVals();
        for(auto&& Qiter: elemvals)
        {
            Oceane::Matrix Bd = Qiter.getBdmat();
            double  J = Qiter.get_dV();
            Oceane::Matrix Bs  = Qiter.getBsmat();

            Oceane::Vector Fvec = Bd*uvec - Bs*phi;
            double delval = Fvec.norm();
            //setting stress gradient

            Oceane::Vector _stsgrad = (-1*J/delval)*(Fvec.transpose()*Bs);
            stsgrad+=_stsgrad;

            //setting strain gradient
            Oceane::Vector _dispgrad = (1*J/delval)*(Fvec.transpose()*Bd);
            dispgrad += _dispgrad;

        }


        auto szStress =stressIndices.size();

        for(size_t i=0; i<szStress; ++i)
            gradient[stressIndices[i]]+= stsgrad[i];
        auto szStrain = strainIndices.size();
        for(size_t i=0; i<szStrain; ++i)
            gradient[strainIndices[i]]+= dispgrad[i];
    }

    for(size_t i =0; i<gradient.size(); ++i)
    {
        grad[i]=gradient[i];
    }
    return  true;
}


} //namespace Oceane
