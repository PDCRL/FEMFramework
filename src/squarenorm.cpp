#include "squarenorm.h"
namespace Oceane {
SquareNorm::SquareNorm(Domain *domain, Oceane::DofManager *dofmanager)
    :m_quadpoints(10)
{
    m_Nvar = dofmanager->getVariableCount();
    for(auto&& iter:domain->getElements())
    {
        m_vals.emplace_back(m_quadpoints,iter,domain);
    }
}

Oceane::Vector SquareNorm::
getDofVector(std::vector<Oceane::Index> indices, const double *x)
{
    Oceane::Vector dofs;
    dofs.resize(indices.size());
    for(size_t i=0; i<indices.size(); ++i) {
       dofs[i]= x[indices[i]];
    }
    return dofs;
}

void SquareNorm::Eval(const double *x, double &result)
{
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
            double delval = cnfn.squaredNorm()*dV;
            result+=delval;
        }
    }
}


bool SquareNorm::Gradient(const double *x, double *grad)
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

            Oceane::Vector stressVec =Bs*phi;
            Oceane::Vector strainVec = Bd*uvec;
            Oceane::Vector Fvec = strainVec - stressVec;
//            double delval = Fvec.squaredNorm();
            //setting stress gradient

            Oceane::Vector _stsgrad = (-2*Bs.transpose()*Fvec)*J;
            stsgrad+=_stsgrad;

            //setting strain gradient
            Oceane::Vector _dispgrad = (2*Bd.transpose()*Fvec)*J;
            dispgrad += _dispgrad;
        }


        auto szStress =stressIndices.size();

        for(size_t i=0; i<szStress; ++i)
            gradient[stressIndices[i]]+= stsgrad[i];
        auto szStrain = strainIndices.size();
        for(size_t i=0; i<szStrain; ++i)
            gradient[strainIndices[i]]+= dispgrad[i];
    }

    for(size_t i=0; i<m_Nvar; ++i)
        grad[i]=gradient[i];
    return true;
}


} //namespace Oceane
