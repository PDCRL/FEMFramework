#include "dofmanager.h"
#include <algorithm>
#include <iterator>
namespace Oceane {

DofManager::DofManager(Domain* domain) : m_domain(domain),
    m_x{}
{
   this->Init();
    std::cout<<"\nDofManager Initiated";
    std::cout<<"\nNo. of Variables\t"<<m_Nvar;

}

DofManager::~DofManager()
{
}



void DofManager::update_x(const double *x)
{
    m_x.assign(x,x+m_x.size());
}

Index DofManager::getVariableCount() {

    return m_Nvar;
}

void DofManager::Init() {

    Index count=0;
    for(auto& iter:m_domain->getNodes())
    {
        if(iter->isFixed()) {
            count+=10;
        }
        else {
            count+=8;
        }

    }

    m_Nvar = count;
    m_x.resize(m_Nvar);
    std::fill(m_x.begin(),m_x.end(),0.0);

    size_t index =-1;

    //distributing strain dof
    //distribute displacement dof
    for(auto& iter:m_domain->getNodes()) {
        iter->addDof(UX,++index);
        iter->addDof(UY,++index);
    }
    m_N_dispDof = index+1;
    //distributing stress dof
    for(auto& iter:m_domain->getNodes())
    {
        iter->addDof(P,++index);
        iter->addDof(PX,++index);
        iter->addDof(PY,++index);
        iter->addDof(PXX,++index);
        iter->addDof(PYY,++index);
        iter->addDof(PXY,++index);
    }


   m_N_stressDof =index+1- m_N_dispDof;


    //distribute Reaction dof
    for(auto& iter:m_domain->getNodes()) {
        if(iter->isFixed()) {
            iter->addDof(RX,++index);
            iter->addDof(RY,++index);
        }
    }

    m_N_reactionDof = index+1-m_N_dispDof-m_N_stressDof;

    //Fill the known values in variable list
    for(auto& iter:m_domain->getNodes())
    {
        if(iter->isFixed())
        {
            auto UxIndex = iter->getDofPos(UX);
            m_x[UxIndex] = 0.0;
            auto UyIndex = iter->getDofPos(UY);
            m_x[UyIndex] = 0.0;
        }
    }
    
    //fill m_x_u and m_x_l
    m_x_u.resize(m_Nvar);
    m_x_l.resize(m_Nvar);

//    std::fill(m_x_l.begin(),m_x_l.end(), -1e19);
//    std::fill(m_x_u.begin(),m_x_u.end(),1e19);

    std::fill(m_x_l.begin(),m_x_l.begin()+m_N_dispDof,-1e-3);
    std::fill(m_x_l.begin()+m_N_dispDof,m_x_l.end(),-1000);
    std::fill(m_x_u.begin(),m_x_u.begin()+m_N_dispDof,1e-3);
    std::fill(m_x_u.begin()+m_N_dispDof,m_x_u.end(),1000);

    //fill initial values.
    std::fill(m_x.begin(),m_x.begin()+m_N_dispDof,0);
    std::fill(m_x.begin()+m_N_dispDof,m_x.end(),1);

    //filling scaling
    m_x_scale.resize(m_Nvar);
    std::fill(m_x_scale.begin(),m_x_scale.begin()+m_N_dispDof,1000);
    std::fill(m_x_scale.begin()+m_N_dispDof,m_x_scale.end(),1);
}

void DofManager::Init_x(std::vector<double> inx)
{
    m_x=inx;
}

const std::vector<double> DofManager::get_x()
{
    return m_x;
}

const std::vector<double>DofManager::get_x_l()
{
    return m_x_l;
}

const std::vector<double>DofManager::get_x_scale()
{
    return m_x_scale;
}

const std::vector<double>DofManager::get_x_u()
{
    return m_x_u;
}
void DofManager::Print() {
    std::cout<<"\\n DOF INFO";
    std::cout<<"TOTAL VARIABLES\t:\t"<<m_x.size()<<"\n";
}

Vector DofManager::getDofVector(std::vector<Index> indices)
{
    Vector vec(indices.size());
    for(size_t i=0; i<indices.size(); ++i)
    {
        vec[i]=m_x[indices[i]];
    }
    return  vec;
}

void DofManager::getDofVector(std::vector<Index> indices, Oceane::Vector& vec)
{
    for(size_t i=0; i<indices.size(); ++i)
    {
        vec[i]=m_x[indices[i]];
    }
}


} // namespace Oceane
