#include "optprob.h"
#include <limits>
#include <random>
#include <ctime>
#include <cstdlib>
#include <fstream>
namespace Oceane
{

OptProb::OptProb():m_Nvar(0),m_objfun(nullptr),
m_Nnz_jacobian(0),m_Nnz_hess_lag(0),m_hasHessian(false)
{
}

int OptProb::getNconstraints()
{
    return m_linearConstraints.size();
}



void OptProb::Init(DofManager *dofmanager)
{
    m_Nvar =dofmanager->getVariableCount();
    m_x=dofmanager->get_x();
    m_x_l =dofmanager->get_x_l();
    m_x_u= dofmanager->get_x_u();
    m_x_scale = dofmanager->get_x_scale();
    m_gradient_val.resize(m_Nvar);
    std::fill(m_gradient_val.begin(),m_gradient_val.end(),0.0);
}

void OptProb::Init(std::shared_ptr<DofManager> dofmanager) {

    m_Nvar = dofmanager->getVariableCount();
    m_x.resize(m_Nvar);
    std::fill(m_x.begin(),m_x.end(),0.0);
    m_x_l.resize(m_Nvar);
    std::fill(m_x_l.begin(),m_x_l.end(),-1e19);
    m_x_u.resize(m_Nvar);
    std::fill(m_x_u.begin(),m_x_u.end(),1e19);
}

void OptProb::Init(int n)
{
    m_Nvar =n;
    m_x.resize(m_Nvar);
    std::fill(m_x.begin(),m_x.end(),1.0);
    m_x_l.resize(m_Nvar);
    std::fill(m_x_l.begin(),m_x_l.end(),-1e19);
    m_x_u.resize(m_Nvar);
    std::fill(m_x_u.begin(),m_x_u.end(),1e19);
    m_gradient_val.resize(m_Nvar);
}

void OptProb::Init(Oceane::OptProb *prob)
{
    this->m_Nvar=prob->m_Nvar;
    this->m_x= prob->m_x;
    this->m_x_l=prob->m_x_l;
    this->m_x_u = prob->m_x_u;
    this->m_objfun=prob->m_objfun;
    this->m_objVal = prob->m_objVal;
    this->m_hasHessian=prob->m_hasHessian;
    this->m_linearConstraints = prob->m_linearConstraints;
    this->m_Nnz_jacobian= prob->m_Nnz_jacobian;
    this->m_Nnz_hess_lag = prob->m_Nnz_hess_lag;
}

void OptProb::InitX(std::vector<double> x)
{
    m_x=x;
}

void OptProb::InitX(double val){
    std::fill(m_x.begin(),m_x.end(),val);
}


void OptProb::setRandomInit()
{

      srand(static_cast<unsigned int>(clock()));
      for (size_t i=0; i < m_x.size(); i++) {

         m_x[i]= double(rand()) / (double(RAND_MAX) + 1.0);

      }

}

void OptProb:: setObjfun(ptrFunction fun) {
    m_objfun =fun;
}

void OptProb::addConstraint(ptrFunction fun) {
    m_linearConstraints.emplace_back(fun);
}

void OptProb::setVal(double val)
{
    std::fill(m_x.begin(),m_x.end(),val);
}

int OptProb::getNnzConstraintJacobian()
{
    int number=0;
    for(auto& iter:m_linearConstraints)
    {
        number+=iter->getNNzGradient();
    }
    return number;
}

ptrFunction OptProb::getObjfun()
{
    return m_objfun;
}

std::vector<ptrFunction> OptProb::getLinearConstraints()
{
    return m_linearConstraints;
}

void OptProb::InitX(Index beg, Index end, double val)
{
    for (;beg<end;++beg)
    {
        m_x[beg]=val;
    }
}

std::vector<double>& OptProb::get_lowerbound_ref()
{
    return m_x_l;
}

std::vector<double>& OptProb::get_upperbound_ref()
{
    return m_x_u;
}

void OptProb::print_result(std::string file)
{
    std::ofstream outfile;
    outfile.open(file);
    std::cout<<"OBJECIVE FUNCTION VALUE\t:"<<m_objVal;
    std::cout<<"\nSOLUTION VALUES\n";
    for(auto& iter:m_x)
        outfile<<"\n"<<iter;
}

std::vector<double> OptProb::getX() const
{
    return m_x;
}

void OptProb::get_linearConstraints_matrix(Oceane::Matrix &mat, std::vector<double>& vec)
{
    auto _m = this->getLinearConstraints().size();
    mat.resize(_m,m_Nvar);
    vec.resize(_m);

    for(size_t iRow=0; iRow<m_linearConstraints.size(); ++iRow)
    {
        auto constraint = m_linearConstraints[iRow];
        vec[iRow] = constraint->getVal();

        std::vector<double> ones(m_Nvar,1.0);
        auto vals = constraint->SparseGradientValues(&ones[0]);
        auto pos = constraint->getPositions();

        Eigen::RowVectorXd row(m_Nvar);
        row.setZero();
        for(size_t i=0; i<vals.size(); ++i)
        {
            row(pos[i]) =vals[i];
        }
        mat.row(iRow)=row;
    }
}

}//namespace Oceane

