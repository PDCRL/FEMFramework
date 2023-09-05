#include "function.h"

namespace Oceane {

GenericFunction::GenericFunction()
{

}

GenericFunction::~GenericFunction()
{

}

void GenericFunction::Eval(const double *x, double &result)
{
    std::cerr<<"No Function is implemented";
}

void GenericFunction::Eval(const std::vector<double> &vec, double &res)
{
    this->Eval(&vec[0],res);
}

double GenericFunction::getU()
{
    return m_u;
}
double GenericFunction::getL()
{
    return m_l;
}


bool GenericFunction::Gradient(const double *, double *)
{
    std::cerr<<"No Gradient Function Evaluation Definition is implemented";
    return false;
}

bool GenericFunction::Gradient(const std::vector<double> &x, std::vector<double>& grad)
{
    this->Gradient(&x[0],&grad[0]);
}

int GenericFunction::getNNzGradient()
{
    return m_pos.size();
}
bool GenericFunction::SparseGradient(const double *x, int *jCol, double *val)
{
    return false;
}

std::vector<double> GenericFunction::SparseGradientValues(const double*x)
{
    std::cerr<<"Don't have any values";

}

std::vector<size_t> GenericFunction::SparseGradientStructure()
{
    std::cerr<<"No structure defined";
}

bool GenericFunction::Hessian(const double *, int *row, int *col, double *val)
{
    return false;
}

int GenericFunction::getNNzHessian()
{
    return false;
}

bool GenericFunction::SparseHessianStructure(std::vector<int> &iRow, std::vector<int> &jCol)
{
    return false;
}

bool GenericFunction::SparseHessianValues(double *vals)
{
    return false;
}

bool GenericFunction::isSparse()
{
    return m_isSparse;
}

void GenericFunction::setSparse()
{
    m_isSparse =true;
}

std::vector<size_t> GenericFunction::getPositions()
{
   return m_pos;
}

void GenericFunction::setPositions(std::vector<size_t> pos)
{
    m_pos = pos;
}

double GenericFunction::getVal()
{
    return m_val;
}

void GenericFunction::setL(double val)
{
    m_l = val;
}


void GenericFunction::setU(double val)
{
    m_u=val;
}

void GenericFunction::setVal(double val)
{
    m_val =val;
}

void GenericFunction::EvalWithGradient(const double *x, double &result, double *grad)
{
    this->Eval(x,result);
    this->Gradient(x,grad);
}
} //namespace Oceane.


