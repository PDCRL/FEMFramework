#ifndef OPTPROB_H
#define OPTPROB_H
#include <memory>
#include <vector>
#include <list>
#include "function.h"
#include "dofmanager.h"
#include "boundary.h"
#include <fstream>
using constraintPtr = std::shared_ptr<Oceane::Boundary::LinearConstraint>;
using ptrFunction = std::shared_ptr<Oceane::GenericFunction>;

namespace Oceane {
class OptProb
{
public:
    OptProb();
    virtual ~OptProb()=default;
    void Init(DofManager* dofmanager);
    void Init(std::shared_ptr<DofManager>);
    void Init(int);
    void Init(Oceane::OptProb* prob);

    std::vector<double> getX() const;
    void InitX(std::vector<double>);
    void InitX(double val);
    void InitX(Index beg, Index end, double val);
    void setRandomInit();
    void setVal(double val);
    int getNconstraints();
    int getNnzConstraintJacobian();

    ptrFunction getObjfun();
    std::vector<ptrFunction> getLinearConstraints();

    void setObjfun(ptrFunction fun);
    void addConstraint(ptrFunction fun);
    void setVariableValues(std::vector<double> initial_x,
                           std::vector<double> x_lower,
                           std::vector<double> x_upper);

    std::vector<double>& get_lowerbound_ref();
    std::vector<double>& get_upperbound_ref();
    void print_result(std::string file);
    void get_linearConstraints_matrix(Oceane::Matrix& mat,std::vector<double>& vec);

protected:
    DofManager* m_dofmgr;
    int m_Nvar;  // No of variables
    ptrFunction m_objfun;
    std::vector<ptrFunction> m_linearConstraints{};
    std::vector<ptrFunction> m_nonlinearConstraints{};
    std::vector<double> m_x;  //initial variables
    std::vector<double> m_x_l; //lower bounds
    std::vector<double> m_x_u;  //upperbounds
    std::vector<double> m_x_scale; //variable scaling.
    double m_objVal{0.0};

    int m_Nnz_jacobian;  // no of nonzeros in jacobians
    int m_Nnz_hess_lag;  // no of nonzeros in hessian
    std::vector<double> m_gradient_val;  //std::vector<double>
    bool m_hasHessian ;
    std::ofstream m_outfile;
};

} //namespace Oceane

#endif // OPTPROB_H
