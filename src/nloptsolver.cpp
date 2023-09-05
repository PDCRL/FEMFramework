#include "nloptsolver.h"
#include <algorithm>
namespace Oceane {
static int nlopt_fcount = 0;
static double ObjectiveFunction(unsigned n, const double *x,
                                    double *grad, void *my_func_data)
{
    ++nlopt_fcount;

//    Oceane::GenericFunction* objfun= (Oceane::GenericFunction*) my_func_data;
    Oceane::GenericFunction* objfun= static_cast<Oceane::GenericFunction*>(my_func_data);
    double res=0;
    objfun->Eval(x,res);
    if(grad)
        objfun->Gradient(x,grad);
    std::cout<<"\nFval("<<nlopt_fcount<<")\t"<<res;
    return res;
}


double ConstraintFunction(unsigned n, const double*x,
                                       double *grad, void *constraint )
{

//    std::cout<<"\n Constraint Evaluated\n";
    Oceane::GenericFunction* function = static_cast<Oceane::GenericFunction*>(constraint);
    double res;
    function->Eval(x,res);
//    std::cout<<res<<"\n";

    if(grad) {
        if(function->isSparse()) {
            auto pos = function->getPositions();
            auto vals=function->SparseGradientValues(x);

            for(size_t i=0; i<n; ++i)
            {
                grad[i]=0;
            }

            for(size_t i=0;i<pos.size(); ++i)
            {
                grad[pos[i]] =vals[i];
            }
        }
        else {
            function->Gradient(x,grad);
        }

    }
    return  res-function->getVal();  // f(x) -c =0;
}


void LinearConstraintFunctions(unsigned m, double* result, unsigned n,
                                 const double* x, double* grad, void * fdata )
{
    Oceane::LinearConstraintData* data = static_cast<Oceane::LinearConstraintData*>(fdata);
    //gradient data
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> A;
    A = data->Kmat;
    std::copy(A.data(),A.data()+m*n,grad);

    Eigen::Map<const Eigen::VectorXd> B(&data->Kvec[0],m); // mapping Kvec into B
    Eigen::Map<const Eigen::VectorXd> X(x,n); //mapping x into X
    Eigen::Map<Eigen::VectorXd> R(result,m); //mapping result into R.

    R= A*X -B;   //constraints evaluations.

}

double NormFunction(unsigned n, const double* x, double* grad, void* fdata)
{
    double value=0.0;
    std::cout<<"\n Constraint value\t:\t" ;
    Oceane::LinearConstraintData* data = static_cast<Oceane::LinearConstraintData*>(fdata);
    Eigen::MatrixXd A = data->Kmat;
    auto m = data->Kvec.size();
    Eigen::Map<const Eigen::VectorXd> B(&data->Kvec[0],m); // mapping Kvec into B
    Eigen::Map<const Eigen::VectorXd> X(x,n); //mapping x into X
    Eigen::VectorXd R = A*X-B;
    value= R.squaredNorm();
    std::cout<<value;
    if(grad) {

    }

}

nloptsolver::nloptsolver()
{

}

nloptsolver::~nloptsolver()
{
    delete m_data;
}


void nloptsolver::Solve()
{
    //Initializing m_data structure
    m_data =new LinearConstraintData ();
    this->get_linearConstraints_matrix(m_data->Kmat,m_data->Kvec);

    m_method=NLOPT_LD_SLSQP;
    double toler = 1e-12;


    auto Nconstr = m_linearConstraints.size();
    double stopval = 1e-8;
    std::vector<double> toleranceVec(Nconstr,1e-8);
    nlopt_opt opt = nlopt_create(m_method,m_Nvar);
    nlopt_set_lower_bounds(opt,&m_x_l[0]);
    nlopt_set_upper_bounds(opt,&m_x_u[0]);

    nlopt_set_min_objective(opt,ObjectiveFunction,m_objfun.get());
    nlopt_add_equality_mconstraint(opt,Nconstr,LinearConstraintFunctions,
                                   m_data,&toleranceVec[0]);
//    nlopt_add_equality_constraint(opt,Nconstr,ConstraintFunction, m_data, toler);

//    nlopt_set_stopval(opt,stopval);
//    double toler = 1e-12;
    nlopt_set_xtol_rel(opt, toler);
    nlopt_set_ftol_rel(opt, toler);
    nlopt_set_ftol_abs(opt, toler);

    std::vector<double> xi = m_x;
    double min_f;
    m_result = nlopt_optimize(opt, &m_x[0], &min_f);
    m_feval_count= nlopt_get_numevals(opt);
    std::cout<<"\n OUTPUT\t"<<m_result;
    if(m_result>0) {
            m_objVal = min_f;
            std::cout<<"\nTotal Function Count"<<nlopt_fcount;
            std::cout<<"\nFEval count\t"<<m_feval_count;
            std::cout<<"\nObjective Function\t"<<m_objVal<<"\n";
            for(auto iter:m_x)
            {
                std::cout<<"\n"<<iter;
            }
    }
    else
        std::cout<<"\nProblem Failed\n";


    nlopt_destroy(opt);
}

void nloptsolver::print_result(std::string file)
{
    if(m_result<0)
        return;

    std::ofstream outfile;
    outfile.open(file);
    std::cout<<"\nALGORITHM\t:\t"<<m_method;
    std::cout<<"\nRESULT\t:\t"<<m_result;
    std::cout<<"\nITERATIONS\t:\t"<<m_feval_count;
    std::cout<<"\nOBJECIVE FUNCTION VALUE\t:\t"<<m_objVal;
    std::cout<<"\nSOLUTION VALUES\n";
    for(auto& iter:m_x)
        outfile<<"\n"<<iter;
}



}//namespace Oceane
