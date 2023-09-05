#include "nloptsolverplus.h"
#include "nloptsolver.h"
#include <algorithm>
namespace Oceane {
static int nlopt_fcount = 0;
double nloptsolverplus::ObjectiveFunction(unsigned n, const double *x,
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

double nloptsolverplus::ConstraintFunction(unsigned int n, const double *x, double *grad, void *fdata)
{
    double res =0.0;
    Oceane::LinearConstraintData* data = static_cast<Oceane::LinearConstraintData*>(fdata);
    std::cout<<"\nConstraint Value\t:\t";
    //gradient data
    int m=data->Kvec.size();
    Eigen::MatrixXd A;
    A = data->Kmat;
    Eigen::Map<const Eigen::VectorXd> B(&data->Kvec[0],m); // mapping Kvec into B
    Eigen::Map<const Eigen::VectorXd> X(x,n); //mapping x into X
    Eigen::VectorXd R(m);

    R= A*X-B;   //constraints evaluations.
    res = R.squaredNorm();
    std::cout<<res;
    if(grad){
        Eigen::Map<Eigen::VectorXd> gradient (grad,n);
        gradient = 2*A.transpose()*R;
    }

    return res;
}


void nloptsolverplus::LinearConstraintFunctions(unsigned m, double* result, unsigned n,
                                 const double* x, double* grad, void * fdata )
{
    Oceane::LinearConstraintData* data = static_cast<Oceane::LinearConstraintData*>(fdata);

    //gradient data
    Eigen::Matrix<double,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> A;
    A = data->Kmat;
    if(grad) {
    std::copy(A.data(),A.data()+m*n,grad);
    }
    Eigen::Map<const Eigen::VectorXd> B(&data->Kvec[0],m); // mapping Kvec into B
    Eigen::Map<const Eigen::VectorXd> X(x,n); //mapping x into X
    Eigen::Map<Eigen::VectorXd> R(result,m); //mapping result into R.

    R= A*X-B;   //constraints evaluations.
}

nloptsolverplus::nloptsolverplus()
{

}

nloptsolverplus::~nloptsolverplus()
{
    delete m_data;
}


void nloptsolverplus::Solve()
{
    //Initializing m_data structure
    m_data =new LinearConstraintData ();
    this->get_linearConstraints_matrix(m_data->Kmat,m_data->Kvec);

    m_method=nlopt::LD_SLSQP;
    nlopt::opt opt(m_method,m_Nvar);
    opt.set_lower_bounds(m_x_l);
    opt.set_upper_bounds(m_x_u);

    //adding objective function
    opt.set_min_objective(ObjectiveFunction, m_objfun.get());

    //adding vector constraints.
    auto Nconstr = m_linearConstraints.size();
//    std::vector<double> toleranceVec(Nconstr,1e-12);
//    opt.add_equality_mconstraint(LinearConstraintFunctions,m_data,toleranceVec);
    //adding nonlinear constraint
    opt.add_equality_constraint(ConstraintFunction,m_data,1e-12);
    //termination
    double stopval = 1e-15;
    double toler = 1e-10;
//    opt.set_stopval(stopval);
//    opt.set_ftol_abs(toler);
//    opt.set_ftol_rel(toler);
//    opt.set_xtol_abs(toler);
    opt.set_xtol_rel(toler);

    //strat optimization
    double min_f=0.0;
    m_result = opt.optimize(m_x,min_f);

    m_feval_count= opt.get_numevals();
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

}

void nloptsolverplus::print_result(std::string file)
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
