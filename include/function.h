#ifndef FUNCTIONBASE_H
#define FUNCTIONBASE_H

#include "matrix.h"
#include <iostream>
#include <vector>

using gradient_t =double*;
using hessian_t = Oceane::Matrix;
namespace Oceane {

class GenericFunction
{
public:
    GenericFunction();
    virtual ~GenericFunction();
    virtual void Eval(const double* x,double& result);
    virtual void Eval(const std::vector<double>& vec, double& res);

    /*if function is sparse, return the structure, i.e positions of variables*/
    virtual bool Gradient(const double* x, double* grad);
    virtual bool Gradient(const std::vector<double>& x,std::vector<double>& grad);

    /*To evaluate both function and gradient together*/
    virtual void EvalWithGradient(const double* x, double& result,double* grad);

    virtual int getNNzGradient(); //get number of nonzeros in graadient
    virtual bool SparseGradient(const double*x, int*JCol,double* vals);
    virtual std::vector<double>  SparseGradientValues(const double* );
    virtual std::vector<size_t> SparseGradientStructure();//return the sparsity of gradient.

    void setPositions(std::vector<size_t> pos);
    std::vector<size_t> getPositions();


    virtual bool Hessian(const double* x,int* row,int*col,double* vals);
    virtual int getNNzHessian(); //get number of nonzeros in hessian.
    virtual bool SparseHessianStructure(std::vector<int>& iRow, std::vector<int>& jCol);
    virtual bool SparseHessianValues(double* vals);

    /* if the value of function is known and feeded to m_val, return the value. Used for
    defining constrained functions*/
    void setVal(double);
    double getVal();
    void setU(double);
    void setL(double);
    double getU();  //getUpperlimitoffunction
    double getL();  //getlowerlimit of function
    bool isSparse();
    void setSparse();

protected:
    std::vector<size_t> m_pos{}; //if function is sparse, positions are feeded.
    double m_val{0};
    double m_u{10e19};         //upper limit of function.
    double m_l{10e-19};         //lower limit of function.
    bool m_isSparse{false};

};


} //namespace Oceane.





#endif // FUNCTIONBASE_H
