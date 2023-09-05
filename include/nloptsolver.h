#ifndef NLOPTSOLVER_H
#define NLOPTSOLVER_H
#include "nlopt.h"
#include "optprob.h"

namespace Oceane {

struct LinearConstraintData{
  int m; //no of functions.
  Oceane::Matrix Kmat;
  std::vector<double> Kvec;
};

class nloptsolver: public OptProb
{


public:

    nloptsolver();
    ~nloptsolver();
    void Solve();
    void print_result(std::string file);
//    static double ObjectiveFunction(unsigned n, const double *x,
//                                  double *grad, void *my_func_data);
//    static double ConstraintFunction(unsigned n, const double* x,
//                                     double *grad, void *constraint);
private:
    nlopt_opt m_opt;
    nlopt_algorithm m_method;
    nlopt_result m_result;
    int m_feval_count{0};
    LinearConstraintData* m_data;
};

} //namespace Oceane

#endif // NLOPTSOLVER_H
