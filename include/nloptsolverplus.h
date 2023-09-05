#ifndef NLOPTSOLVERPLUS_H
#define NLOPTSOLVERPLUS_H

#include "nlopt.hpp"
#include "optprob.h"

namespace Oceane {
struct LinearConstraintData;

class nloptsolverplus: public OptProb
{


public:

    nloptsolverplus();
    ~nloptsolverplus();
    void Solve();
    void print_result(std::string file);
    static double ObjectiveFunction(unsigned n, const double *x,
                                  double *grad, void *my_func_data);
    static void LinearConstraintFunctions(unsigned m, double* result, unsigned n,
                                            const double* x, double* grad, void * fdata);
    static double ConstraintFunction(unsigned n, const double*x,
                                           double *grad, void *constraint );
    private:

    nlopt::algorithm m_method;
    nlopt::result m_result;
    int m_feval_count{0};
    LinearConstraintData* m_data;
};

} //namespace Oceane


#endif // NLOPTSOLVERPLUS_H
