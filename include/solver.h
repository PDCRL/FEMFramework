#ifndef SOLVER_H
#define SOLVER_H

#include "optprob.h"

namespace Oceane {

template <typename Optimizer>
class Solver : public OptProb
{
public:
    Solver() {
        m_optimizer =static_cast<Optimizer*>(this);
             }
    ~Solver() { }

    void Solve(){m_optimizer->Solve();}
    void setOptions(){m_optimizer->setOptions();}
protected:
    Optimizer* m_optimizer;
};

} // namespace Oceane

#endif // SOLVER_H
