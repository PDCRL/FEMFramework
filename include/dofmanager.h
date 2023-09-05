#ifndef DOFMANAGER_H
#define DOFMANAGER_H
#include <vector>
#include <list>
#include "domain.h"
#include "fwd.h"
#include "matrix.h"
namespace Oceane {
class DofManager
{
public:

    DofManager(Domain* domain);
    ~DofManager();

    //coordinate with nodes and fill dof information.
    void Init();
    void Init_x(std::vector<double>);
    void update_x(const double*);
    void update_x(std::vector<double>);
    Index getVariableCount();
    const std::vector<double> get_x();
    const std::vector<double> get_x_u();
    const std::vector<double> get_x_l();
    const std::vector<double> get_x_scale();
    int getN_objvars();
    void Print();
    Vector getDofVector(std::vector<Index> indice);
    void getDofVector(std::vector<Index>,Oceane::Vector&);


private:
    Domain* m_domain;
    Index m_Nvar;
    std::vector<double> m_x;
    std::vector<double> m_x_u;
    std::vector<double> m_x_l;
    std::vector<double> m_x_scale;
    int m_N_objvars;
    int m_N_stressDof;
    int m_N_dispDof;
    int m_N_reactionDof;
};

} // namespace Oceane

#endif // DOFMANAGER_H
