#ifndef FEVALUES_H
#define FEVALUES_H
#include <vector>
#include "element.h"
#include "node.h"
#include "quadrature.h"
#include "domain.h"
namespace Oceane    {

class QuadPointVals;
class FEvalues
{
public:
public:
    FEvalues(int quadpts,Elementptr element,Domain* domain);
    void operator()();
    void Init(int quadpoint);
    const std::vector<QuadPointVals>& getQuadPointVals() const;
    const std::vector<Index>& getStressDofIndices();
    const std::vector<Index>& getStrainDofIndices();
    void getElementId();
protected:
    int        m_quadpts;
    Elementptr m_element;
    Domain* m_domain;
    std::vector<QuadPointVals> m_quadPointVals;
    std::vector<Index> m_stressDof;
    std::vector<Index> m_strainDof;
};

class QuadPointVals
{
public:
    QuadPointVals(Oceane::Matrix Bs,Oceane::Matrix Bd,
                  double detJ,double dV);
    const Oceane::Matrix& getBsmat();
    const Oceane::Matrix& getBdmat();
    double get_dV();
    double get_detJ();
private:
    Oceane::Matrix _bsmat;
    Oceane::Matrix _bdmat;
    double _detJ;
    double _dV;
};

} // namespace Oceane

#endif // FEVALUES_H
