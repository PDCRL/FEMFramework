#include "fevalues.h"

namespace Oceane {

QuadPointVals::QuadPointVals(Oceane::Matrix Bs,Oceane::Matrix Bd,
              double detJ,double dV)
    :_bsmat(Bs),_bdmat(Bd),_detJ(detJ),_dV(dV){}

double QuadPointVals::get_dV() {
    return _dV;
}

double QuadPointVals::get_detJ() {
    return  _detJ;
}

const Oceane::Matrix&  QuadPointVals::getBdmat()
{
    return  _bdmat;
}

const Oceane::Matrix&  QuadPointVals::getBsmat()
{
    return  _bsmat;
}

FEvalues::FEvalues(int quadpts,Elementptr elem,Domain* domain)
    :m_quadpts(quadpts),m_element(elem),m_domain(domain)
{

    Init(quadpts);
    m_stressDof=m_element->get_StressDof_indices();
    m_strainDof=m_element->get_displacementDof_indices();
}

void FEvalues::Init(int quadpts)
{
    Oceane::Quadrature quad(quadpts);
    auto nodes = quad.get_nodes();
    auto weight = quad.get_weight();

    for(size_t i=0; i<nodes.size(); ++i) {
        for(size_t j=0; j<nodes.size(); ++j) {
            double x,y;
            x=nodes[i];
            y=nodes[j];

            Oceane::Matrix Bs(3,24),Bd(3,8);
            double detJ=0; double dV=0.0;

            m_element->getBsmat_at(x,y,Bs,detJ);
            m_element->getBdmat_at(x,y,Bd,detJ);
            Eigen::IOFormat OctaveFmt(Eigen::StreamPrecision, 0, ", ", ";\n", "", "", "[", "]");
            dV= weight[i]*weight[j]*detJ;
//            std::cout<<"\nBsmat\n"<<Bs.format(OctaveFmt)<<"\nBdmat\n"<<Bd.format(OctaveFmt);
//            std::cout<<"\n dV\n"<<dV;

            Oceane::Matrix C = m_domain->getConstitutiveFunction()->getComplianceMatrix();
            assert(C.cols()==Bs.rows());
            Bs = C*Bs;
            m_quadPointVals.emplace_back(Bs,Bd,detJ,dV);
        }
    }
}

const std::vector<Index>& FEvalues::getStrainDofIndices()
{
    return m_strainDof;
}

const std::vector<Index>& FEvalues::getStressDofIndices()
{
    return m_stressDof;
}


void FEvalues::operator()()
{


}

const std::vector<QuadPointVals>& FEvalues::getQuadPointVals() const
{
    return m_quadPointVals;
}

void FEvalues::getElementId()
{
    m_element->getIndex();
}

} // namespace Oceane
