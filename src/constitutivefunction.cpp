#include "constitutivefunction.h"

namespace Oceane {

ConstitutiveFunction::ConstitutiveFunction():
    m_name("STEEL"),
    m_poisson(0.3),
    m_modulus(2e5)
{

}

ConstitutiveFunction::~ConstitutiveFunction()
{

}

PlaneStress::PlaneStress()
{

}

PlaneStress::~PlaneStress()
{

}

Oceane::Matrix PlaneStress::getComplianceMatrix()
{
    auto E = m_modulus;
    auto n = m_poisson;
    auto G = E/(2*(1+n));
    Eigen::Matrix3d mat;
    mat<< 1/E,-n/E,0,
            -n/E,1/E,0,
            0,0,1/G;
    return  mat;
}

Oceane::Matrix PlaneStress::getStiffnessMatrix()
{
    auto E =m_modulus;
    auto n= m_poisson;
    Eigen::Matrix3d mat;
    mat<<1.0,n,0.0,
            n,1.0,0.0,
            0.0,0.0,(1-n)/2;
    return E/(1-n*n)*mat;
}

} // namespace Oceane
