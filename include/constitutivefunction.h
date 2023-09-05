#ifndef CONSTITUTIVEFUNCTION_H
#define CONSTITUTIVEFUNCTION_H
#include<string>
#include<memory>
#include "matrix.h"
#include <vector>
namespace Oceane {

class ConstitutiveFunction
{
public:
    ConstitutiveFunction();
    virtual ~ConstitutiveFunction();
    virtual Oceane::Matrix getComplianceMatrix(){}
    virtual Oceane::Matrix getStiffnessMatrix(){}
    virtual Oceane::Vector getConstitutiveFunction(const std::vector<double>& stress,
                                                   const std::vector<double>& strain){}
    virtual Oceane::Vector getConstitutiveFunction(const Oceane::Vector& stress,
                                                   const Oceane::Vector& strain){}

protected:
    std::string m_name;
    double m_poisson;
    double m_modulus;
};

class PlaneStress:public ConstitutiveFunction
{
public:
    PlaneStress();
    ~PlaneStress();
    Oceane::Matrix getComplianceMatrix();
    Oceane::Matrix getStiffnessMatrix();



};

using constitutiveFnPtr = std::shared_ptr<ConstitutiveFunction>;

} // namespace Oceane

#endif // CONSTITUTIVEFUNCTION_H
