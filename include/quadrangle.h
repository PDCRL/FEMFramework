#ifndef QUADRANGLE_H
#define QUADRANGLE_H
#include "boundary.h"
namespace Oceane {
namespace Boundary {

namespace Quad {
/* Boundary condition Txx. constructor require Edge and node no. indicating respective node*/
class Txx: public LinearConstraint
{
public:
    Txx(Edgeptr edge, Index node);
    ~Txx()=default;

    void Eval(const double* x, double& result);
    std::vector<double> SparseGradient(const double*)
    {
        std::vector<double> gradvals;
        gradvals.push_back(m_edge->getNx());
        gradvals.push_back(-(m_edge->getNy()));
        return gradvals;
    }
};


/* Boundary condition Tyy. constructor require Edge and node no. indicating respective node*/
class Tyy:public LinearConstraint
{
public:
    Tyy(Edgeptr edge,Index node);
    ~Tyy()=default;
    void Eval(const double* x,double& result)
    {

       auto Pxx = x[m_pos[0]];
       double Pxy = x[m_pos[1]];
       result = Pxx* m_edge->getNy() - Pxy*m_edge->getNx();
    }

    std::vector<double> SparseGradient(const double*)
    {
        std::vector<double> gradvals;
        gradvals.push_back(m_edge->getNy());
        gradvals.push_back(-(m_edge->getNx()));
        return gradvals;
    }
};


class FEVal;
using FEValPtr = std::shared_ptr<FEVal>;
/* constructor requires Edgeptr, FEVal  */
class Fxx: public LinearConstraint
{
public:
    Fxx(Edgeptr edge,FEValPtr vals):LinearConstraint (edge){

        if(m_edge->getNx() == 1.0 )
        {
            auto p1 = m_edge-> getNode(0)->getDofPos(PYY);
            auto p2 = m_edge-> getNode(1)->getDofPos(PYY);
            m_pos.push_back(p1);
            m_pos.push_back(p2);
        }

        else if(m_edge->getNy() == 1.0)
        {
            auto p1 = m_edge-> getNode(0)->getDofPos(PXY);
            auto p2 = m_edge-> getNode(1)->getDofPos(PXY);
            m_pos.push_back(p1);
            m_pos.push_back(p2);
        }

        else
        {
            auto p1 =m_edge-> getNode(0)->getDofPos();
            m_pos.insert(m_pos.end(),p1.begin(),p1.end());
            auto p2 =m_edge-> getNode(1)->getDofPos();
            m_pos.insert(m_pos.end(),p2.begin(),p2.end());
        }

        m_fe = vals;
    }

    Fxx(Edgeptr edge):LinearConstraint (edge){

        if(m_edge->getNx() == 1.0 )
        {
            auto p1 = m_edge-> getNode(0)->getDofPos(PYY);
            auto p2 = m_edge-> getNode(1)->getDofPos(PYY);
            m_pos.push_back(p1);
            m_pos.push_back(p2);
        }

        else if(m_edge->getNy() == 1.0)
        {
            auto p1 = m_edge-> getNode(0)->getDofPos(PXY);
            auto p2 = m_edge-> getNode(1)->getDofPos(PXY);
            m_pos.push_back(p1);
            m_pos.push_back(p2);
        }

        else
        {
            auto p1 =m_edge-> getNode(0)->getDofPos();
            m_pos.insert(m_pos.end(),p1.begin(),p1.end());
            auto p2 =m_edge-> getNode(1)->getDofPos();
            m_pos.insert(m_pos.end(),p2.begin(),p2.end());
        }
    }


    ~Fxx()=default;
    void Eval(const double*, double& result);
    std::vector<double> SparseGradient(const double*);
private:
    FEValPtr m_fe;
};



class Fyy : public LinearConstraint
{
public:
    Fyy(Edgeptr edge):LinearConstraint(edge){}
    ~Fyy(){}
    void impl_eval(const double* x,double& result);
    virtual void impl_grad(const double*,gradient_t);
};



class Moment :public Oceane::GenericFunction
{
public:
    Moment(){}
    ~Moment(){}
    void impl_eval(const double* x,double& result);
    virtual void impl_grad(const double*,gradient_t);
private:
};

}//namespace Quad
} //namespace Boundary
}//namespace Oceane



#endif // QUADRANGLE_H
