#ifndef RECTANGLE_H
#define RECTANGLE_H
#include  "boundary.h"
namespace Oceane {
namespace Boundary
{
namespace Rectangle {

namespace Traction {

/* Boundary condition Txx. constructor require Edge and node no. indicating respective node*
    in constructor
        m_gradvals should be assigned.
*/
class Txx:public Boundary::LinearConstraint
{
public:
    Txx(Edgeptr edge, Index node);
    ~Txx()=default;
    void Eval(const double* x, double& result);
};


/* Boundary condition Tyy. constructor require Edge and node no. indicating respective node*/
class Tyy:public LinearConstraint
{
public:
    Tyy(Edgeptr edge,Index node);
    ~Tyy()=default;
    void Eval(const double* x,double& result);
};


class Fxx: public LinearConstraint
{
public:
    Fxx(Edgeptr edge);
    ~Fxx()=default;
    void Eval(const double* x, double& result);
};

class Fyy:public LinearConstraint
{
public:
    Fyy(Edgeptr edge);
    ~Fyy()=default;
    void Eval(const double* x, double& result);
};

class Moment:public LinearConstraint
{
public:
    Moment(Edgeptr edge);
    ~Moment()=default;
    void Eval(double const* x, double& result);
private:
    double x1,x2,y1,y2;
};

} //namespace Traction

namespace Displacement {

class Txx:public Boundary::LinearConstraint
{
public:
    Txx(Edgeptr edge, Index node);
    void Eval(const double* x, double& result);
    ~Txx()=default;
};


/* Boundary condition Tyy. constructor require Edge and node no. indicating respective node*/
class Tyy:public LinearConstraint
{
public:
    Tyy(Edgeptr edge,Index node);
    ~Tyy()=default;
    void Eval(const double* x,double& result);
};


class Fxx: public LinearConstraint
{
public:
    Fxx(Edgeptr edge);
    void Eval(const double* x, double& result);
    ~Fxx()=default;
};

class Fyy:public LinearConstraint
{
public:
    Fyy(Edgeptr edge);
    void Eval(const double* x, double& result);
    ~Fyy()=default;
};

class Moment:public LinearConstraint
{
public:
    Moment(Edgeptr edge);
    void Eval(double const* x, double& result);
    ~Moment()=default;
private:
    double x1,x2,y1,y2,length;
};



}//namespace Displacement

} //namespace Rectangle
} //namespace Boundary

} //namespace Oceane

#endif // RECTANGLE_H
