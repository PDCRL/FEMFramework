#ifndef RECT4_H
#define RECT4_H
#include"matrix.h"
#include "element.h"

namespace Oceane {

class Rect4:public Element
{
public:
    Rect4(Index id,std::vector<Nodeptr> nodes);
    ~Rect4()=default;
    void getBsmat_at(double x, double y, Matrix& bsmat, double& detJ);
    void getBdmat_at(double x, double y, Matrix& bdmat, double& detJ);
    void getBmats_at(double x, double y, Matrix& bdmat, Matrix& bsmat, double& detJ);
private:
    Oceane::Vector m_xcoord;
    Oceane::Vector m_ycoord;

};

} // namespace Oceane

#endif // RECT4_H
