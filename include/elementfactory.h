#ifndef ELEMENTFACTORY_H
#define ELEMENTFACTORY_H

#include "element.h"
#include "quad4.h"
#include "rect4.h"

namespace Oceane {



class ElementFactory
{
public:
    ElementFactory();
    Elementptr createElement(std::string type,Index id, std::vector<Nodeptr> nodes);
    Elementptr createElement(std::string type, Index id,std::vector<Index> nodeids);
};


} // namespace Oceane

#endif // ELEMENTFACTORY_H
