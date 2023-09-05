#include "elementfactory.h"

namespace Oceane {

ElementFactory::ElementFactory()
{

}

Elementptr ElementFactory::createElement(std::string type,Index id, std::vector<Nodeptr> nodes)
{
    if(type=="CPS4")
        return std::make_shared<CPS4>(id,nodes);
    if(type=="CPS4P")
        return  std::make_shared<CPS4P>(id,nodes);
    if(type=="RECT4")
        return std::make_shared<Rect4>(id,nodes);
}


} // namespace Oceane
