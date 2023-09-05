#ifndef MATERIAL_H
#define MATERIAL_H

#include <memory>
#include "matrix.h"

namespace Oceane {

class Material
{
public:
    Material();
    virtual ~Material();

};

using Materialptr = std::shared_ptr<Material>;

} // namespace Oceane
#endif // MATERIAL_H
