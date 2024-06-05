#ifndef INITIALIZE_H
#define INITIALIZE_H

#include "constraints.h"

#include <vector>
#include <string>

class Initialize {
public:
    Initialize();

    void startInitializing();

private:
    std::string filename;
};

#endif // READ_H