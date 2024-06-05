#include "read.h"
#include "initialize.h"
#include "solve.h"
#include <iostream>

int main() 
{
    Read reader;
    reader.startreading();
    Initialize init;
    init.startInitializing();
    Solve solver1;
    solver1.Solver();

    return 0;
}
