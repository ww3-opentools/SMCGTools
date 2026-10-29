/* This SMCGradn.h header declared a group of parameters for the SMCGradn.cpp 
 * module to calculate gradients of a variable on SMC grid.  JGLi06Jun2026  */

#pragma once
#include <vector>

struct Gradnts {

    int Nmbr; 
    std::vector<float> varb, dvdx, dvdy;
    Gradnts(int nmbr): Nmbr(nmbr), 
        varb(nmbr), dvdx(nmbr), dvdy(nmbr) {}

};

int SMCGradn( Gradnts& gdnt, bool L0r1 = true );

// End of SMCGradn.h.

