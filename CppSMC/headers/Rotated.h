/* This Rotated.h header file declare a group of working vectors for the 
 * conversion functions between standard and rotated lat-lon grids.         */ 

#pragma once
#include <vector>

// namespace Rtd {     // Only use the namespace scope if clash may occur.

struct RotateArrays {

    float polat, polon;
    std::vector<float>  slat, slon, elat, elon, angd;
    RotateArrays(int n): slat(n), slon(n), elat(n), elon(n), angd(n) {}

};

int lltoeqangled(RotateArrays& wk, int Numbr);
int eqtollangled(RotateArrays& wk, int Numbr);

// }    // Suspended namespace scope.

