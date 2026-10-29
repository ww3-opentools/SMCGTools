#ifndef CPP_SMC1R3_H
#define CPP_SMC1R3_H

#include <string>

/***************************************************************************
 *  The header declares dimension constants for SMC1R3 grid.
 *  To be used by #include "SMC1R3.h"       JGLi23Jul2026  
 ***************************************************************************/

// namespace SMCGrid {
    const long int  NCM = 59000L, NFM = 60000L; // Cell & face max numbers.
    const int  NLV = 3, NPol = 2;   // Number of levels and polar cells.
    constexpr int nlv1= NLV + 1;    // nlv1 for count line dimenion = NLV + 1.
    const float zlon = 0.0, zlat = 0.0, dlon = 0.3515625, dlat = 0.2500;
    bool  Arctic = true;

    std::string indir = "/data/users/jianguo.li/SMCGTools/DatSWE/"; 
    std::string model = "SMC1R3";  // C++ style string.
// }

/* End of SMC1R3.h  */

#endif

