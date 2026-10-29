#ifndef CPP_SMC61250_H
#define CPP_SMC61250_H

#include <string>

/***************************************************************************
 *  The header declares dimension constants for SMC61250 grid.
 *  To be used by #include "SMC61250.h"       JGLi28Apr2026  
 ***************************************************************************/

// namespace SMCGrid {
    const long int  NCM = 250000L, NFM = 290000L; // Cell & face max numbers.
    const int  NLV = 4, NPol = 1;   // Number of levels and polar cells.
    constexpr int nlv1= NLV + 1;    // nlv1 for count line dimenion = NLV + 1.
    const float zlon = 0.0, zlat = 0.0, dlon = 0.087890625, dlat = 0.058593750;
    bool  Arctic = true;

    std::string indir = "/data/users/jianguo.li/SMCGTools/DatGMC/"; 
    std::string model = "SMC61250";  // C++ style string.
// }

/* End of SMC61250.h  */

#endif

