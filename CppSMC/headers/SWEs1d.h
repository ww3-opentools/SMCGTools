#ifndef CPP_SWEs1d_H
#define CPP_SWEs1d_H

#include <string>

/***************************************************************************
 *  The header declares dimension constants for SWEs1d grid.
 *  To be used by #include "SWEs1d.h".        JGLi28Apr2026  
 ***************************************************************************/

// namespace SMCGrid {
    const int  NCM = 37000, NFM = 38000; // Cell & face max numbers.
    const int  NLV = 1, NPol = 2;   // Number of levels and polar cells.
    constexpr int nlv1= NLV + 1;    // nlv1 for count line dimenion = NLV + 1.
    const float zlon = 0.0, zlat = 0.0, dlon = 1.40625, dlat = 1.00;
    bool  Arctic = true;

    std::string indir = "/data/users/jianguo.li/MFTModel/SWEs1d/DatGMC/"; 
    std::string model = "SWEs1d";  // C++ style string.
// }

/* End of SWEs1d.h  */

#endif

