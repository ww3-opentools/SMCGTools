#ifndef CPP_SMC1d_H
#define CPP_SMC1d_H

#include <string>

/***************************************************************************
 *  The header declares dimension constants for SMC1d grid.
 *  To be used by #include "SMC1d.h".        JGLi16Jun2026  
 ***************************************************************************/

// namespace SMCGrid {
    const int  NCM = 36200, NFM = 36400; // Cell & face max numbers.
    const int  NLV = 1, NPol = 2;   // Number of levels and polar cells.
    constexpr int nlv1= NLV + 1;    // nlv1 for count line dimenion = NLV + 1.
    const float zlon = 0.0, zlat = -0.5, dlon = 1.40625, dlat = 1.00;
    bool  Arctic = false;

    std::string indir = "/data/users/jianguo.li/SMCGTools/DatS1d/"; 
    std::string model = "SMC1d";  // C++ style string.
// }

/* End of SMC1d.h  */

#endif

