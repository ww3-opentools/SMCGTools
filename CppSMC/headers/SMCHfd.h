#ifndef CPP_SMCHfd_H
#define CPP_SMCHfd_H

#include <string>

/***************************************************************************
 *  The header declares dimension constants for SMCHfd grid.
 *  To be used by #include "SMCHfd.h".        JGLi11Nov2026  
 ***************************************************************************/

// namespace SMCGrid {
    const int  NCM = 145000, NFM = 146000; // Cell & face max numbers.
    const int  NLV = 1, NPol = 2;   // Number of levels and polar cells.
    constexpr int nlv1 = NLV + 1;   // nlv1 for count line dimenion = NLV + 1.
    const float zlon = 0.0, zlat = -0.5, dlon = 0.703125, dlat = 0.50;
    bool  Arctic = true;

    std::string indir = "/data/users/jianguo.li/SMCGTools/DatHfd/"; 
    std::string model = "SMCHfd";  // C++ style string.
// }

/* End of SMCHfd.h  */

#endif

