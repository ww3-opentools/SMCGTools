/* The header file WvCgrpKs.h to declare the WvCgrpKs struct and its related 
 * functions, including extracting wave number k and group speed at any sea 
 * point.                                                    JGLi06Jul2026  */

#ifndef WvCgrpKs_H
#define WvCgrpKs_H

#include <vector>
#include <cstddef>

struct WvCgrpKs {
    std::size_t nsea, nfrq;
    std::vector<std::vector<float>> WnmK, Cgrp, Rfrc;

    // Dimension initialization constructor.
    WvCgrpKs(std::size_t NSea, std::size_t NFrq ):
        nsea(NSea), nfrq(NFrq),
        WnmK(NSea, std::vector<float>(NFrq+1)),
        Cgrp(NSea, std::vector<float>(NFrq+1)),
        Rfrc(NSea, std::vector<float>(NFrq+1)) {}

};

// Function to initialize wave number and group speed.
int WvCGrpKIni( WvCgrpKs& CgK );

#endif

