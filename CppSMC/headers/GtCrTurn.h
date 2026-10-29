/* This GtCrTurn.h header file declares a group of parameters used for the  
 * great circle turning and refraction at all sea points.    JGLi03Jul2026  */

#pragma once
#include <vector>
#include "WvSpectr.h"
#include "WvCgrpKs.h"

struct GtCrTn {
    std::size_t nsea, ndrc;
    std::vector<float> cthr, dhdx, dhdy;
    std::vector<std::vector<float>> DHLmt, DHThn;
    GtCrTn(std::size_t NSea, std::size_t NDrc): 
        nsea(NSea), ndrc(NDrc), 
        cthr(NSea), dhdx(NSea), dhdy(NSea),
        DHLmt(NSea, std::vector<float>(NDrc)), 
        DHThn(NSea, std::vector<float>(NDrc)) {}
};

int GtCrTIni( GtCrTn& Cfr );
int GtCrTurn( WvSpectr& Spc, GtCrTn& Cfr, WvCgrpKs& CgK, float DTLc );

// End of GtCrTurn.h.

