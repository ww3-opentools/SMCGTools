/* This Arctics.h header file declare a group of parameters used for polar 
 * part(s) and definition of map-east reference direction.   JGLi06Jun2026  */

#pragma once
#include <vector>
#include "WvSpectr.h"
#include "GtCrTurn.h"

struct Arctics {

    int NCel, NGlo, NArc, NGlB, NArB, NPol, South, North;
    float PLat=0.0, PLon=-180.0, PCRd;
    std::vector<float> Angld, ELats;
    std::vector<int>  MBGlo, MBArc; 

};

const Arctics&  getArctics();
int ArctAngd( bool TstOut = false );
int Specturn( std::vector<float>& spct, std::size_t NFrq, std::size_t NDrc,
              float AnglD );
int ArctBndy( WvSpectr& Spc, const Arctics& Arc, int NT, bool TstOut = false );
int ArcScalr( std::vector<float>& Sclr, const Arctics& Arc );
int ArcVectr( std::vector<float>& Uctr, 
              std::vector<float>& Vctr, const Arctics& Arc );

// End of Arctics.h.

