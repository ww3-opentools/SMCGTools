/* The SWHeight module integrates wave spectral to the SWHt for every grid   
 * point and store it into a vector variable for printing.   JGLi08Jul2026  */

#include <omp.h>
#include <iostream>
using namespace std;
#include <cmath>
#include <string>
#include <vector>
#include <algorithm>

#include "SMCGrid.h"
#include "Working.h"

#include "Spectrum.h"
#include "WvSpectr.h"
#include "SMConsts.h"
using namespace SMConsts;

// Function to initialize Great circle turning parameters.
std::vector<float> SWHeight( WvSpectr& Spc ) {
    // Fetch the spectral array parameters. 
    auto NSea = Spc.NSeapt();
    auto NFrq = Spc.NFrqcy();
    auto NDrc = Spc.NDrctn();
    auto MSpc = Spc.MSpctr();

    // Declare SWHs vector and initialize to be zero.
    std::vector<float> SWHs(NSea);
    SWHs.assign(NSea, 0.0f);

    // Fetch the spectrum struct for spectral definitions.
    const auto& spm = getSpectrum();   
    auto DThta = spm.DThta;
    auto DSgma = spm.dsgma;
    auto ndrcs = spm.ndrc;
    auto nfrqs = spm.nfrq;
    if( ndrcs != NDrc || nfrqs != NFrq ) {
        cout << " input Spc not matching spm parameters \n"; 
        return SWHs;
    }

    // Integrate spectrum at every point.
    #pragma omp parallel for default(none) \
         shared(Spc, NSea, NFrq, NDrc, DThta, DSgma, SWHs) 
    for(int n=0; n<NSea; n++) {
        auto& spct = Spc[n];
        float engy = 0.0f;
        for(int k=0; k<NFrq; k++) {
            for(int i=0; i<NDrc; i++) {
                engy += spct[k*NDrc + i]*DSgma[k];
            }
        }
        float swht = sqrt( abs(engy*DThta) );
        SWHs[n] = copysign(swht, engy); 
    }
    return SWHs;
  
} // End of function SWHeight.cpp module.


