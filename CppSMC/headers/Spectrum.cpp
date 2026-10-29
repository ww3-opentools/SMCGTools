/****************************************************************************
 *  This Spectrum.cpp program set up the spectral related parameters, which 
 *  is set once and quoted by other modules.       JGLi03Jul2026  
 ****************************************************************************/

#include <iostream>
#include <iomanip>
using namespace std;
#include <stdexcept>
#include <cmath>

#include "SMConsts.h"
using namespace SMConsts;

#include "Spectrum.h"

namespace {     // Anonymous namespace only available to this file.
    struct Spectrum spm;
    bool initialized = false;
}

/* A function to set up spectrum related parameters.  */
void setSpectrum( std::size_t NFrq, std::size_t NDrc ) { 
    if( initialized ) return;

    // Resize spectral variables.
    spm.nfrq = NFrq;
    spm.ndrc = NDrc;
    spm.sigma.resize(NFrq+1);
    spm.dsgma.resize(NFrq+1);
    spm.theta.resize(NDrc);
    spm.cthta.resize(NDrc);
    spm.sthta.resize(NDrc);

    // Define radian frequcies.
    float dsfct = 0.5f*(spm.Frqxf - 1.0f/spm.Frqxf);
    if( NFrq == 1 ) {
        spm.sigma[0] = static_cast<float>( spm.Frqc1 * TPie ); 
        spm.sigma[1] = static_cast<float>( spm.Frqc1 * TPie * spm.Frqxf ); 
        spm.dsgma[0] = 1.0f;
        spm.dsgma[1] = 1.0f;
    }
    else {
        spm.sigma[0] = static_cast<float>( spm.Frqc0 * TPie );
        spm.dsgma[0] = spm.sigma[0] * dsfct;
        for(std::size_t n=0; n<NFrq; n++ ) {
            spm.sigma[n+1] = spm.sigma[n] * spm.Frqxf;
            spm.dsgma[n+1] = spm.sigma[n+1] * dsfct;
        }
    }
     
    // Define directional bins and trianglar function values. 
    if( NDrc == 1 ) {
        cout <<" *** Only one direction bin! *** \n";
        spm.Thta0 = 0.0f;
        spm.DThta = 1.0f;   
        spm.theta[0] = 0.0f;
        spm.cthta[0] = 1.0f;
        spm.sthta[0] = 0.0f;
    }
    else {
        spm.DThta = static_cast<float>( TPie/NDrc );
        spm.Thta0 = 0.5f * spm.DThta;
        for(std::size_t m=0; m<NDrc; m++ ) {
            spm.theta[m] = spm.Thta0 + spm.DThta*m; 
            spm.cthta[m] = cos(spm.theta[m]);
            spm.sthta[m] = sin(spm.theta[m]);
        }
    }
    // Reset initialized flag. 
    initialized = true;

} // End of function setSpectrum.


/* A function to fetch struct Spectrum for shared use. */
const Spectrum& getSpectrum() {
    if (!initialized) {
        throw std::runtime_error("*** Spectrum data not initialized.");
    }
    // Return the Spectrum spm data, otherwise.
    return spm;

} // End of getSpectrum function.

// End of Spectrum.cpp module.

