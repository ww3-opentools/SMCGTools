/* The CGrpKIni module initializes the wave number K and group speed fields 
 * at all sea points and frequencies and store in WvNumbrK.   JGLi09Jun2026 */

#include <omp.h>
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "SMCGrid.h"
#include "Working.h"
#include "Spectrum.h"
#include "WvCgrpKs.h"

#include "SMConsts.h"
using namespace SMConsts;

// Function to initialize SMC grid cell wave numbers and group speed.
int WvCGrpKIni( WvCgrpKs& CgK ) {
    
    // Work out wave number and group speed dimensions.
    auto NSea = CgK.nsea;
    auto NFrq = CgK.nfrq;
    
    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    auto NCel = grd.ncel;
    cout << " Total number of cells grd.ncel = " << NCel << endl;
    auto zlon = grd.zlon;
    auto zlat = grd.zlat;
    auto dlon = grd.dlon;
    auto dlat = grd.dlat;

    //  Check nlvl and npol from grd struct.
    auto NLvl = grd.nlvl; 
    auto NPol = grd.npol;
    cout << " Number of levels = " << NLvl << endl;
    cout << " Number of poles  = " << NPol << endl;

    // Factor of maxium level cell size.
    int MFct = 1 << (NLvl - 1);
    cout << " Maximum level MFC = " << MFct << endl;
 
    // Fetch NSea and NFrq from Spectrum spm strucut.
    const auto& spm = getSpectrum();   // Read only spectral parameters.
    auto sgma = spm.sigma;

    // Check CgK dimensions matching with grd and spm ones.
    if( NSea < NCel +10 ) cout <<" *** Warning NSea less than NCel +10 *** \n";
    if( NFrq < spm.nfrq ) cout <<" *** Warning NFrq less than spm.nfrq *** \n";

    // Assign 0.0 to working variable of required dimensions.
    // Note assign first resize the vector if not present or overwrite it.
    cnst = static_cast<float>( 0.5*WtrDepthMn );
    wrk0.assign(NSea, cnst); 

    // Pass cell water depth to wrk0, minimum depth WtrDepthMn. 
    cnst1 = static_cast<float>( WtrDepthMn );
    for(int n=0; n<NCel; n++) {
        wrk0[n] = grd.CEL[n][4];
        if( wrk0[n] < cnst1 ) wrk0[n] = cnst1;
    }

    // Precision for wave number approximation.
    cnst0 = 1.0E-6f;
    cout << " Calculating WnmK with precision = " << cnst0 << endl;

    // Nested loops over frequency and sea points.
    int itk, itkmx = 0;
    for(int n=0; n<NCel; n++) {
        cnst3 = wrk0[n];
        for(int k=0; k<NFrq+1; k++) {
            cnst  = sgma[k]; 
            cnst4 = static_cast<float>( cnst*cnst/Gravty );

            // Intitial values for wave number K iteration.
            cnst1 = cnst4/tanh( cnst4*cnst3 );
            cnst2 = cnst4/tanh( cnst1*cnst3 );
            itk = 0;
            while( abs( cnst2-cnst1 ) > cnst0 ) {
                cnst1 = cnst2;
                cnst2 = cnst4/tanh( cnst1*cnst3 );
                itk += 1;
            }
            // Save the final wave number K value.
            CgK.WnmK[n][k] = static_cast<float>( cnst2 );
            itkmx = max(itk, itkmx);

            // Group speed for this frequency.
            cnst5 = cnst2*cnst3;
            cnst6 = 1.0/cosh( cnst5 );
            CgK.Cgrp[n][k] = (Gravty*0.5/cnst)*(cnst6*cnst6*cnst5 + tanh(cnst5));

            // Refraction coefficient uses minimum depth of 30 m.
            cnst1 = max(30.0f, cnst3);
            CgK.Rfrc[n][k] = cnst/sinh(2.0*cnst2*cnst1);
            
        } // End of frequency NFrq+1 loop.
//      if((n % 5000) == 0) cout <<" Wnmk calculated for nsea = "<< n << endl;

    } // End of cell loops NCel.
    cout << " Maximum K iteratation number itkmx = " << itkmx << endl;

    // Assign Cgrp and Rfrc = 0 for boundary cells n>=NCel.
    if( NSea > NCel ) {
        for(int n=NCel; n<NSea; n++) {
            std::fill(CgK.Cgrp[n].begin(), CgK.Cgrp[n].begin() + NFrq+1, 0.0f);
            std::fill(CgK.Rfrc[n].begin(), CgK.Rfrc[n].begin() + NFrq+1, 0.0f);
        }
        cout << " Boundary cell Cgrp Rfrc set to zero CgK.Cgrp[NCel][NFrq] = " 
                                              << CgK.Cgrp[NCel][NFrq] << endl;
    } 

    return(0); 

} // End of function WvCGrpKIni.cpp.

