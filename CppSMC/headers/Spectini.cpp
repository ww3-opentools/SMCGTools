/* The Spectini module initializes the wave spectral field and stores it in 
 * the class WvSpectr storage space.                          JGLi09Jun2026 */

#include <omp.h>
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "WvSpectr.h"
#include "Spectrum.h"
#include "SMCGrid.h"
#include "Working.h"
#include "Propgtn.h"

#include "SMConsts.h"
using namespace SMConsts;

// Function to initialize SMC grid cell wave spectral field. 
int Spectini( WvSpectr& Spc ) {
    
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
 
    // Fetch NDrc and NFrq from Spc class.
    auto NSea = Spc.NSeapt();
    auto NFrq = Spc.NFrqcy();
    auto NDrc = Spc.NDrctn();
    auto MSpc = Spc.MSpctr();

    if( NSea != NCel ) cout <<" *** Nsea not matching NCel *** \n";
    cout << " Total number of NSea points = " << NSea << endl;

    // Fetch NDrc and DThta from Spectrum spm strucut.
    const auto& spm = getSpectrum();   // Read only spectral parameters.
    auto DThta = spm.DThta;
    auto Theta = spm.theta;
    auto CThta = spm.cthta;
    auto SThta = spm.sthta;
    auto DSgma = spm.dsgma;
    auto Sigma = spm.sigma;

    if( spm.ndrc != NDrc or spm.nfrq != NFrq ) {
        cout << NDrc <<" from Spc not equal to spm.ndrc = " << spm.ndrc << endl;
        cout << NFrq <<" from Spc not equal to spm.nfrq = " << spm.nfrq << endl;
        exit(1);
    }

    // Assign 0.0 to working variable of required dimensions.
    // Note assign first resize the vector if not present or overwrite it.
    wrk1.assign(NDrc, 0.0f); 
    wrk2.assign(NDrc, 0.0f); 
    wrk3.assign(MSpc, 0.0f); 
    wrk4.assign(MSpc, 0.0f); 

    // Initial wave spectra at 45 degree direction.
    float dth = (float)TPie/(float)NDrc;
    cnst = 128.0f/(float)Pie;
    cnst6 = 0.0f;
    for(int k=0; k<NDrc; k++) {
        cnst1 = cos(k*dth + Pie/4.0f);
        cnst2 = cos(k*dth - Pie/4.0f);
        if( cnst1 > 0.0 ) wrk1[k]=cnst*cnst1*cnst1;
        if( cnst2 > 0.0 ) wrk2[k]=cnst*cnst2*cnst2;
        cnst6 += wrk1[k];
    }
    float SWH0 = sqrt( cnst6*dth );

    // Frequecy variable will be defined using the well-developed 
    // Pierson-Moskowitz (P-M) spectrum but normalised for the selected 
    // frequency range with a weighting factor, E0.
    // as E0*(x^5)*EXP(-beta*x^4), where x=Ompeak/omega and
    // beta = 1.25 is a constant. 
    std::vector<float> SpctPM(NFrq);
    float E0PM;
    // Single frequency option.
    if( NFrq == 1 ) {
        wrk3 = wrk1;
        wrk4 = wrk2;
    }
    else {
        cnst0 = 0.0f;
        cnst2 = 1.25f;
        cnst3 = Sigma[NFrq/2];
        for(int m=0; m<NFrq; m++) {
            cnst1 = cnst3/Sigma[m];
            cnst4 = cnst1*cnst1*cnst1*cnst1;
            SpctPM[m] = (cnst1*cnst4)*exp(-cnst2*cnst4);
            cnst0 += SpctPM[m]*DSgma[m];
        }
        // Normalisation factor for selected range.
        E0PM = 1.0f/cnst0;

        for(int m=0; m<NFrq; m++) {
            for(int i=0; i<NDrc; i++) {
                auto ks = i + m*NDrc;
                wrk3[ks] = wrk1[i]*E0PM*SpctPM[m];
                wrk4[ks] = wrk2[i]*E0PM*SpctPM[m];
            }
        }
    } // if-else block.

    // Assign spectral to Spc if points are within belt range.
    int jnor = static_cast<int>( ( 45.0f - zlat)/dlat );
    int jnnd = static_cast<int>( ( 75.0f - zlat)/dlat );
    int jpol = static_cast<int>( ( 86.0f - zlat)/dlat );

    // Add Atlantic round patch for comparison with Arctic one.
    int ijk = static_cast<int>( (330.0f - zlon)/(dlon*MFct) )*MFct + MFct/2;
    int lmn = static_cast<int>( (  5.0f - zlat)/(dlat*MFct) )*MFct + MFct/2;
    cout <<" Atlanctic patch centre i, j = "<< ijk <<" "<< lmn << endl;

    celtyp cell;
    for (int n=0; n<NCel; n++) {
        cell = grd.CEL[n];

        if( (cell[1] >= jnor and cell[1] < jnnd) or cell[1] < -jpol ) {
            Spc[n] = wrk3; 
        }
        if( (cell[1] < -jnor and cell[1] > -jnnd) or cell[1] >= jpol ) {
            Spc[n] = wrk4; 
        }
        // Atlantic patch within a cirecle of 5 degree.
        int kk = cell[0] + cell[2]/2 - ijk;
        int nn = cell[1] + cell[3]/2 - lmn;
        cnst3 = dlon*kk/5.0f;
        cnst4 = dlat*nn/5.0f;
        if( (cnst3*cnst3 + cnst4*cnst4) < 1.01f ) {
            Spc[n] = wrk4; 
        }
    } // End of n<NCel loop.

    return(0); 

    // End of Spctrini function.
}

