/* The GtCrTurn module performs the directional bin turning for rotations in 
 * great circle turning and refraction induced deflection.   JGLi29Jun2026  */

#include <omp.h>
#include <iostream>
using namespace std;
#include <cmath>
#include <string>
#include <vector>
#include <optional>
#include <algorithm>

#include "SMCGrid.h"
#include "Arctics.h"
#include "Working.h"

#include "SMCGradn.h"
#include "Spectrum.h"
#include "WvSpectr.h"
#include "WvCgrpKs.h"
#include "GtCrTurn.h"
#include "SMConsts.h"
using namespace SMConsts;

// Function to initialize Great circle turning parameters.
int GtCrTIni( GtCrTn& Cfr ) {
    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    auto NCel = grd.ncel;
    auto NGlo = grd.nglo;
    auto NPol = grd.npol;
    auto NSea = Cfr.nsea;
    auto NDrc = Cfr.ndrc;
    if( NSea != NCel ) {
        cout << NSea <<" from Cfr.nsea not equal to grd.ncel = " << NCel << endl;
        exit(1);
    }

    // Fetch the Arctic parameters and store in optional variables.
    std::vector<float> ELats, Angld;
    std::size_t NArc;
    if( grd.Arctic ) {
        const Arctics& Arc = getArctics();
        NArc = Arc.NArc;
        ELats = Arc.ELats;
        Angld = Arc.Angld;
        cout << " * NArc and NGlo = " << Arc.NArc <<" "<< Arc.NGlo << endl;
        cout << " * NGlB and NArB = " << Arc.NGlB <<" "<< Arc.NArB << endl;
        if( Arc.NGlo != NGlo ) cout <<" *** Warning: Arc.NGlo not matching grd.nglo ! \n";
    }

    // Define great circle turning factor for all except for polar cells.
    auto zlat = grd.zlat;
    auto dlat = grd.dlat;
    float hfdlat = 0.5*dlat;
    float rearthinv = 1.0/REarth;

    #pragma omp parallel for default(none) \
         shared(NGlo, Cfr, grd, zlat, dlat, hfdlat, rearthinv) 
    for(int i=0; i<NGlo; i++) {
        celtyp cell = grd.CEL[i];
        float radlat = (zlat + (2*cell[1] + cell[3])*hfdlat)*D2Rad;
        Cfr.cthr[i] = -rearthinv * tan(radlat);
    }

    if( grd.Arctic ) {
         #pragma omp parallel for default(none) \
              shared(NGlo, NArc, Cfr, ELats, rearthinv)
        for(int i=0; i<NArc; i++) {
            float radlat = ELats[i]*D2Rad;
            Cfr.cthr[i+NGlo] = -rearthinv * tan(radlat);
        }
    }
    cout << " GtCrTn Cfr.cthr defined with NSea = " << NSea << endl;

    // Fetch NDrc and DThta from Spectrum spm strucut.
    const auto& spm = getSpectrum();   // Read only spectral parameters.
    auto DThta = spm.DThta;
    auto Theta = spm.theta;
    auto CThta = spm.cthta;
    auto SThta = spm.sthta;

    if( spm.ndrc != NDrc ) {
        cout << NDrc <<" from Cfr.ndrc not equal to spm.ndrc = " << spm.ndrc << endl;
        exit(1);
    }

    // Calculate SMC grid bathymetry gradients.
    Gradnts Gdnt(NCel);  // Allocate Gradnts struct.
 
    // Pass the cell depth to Gdnt.varb for floor gradient evaluation.
    #pragma omp parallel for default(none) \
         shared(NCel, Gdnt, grd)
    for(int i=0; i<NCel; i++ ) {
        Gdnt.varb[i] = grd.CEL[i][4];   //Default integer conversion to float.
    }

    cout << " Calling SMCGradn module ... \n";
    SMCGradn( Gdnt );
    
    // Save the bathymetry gradient to Cfr.dhdx and dhdy.
    Cfr.dhdx = Gdnt.dvdx;
    Cfr.dhdy = Gdnt.dvdy;
    cout << " Bathymetry gradient saved to Cfr.dhdx/y. \n"; 

    // Calculate refraction rotation limit by bathymetry gradient.
    float Angmx = Rotatn_DMx*D2Rad;   // Maximum refraction rotation angle in radian.
    #pragma omp parallel for default(none) \
         shared(NGlo, NDrc, Cfr, grd, SThta, CThta, Theta, DThta, Angmx) \
         private(cnst1, cnst2, cnst3, cnst4, cnst5, cnst6)
    for(int n=0; n<NGlo; n++ ) {
        cnst1 = Cfr.dhdx[n];
        cnst2 = Cfr.dhdy[n];
        cnst4 = sqrt( cnst1*cnst1 + cnst2*cnst2 );
        for(int i=0; i<NDrc; i++) {
            Cfr.DHThn[n][i] = cnst1*SThta[i] - cnst2*CThta[i];
            if( cnst4 > 1.0E-6f ) {
                cnst5 = cnst1*CThta[i] + cnst2*SThta[i];
                cnst6 = acos( -cnst5/cnst4 );
                cnst3 = std::min( cnst6, static_cast<float>(Pie-cnst6) );
                Cfr.DHLmt[n][i] = std::min( Angmx, cnst3 )/DThta;
            }
            else   Cfr.DHLmt[n][i] = 0.0;
        }
    } // End of NGlo loop.

    if( grd.Arctic ) {
         #pragma omp parallel for default(none) \
              shared(NGlo, NDrc, NArc, Cfr, ELats, Angld, Theta, DThta, Angmx) \
              private(cnst0, cnst1, cnst2, cnst3, cnst4, cnst5, cnst6)
        for(int n=0; n<NArc; n++ ) {
        cnst1 = Cfr.dhdx[n+NGlo];
        cnst2 = Cfr.dhdy[n+NGlo]; 
        cnst4 = sqrt( cnst1*cnst1 + cnst2*cnst2 );
        for(int i=0; i<NDrc; i++) {
            cnst0 = Theta[i] - Angld[n]*(float)D2Rad;
            Cfr.DHThn[n+NGlo][i] = cnst1*sin(cnst0) - cnst2*cos(cnst0);
            if( cnst4 > 1.0E-6f ) {
                cnst5 = cnst1*cos(cnst0) + cnst2*sin(cnst0);
                cnst6 = acos( -cnst5/cnst4 );
                cnst3 = std::min( cnst6, static_cast<float>(Pie-cnst6) );
                Cfr.DHLmt[n+NGlo][i] = std::min( Angmx, cnst3 )/DThta;
            }
            else   Cfr.DHLmt[n+NGlo][i] = 0.0;
        }
        } // End of NArc loop.
    } // End if( Arctic ) block.
    cout << " Bathymetry refraction Cfr.DHLmt and DHTha defined. \n"; 
    
    return(0);
  
} // End of function GtCrTIni. 


// Function to rotate directional bins by a given angle.  
int GtCrTurn( WvSpectr& Spc, GtCrTn& Cfr, WvCgrpKs& CgK, float DTLc ) {
//  cout <<" GtCrTurn starts .... \n";
    //  Find out spectral dimensions.
    auto NSea = Spc.NSeapt();
    auto NDrc = Spc.NDrctn();
    auto NFrq = Spc.NFrqcy();
    auto MSpc = Spc.MSpctr();

    // Use a copy of spectrum for GtCrTurning.
    std::vector<float> spct(MSpc);

//  cout << " Number of sea points in the grid = " << NSea << endl;
//  cout << " Number of spectral direction bin = " << NDrc << endl;
//  cout << " Number of spectral frequency bin = " << NFrq << endl;

    if( NDrc == 0 ) {
        cout << " *** No directional dimension! \n"; 
        exit(1);
    }
    // Use Spectrum parameters.
    const auto& spm = getSpectrum();   // Read only spectral parameters.
//  cout << " * Spectral parameters set up spm.DThta = " << spm.DThta << endl;
    auto CThta = spm.cthta;
    auto SThta = spm.sthta;
    float cthr0 = DTLc/spm.DThta; 
    auto ndrc = static_cast<int>( NDrc );
    int j, k;

    // Loop over all Nsea point in wave spectra Spc.
    for(int ns=0; ns<NSea; ns++) {

        // Reference to the spectrum MSpc elements for site isea.
        spct = Spc[ns];   

    // Loop over all frequency range.
    for (int m=0; m<NFrq; m++) {

        // Copy directional bins into wrk0[0:NDrc]. 
        std::copy(spct.begin() + NDrc*m, 
                  spct.begin() + NDrc*(m+1), wrk0.begin());
        wrk1.assign(NDrc, 0.0f); 
        float cflgc = cthr0*Cfr.cthr[ns]*CgK.Cgrp[ns][m];
        float cflrf = cthr0*CgK.Rfrc[ns][m];

        // Loop over all directional bins for fraction or s + 10 boundary cells.
//      The OpenMP directive lines suspended as it makes the function delay.
//      #pragma omp parallel for default(none) \
//           shared(m, ns, ndrc, Cfr, cflgc, cflrf, CThta, wrk0, wrk1) \
//           private(cnst1, cnst2, cnst3, cnst4, cnst5, cnst6, j, k)
        for (int i=0; i<ndrc; i++) {

            // Work out bathymetry refraction angle for this direction bin.
            cnst3 = cflrf*Cfr.DHThn[ns][i];
            cnst4 = copysign( min( abs(cnst3), Cfr.DHLmt[ns][i] ), cnst3);

            // Copy the GCT and Refraction combined turning Courant number.
            cnst6 = cflgc*CThta[i] + cnst4;  

            // Work out absolute turning bin numbers.
            cnst5 = abs( cnst6 );
            int L = static_cast<int>( cnst5 );
            cnst1 = cnst5 - L;
            // Remove full circle turning if any.
            L %= ndrc; 

            // Partitioning coef for neighbouring bins.
            cnst2 = 1.0f - cnst1;
            
            // Positive or anti-clock turning.
            if( cnst6 > 0.0 ) {
                j = i + L;
                k = j + 1;
                if( j >= ndrc ) j -= ndrc;
                if( k >= ndrc ) k -= ndrc;
            }
            else {
                j = i - L;
                k = j - 1;
                if( j < 0 ) j += ndrc;
                if( k < 0 ) k += ndrc;
            }

            // Partitioning i bin into j and k bins.
//          #pragma omp atomic update
            wrk1[j] += wrk0[i]*cnst2;
//          #pragma omp atomic update
            wrk1[k] += wrk0[i]*cnst1;

        }  // End of directional bin i loop. 

        // Copy rotated directional bins wrk1 back to Spc directly. 
        std::copy(wrk1.begin(), wrk1.begin() + NDrc, Spc[ns].begin() + NDrc*m);

    }   // End of frequence n-loop
    }   // End of sea point ns-loop

//  cout <<" GtCrTurn completed. \n";

    return(0); 

    // End of GtCrTurn function.
}

// End of module GtCrTurn.cpp file.

