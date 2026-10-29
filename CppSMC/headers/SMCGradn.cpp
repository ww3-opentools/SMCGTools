/* The SMCGradn.cpp module calculates gradients of any variable defined on a 
 * SMC grid with given cell and side arrays.                 JGLi08Jun2026  */

#include <omp.h>
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "SMCGrid.h"
#include "Working.h"
#include "SMCGradn.h"

#include "SMConsts.h"
using namespace SMConsts;

// Function to calculate variable gradient on a SMC grid. 
int SMCGradn( Gradnts& Gdnt, bool L0r1 ) {
    
    //  Reset Gdnt.Nmbr if not matching Gdnt.varb dimension.
    if( Gdnt.varb.size() != Gdnt.Nmbr ) {
        cout << " Gdnt.Nmbr = " << Gdnt.Nmbr << endl;
        Gdnt.Nmbr = Gdnt.varb.size();
        cout << " reset to Gdnt.varb.size() = " << Gdnt.Nmbr << endl;
    }

    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    const auto& gsd = getSMCSID();   // Get gsd as read only const reference.
    auto NPol = grd.npol;
    cout << " Number of polar cells grd.npol = " << NPol << endl;

    // Total cell number should match Gdnt.Nmbr
    auto NCel = grd.ncel;
    if( Gdnt.Nmbr != NCel ) {
        cout << " Gdnt.Nmbr = " << Gdnt.Nmbr << endl;
        cout << " not matchin grd.ncel = " << NCel << endl;
        return(1);
    }
 
    // Confirm boundary condition to be first (true) or second (false) type.
    if( L0r1 ) cout << " First type boundary with zero boundary cells. \n";
    else  cout << " Second type boundary condition with zero-gradient. \n";
 
    // Assign 0.0 to working variable of dimension NCel + 10 boundary cells.
    // Note assign first resize the vector if not present or overwrite it.
    double small = 1.0E-6;
    int NCB = NCel + 10;
    wrk0.assign(NCB, 0.0f); 
    wrk1.assign(NCB, 0.0f); 
    wrk2.assign(NCB, 0.0f); 
    wrk3.assign(NCB, 0.0f); 
    wrk4.assign(NCB, 0.0f); 
    
    // Copy Gdnt.varb to wrk1[0:NCel] and set boundary values zero.
    std::copy(Gdnt.varb.begin(), Gdnt.varb.begin() + NCel, wrk0.begin());
    for (int i=NCel; i<NCB; i++) wrk0[i] = 0.0;

    cout << " Grnt.varb passed to wrk1[0:NCel] \n";
    cout << " wrk0[0] = " << wrk0[0] <<"\n";
    cout << " wrk0[NCel-1] = " << wrk0[NCel-1] <<"\n";
    cout << " wrk0[NCB -1] = " << wrk0[NCB -1] <<"\n";

    // Fetch some constants for gradient calculation.
    auto zlon = grd.zlon;
    auto zlat = grd.zlat;
    auto dlon = grd.dlon;
    auto dlat = grd.dlat;
    double DX0I = 1.0/(dlon*D2Rad*REarth);
    double DY0I = 1.0/(dlat*D2Rad*REarth);
    auto NUFc = gsd.nufc;
    auto NVFc = gsd.nvfc;

    #pragma omp parallel for default(none) \
         shared(NUFc, NCel, L0r1, grd, gsd, wrk0, wrk1) \
         private(cnst1, cnst2, cnst5)
    for (int i=0; i<NUFc; i++) {
        sidtyp sidu = gsd.ISD[i];
        size_t L = static_cast<size_t>( sidu[4] ); 
        //ISD[3:6] already converted to C++ index from 0 to NCel+10;
        size_t M = static_cast<size_t>( sidu[5] ); 

        // Only calculate boundary gradient if L0r1 = true.
        if( L0r1 || (L < NCel && M < NCel) ) {
            celtyp cell = grd.CEL[L];
            celtyp celm = grd.CEL[M];
            cnst1 = sidu[2];
            cnst2 = cell[2] + celm[2];
            cnst5 = cnst1*(wrk0[M] - wrk0[L])/cnst2;
            // Store flux into inner cells only.
            if( L < NCel ) {
                // Atomically update the element.
                #pragma omp atomic update
                wrk1[L] += cnst5;
            }
            if( M < NCel ) {
                // Atomically update the element.
                #pragma omp atomic update
                wrk1[M] += cnst5;
            }
        }  // End of if( L0r1 etc.
    }   // End of i-loop

    // Calculate cell centre x-gradient except for polar cells.
    #pragma omp parallel for default(none) \
         shared(NCel, NPol, Gdnt, grd, wrk1, DX0I) \
         private(cnst)
    for (int n=0; n<NCel-NPol; n++) {
        celtyp cell = grd.CEL[n];
        cnst = grd.CLat[n];
        Gdnt.dvdx[n] = DX0I*wrk1[n]/(cnst*cell[3]);
    }

    #pragma omp parallel for default(none) \
         shared(NVFc, NCel, L0r1, grd, gsd, wrk0, wrk2) \
         private(cnst3, cnst4, cnst6)
    for (int j=0; j<NVFc; j++) {
        sidtyp sidv = gsd.JSD[j];
        size_t L = static_cast<size_t>( sidv[4] ); 
        //ISD[3:6] already converted to C++ index from 0 to NCel+10;
        size_t M = static_cast<size_t>( sidv[5] ); 

        // Only calculate boundary gradient if L0r1 = true.
        if( L0r1 || (L < NCel && M < NCel) ) {
            celtyp cell = grd.CEL[L];
            celtyp celm = grd.CEL[M];
            cnst3 = sidv[2];
            cnst4 = cell[3] + celm[3];
            cnst6 = cnst3*(wrk0[M] - wrk0[L])/cnst4;
            // Store flux into inner cells only.
            if( L < NCel ) {
                // Atomically update the element.
                #pragma omp atomic update
                wrk2[L] += cnst6;
            }
            if( M < NCel ) {
                // Atomically update the element.
                #pragma omp atomic update
                wrk2[M] += cnst6;
            }
        }   // End of if( L0r1 etc ).
    }  // End of j-loop

    // Calculate cell-centre y-gradient except for polar cells.
    #pragma omp parallel for default(none) \
         shared(NCel, NPol, Gdnt, grd, wrk2, DY0I) 
    for (int k=0; k<NCel-NPol; k++) {
        celtyp cell = grd.CEL[k];
        // Cell x-size is used to cancell y-face factors in wrk2.
        Gdnt.dvdy[k] = DY0I*wrk2[k]/(float)cell[2];
    }  // End of k-loop

    // Set polar cell gradient to be zero, if any.
    if( NPol > 0 ) {
        for (int m=0; m<NPol; m++) {
            Gdnt.dvdx[NCel-NPol+m] = 0.0;
            Gdnt.dvdy[NCel-NPol+m] = 0.0;
        }
    }

    return(0); 

    // End of SMCGradn function.
}

