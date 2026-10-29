/* The SMCUNO2y module calculates mid-fluxes on v-side faces on a SMC grid. 
 * It works on portions of v-faces for multi-resolution use.  JGLi10Jun2026 */

#include <omp.h>
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "SMCGrid.h"
#include "Arctics.h"
#include "Working.h"
#include "Rotated.h"
#include "Propgtn.h"

#include "SMConsts.h"
using namespace SMConsts;

// Function to calculate variable gradient on a SMC grid. 
int SMCUNO2y( int Lvl, Propxy& Pro, bool L0r1 ) {
    
    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    const auto& gsd = getSMCSID();   // Get gsd as read only const reference.
    auto NCel = grd.ncel;
    auto NPol = grd.npol;
//  cout << " Total number of cells grd.ncel = " << NCel << endl;

    //  Check jvf < nvf for u-face loop bounds.
    auto jvf = gsd.jsdlv[Lvl - 1];
    auto nvf = gsd.jsdlv[Lvl];
//  cout <<" Lvl jvf nvf = "<< Lvl <<" "<< jvf <<" "<< nvf << endl;
    if( jvf > nvf ) {
        cout << " * V-side count jvf = " << jvf << endl;
        cout << " is larger than nvf = " << nvf << endl;
        exit(1);
    }

// Bug searching check of Pro.fy values at size 4 cell ijk=52600.
//  int ijk = 52600;

// Output Pro.fy value before Lvl = 3 calculation.
//  if( Lvl == 3 ) {
//  cout <<" Before Lvl = "<< Lvl <<" Pro.fy[ijk] = "<< Pro.fy[ijk] << endl;
//  }

    // Confirm boundary condition to be first (true) or second (false) type.
//  if( L0r1 ) cout << " First type boundary with zero boundary cells. \n";
//  else  cout << " Second type boundary condition with zero-gradient. \n";
 
    // Set fraction timestep step for given level.
    int FTS = 1 << (Lvl - 1);

    // Setup diffusion coefficient with fractional timestep.
    // Factor 2.0 is added to cancell the doubled distance in gradient.
    float ADFr = Pro.ADfv*FTS*2.0;

    // Assign 0.0 to working variable of dimension NCel + 10 boundary cells.
    // Note assign first resize the vector if not present or overwrite it.
    auto NC10 = NCel + 10;
    auto NVFc = gsd.nvfc; 
    wrk0.assign(NC10, 0.0f); 
    wrk1.assign(NVFc, 0.0f); 
    wrk2.assign(NVFc, 0.0f); 
    
    // Copy Pro.cd to wrk1[0:NCel] and set 10 boundary values zero.
    std::copy(Pro.cv.begin(), Pro.cv.begin() + NCel, wrk0.begin());
//  cout << " Pro.cd passed to wrk0[NCel] = "<< wrk0[NCel] <<"\n";

    #pragma omp parallel for default(none) \
         shared(jvf, nvf, Pro, grd, gsd, L0r1, NCel, NPol, FTS, ADFr, wrk0, wrk1, wrk2) \
         private(cnst1, cnst2, cnst3, cnst4, cnst5, cnst6)
    for (int j=jvf; j<nvf; j++) {
        sidtyp sdv = gsd.JSD[j];
        size_t K = static_cast<size_t>( sdv[3] );
        size_t L = static_cast<size_t>( sdv[4] ); 
        size_t M = static_cast<size_t>( sdv[5] ); 
        size_t N = static_cast<size_t>( sdv[6] ); 

        // Main timestep Courier number multiplied with FTS.
        // Exclude Polar and boundary cell average, using inner cells.
        if( L0r1 ) cnst6 = Pro.vv[j]*FTS;
        else {
            cnst6 = 0.5*(Pro.vc[L] + Pro.vc[M])*FTS;
            if( L >= NCel-NPol ) cnst6 = Pro.vc[M]*FTS;
            if( M >= NCel-NPol ) cnst6 = Pro.vc[L]*FTS;
        }

        // JSD size with cosine of latitude factor.
        cnst1 = gsd.CLatv[j]*sdv[2];

        // Calculate face fluxes including boundary faces. 
        celtyp cell = grd.CEL[L];
        celtyp celm = grd.CEL[M];
        cnst5 = (wrk0[M] - wrk0[L])/(cell[3]+celm[3]); 
            
        // For positive velocity case.
        if( cnst6 >= 0.0 ) {
            celtyp celu = grd.CEL[K];
            cnst4 = (wrk0[L] - wrk0[K])/(cell[3] + celu[3]); 
            // Minimum gradient with cnst5 sigin.
            cnst2 = copysign( min( abs(cnst4), abs(cnst5) ), cnst5 );
            wrk1[j] = (wrk0[L] + cnst2*(cell[3] - cnst6))*cnst1;
        }
        // For negative velocity case.
        else { 
            celtyp celu = grd.CEL[N];
            cnst4 = (wrk0[N] - wrk0[M])/(celu[3] + celm[3]);
            // Minimum gradient with cnst5 sigin.
            cnst2 = copysign( min( abs(cnst4), abs(cnst5) ), cnst5 );
            wrk1[j] = (wrk0[M] - cnst2*(celm[3] + cnst6))*cnst1;
        }

        // Diffusion flux always use L-M gradient.
        wrk2[j] = ADFr*cnst5*cnst1;

        // Save y-fluxes into net-flux arrays.
        cnst3 = wrk1[j]*cnst6 - wrk2[j];
        if( L < NCel ) {
            // Atomically update the element.
            #pragma omp atomic update
            Pro.fy[L] -= cnst3; 
        }
        if( M < NCel ) {
            // Atomically update the element.
            #pragma omp atomic update
            Pro.fy[M] += cnst3; 
        }

    }   // End of j-loop

// Output Pro.fy value after Lvl = 2 calculation.
//  if( Lvl == 2 ) {
//  cout <<" After  Lvl = "<< Lvl <<" Pro.fy[ijk] = "<< Pro.fy[ijk] << endl;
//  }

    // Update cv for y-fluxes. Note cv was copied from Prx.cv. 
    auto jnc = grd.ncelv[Lvl - 1];
    auto nnc = grd.ncelv[Lvl];
//  cout <<" Lvl jnc nnc = "<< Lvl <<" "<< jnc <<" "<< nnc << endl;

    #pragma omp parallel for default(none) \
         shared(jnc, nnc, Pro, grd) 
    for (int n=jnc; n<nnc; n++) {
        Pro.cv[n] = Pro.cd[n] + (Pro.fy[n]*grd.CAra[n]/grd.CLat[n]);
        // Reset Pro.fy[n] to zero.
        Pro.fy[n] = 0.0;
    }

    return(0); 

}   // End of SMCUNO2y function.

