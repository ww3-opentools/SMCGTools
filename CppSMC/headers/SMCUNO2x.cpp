/* The SMCUNO2x module calculates mid-fluxes on u-side faces on a SMC grid. 
 * It works on portions of u-faces for multi-resolution use.  JGLi09Jun2026 */

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
int SMCUNO2x( int Lvl, Propxy& Pro, bool L0r1 ) {
    
    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    const auto& gsd = getSMCSID();   // Get gsd as read only const reference.
    auto NCel = grd.ncel;
    auto NPol = grd.npol;
//  cout << " Total number of cells grd.ncel = " << NCel << endl;

    //  Check iuf < nuf for u-face loop bounds.
    auto iuf = gsd.isdlv[Lvl -1];
    auto nuf = gsd.isdlv[Lvl];
//  cout <<" Lvl iuf nuf = "<< Lvl <<" "<< iuf <<" "<< nuf << endl;
    if( iuf > nuf ) {
        cout << " * U-side count iuf = " << iuf << endl;
        cout << " is larger than nuf = " << nuf << endl;
        exit(1);
    }

    // Confirm boundary condition to be first (true) or second (false) type.
//  if( L0r1 ) cout << " First type boundary with zero boundary cells. \n";
//  else  cout << " Second type boundary condition with zero-gradient. \n";
 
    // Set fraction timestep step for given level.
    int FTS = 1 << (Lvl - 1);
 
    // Setup diffusion coefficient with fractional timestep.
    float ADFr = Pro.ADfu*FTS*2.0;

    // Assign 0.0 to working variable of dimension NCel + 10 boundary cells.
    // Note assign first resize the vector if not present or overwrite it.
    auto NC10 = NCel + 10;
    auto NUFc = gsd.nufc; 
    wrk0.assign(NC10, 0.0f); 
    wrk1.assign(NUFc, 0.0f); 
    wrk2.assign(NUFc, 0.0f); 
    
    // Copy Pro.cv to wrk1[0:NCel] and set 10 boundary values zero.
    std::copy(Pro.cv.begin(), Pro.cv.begin() + NCel, wrk0.begin());
//  cout << " Pro.cv passed to wrk0[NCel] = " << wrk0[NCel] <<"\n";

    #pragma omp parallel for default(none) \
         shared(iuf, nuf, Pro, grd, gsd, L0r1, NCel, NPol, FTS, ADFr, wrk0, wrk1, wrk2) \
         private(cnst1, cnst2, cnst3, cnst4, cnst5, cnst6)  
    for (int i=iuf; i<nuf; i++) {
        sidtyp sdu = gsd.ISD[i];
        size_t K = static_cast<size_t>( sdu[3] );
        size_t L = static_cast<size_t>( sdu[4] ); 
        size_t M = static_cast<size_t>( sdu[5] ); 
        size_t N = static_cast<size_t>( sdu[6] ); 

        // Main timestep Courier number multiplied with FTS.
        // Use bool L0r1 to decide face or central Courant number.
        if( L0r1 ) cnst6 = Pro.uu[i]*FTS;
        else {
            cnst6 = 0.5*(Pro.uc[L] + Pro.uc[M])*FTS;
            if( L >= NCel-NPol ) cnst6 = Pro.uc[M]*FTS;
            if( M >= NCel-NPol ) cnst6 = Pro.uc[L]*FTS;
        }

        // Calculate face fluxes includeing boundary faces.
        celtyp cell = grd.CEL[L];
        celtyp celm = grd.CEL[M];
        cnst5 = (wrk0[M] - wrk0[L])/(cell[2]+celm[2]); 
            
        // For positive velocity case.
        if( cnst6 >= 0.0 ) {
            celtyp celu = grd.CEL[K];
            cnst4 = (wrk0[L] - wrk0[K])/(cell[2] + celu[2]); 
            // Minimum gradient with cnst5 sigin.
            cnst2 = copysign( min( abs(cnst4), abs(cnst5) ), cnst5 );
            wrk1[i] = (wrk0[L] + cnst2*(cell[2] - cnst6))*sdu[2];
        }
        // For negative velocity case.
        else { 
            celtyp celu = grd.CEL[N];
            cnst4 = (wrk0[N] - wrk0[M])/(celu[2] + celm[2]);
            // Minimum gradient with cnst5 sigin.
            cnst2 = copysign( min( abs(cnst4), abs(cnst5) ), cnst5 );
            wrk1[i] = (wrk0[M] - cnst2*(celm[2] + cnst6))*sdu[2];
        }

        // Use CLat of inner cell if L or M is boundary cell.
        // Diffusion flux always use L-M gradient.
        auto ij = min(L, M);
        cnst1 = 1.0/(grd.CLat[ij]);
        wrk2[i] = ADFr*cnst1*cnst1*cnst5*sdu[2];

        // Save fluxes into net-flux array, inner cells only.
        cnst3 = wrk1[i]*cnst6 - wrk2[i];
        if( L < NCel ) {
            // Atomically update the element.
            #pragma omp atomic update
            Pro.fx[L] -= cnst3; 
            #pragma omp atomic update
            Pro.fxad[L] -= (wrk1[i]*Pro.uc[L]*FTS - wrk2[i]); 
        }
        if( M < NCel ) {
            // Atomically update the element.
            #pragma omp atomic update
            Pro.fx[M] += cnst3; 
            #pragma omp atomic update
            Pro.fxad[M] += (wrk1[i]*Pro.uc[M]*FTS - wrk2[i]);
        }
        
    }   // End of i-loop

    // Update cv, cd for x-fluxes. cd for y-propagation, cv for final update.
    auto inc = grd.ncelv[Lvl -1];
    auto nnc = grd.ncelv[Lvl];
//  cout <<" Lvl inc nnc = "<< Lvl <<" "<< inc <<" "<< nnc << endl;

    #pragma omp parallel for default(none) \
         shared(inc, nnc, Pro, grd) 
    for (int n=inc; n<nnc; n++) {
        Pro.cd[n] = Pro.cv[n] + Pro.fx[n]*grd.CAra[n];
        Pro.cv[n] += Pro.fxad[n]*grd.CAra[n];
        // Reset Pro.fx and fxad[n] to zero.
        Pro.fx[n] = 0.0;
        Pro.fxad[n] = 0.0;
    }
//  cout << " Pro.cd updated Pro.cd[0] = " << Pro.cd[0] <<"\n";

    return(0); 

}   // End of SMCUNO2x function.

