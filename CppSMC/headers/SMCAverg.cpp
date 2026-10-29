/* The SMCAverg module smoothes varaible on a SMC grid by 1-2-1 average over
 * both x and y-dimensions. It takes a vector variable reference input and 
 * replaces original variable with the averaged values.      JGLi12Jun2026  */ 

#include <omp.h>
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "SMCGrid.h"
#include "Working.h"

// Function to smooth variable on a SMC grid by 1-2-1 average. 
int SMCAverg( vector<float>& varb ) {
    //  Find out input variable size.
    auto NAvr = varb.size();
//  cout << " Input varabel varb.size() = " << NAvr << endl;

    // Fetch SMC grid cell and face arrays and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    const auto& gsd = getSMCSID();   // Get gsd as read only const reference.
    auto NCel = grd.ncel;
//  cout << " Number of SMC grid cells NCel = " << NCel << endl;
    if( NCel > NAvr ) {
        cout <<" *** "<< NCel <<" NCel larger than NAvr "<< NAvr << '\n'; 
        exit(1);
    }

    auto NPol = grd.npol;
    auto NUFc = gsd.nufc;
    auto NVFc = gsd.nvfc;
//  cout << NPol <<" polar cell value will be kept unchanged \n";

    // Assign 0.0 to working variable of dimension NCel + 10 boundary cells.
    // Note assign first resize the vector if not present or overwrite it.
    int NCB = NCel + 10;
    wrk0.assign(NCB, 0.0f); 
    wrk1.assign(NCB, 0.0f); 
    wrk2.assign(NCB, 0.0f); 
    wrk3.assign(NCB, 0.0f); 
    wrk4.assign(NCB, 0.0f); 
    
    // Copy vector varb[0:NCel] to wrk0[0:NCel]. 
    std::copy(varb.begin(), varb.begin() + NCel, wrk0.begin());

//  cout << " varb passed to wrk0[0:NCel] \n";
//  cout << " wrk0[0] = " << wrk0[0] <<"\n";
//  cout << " wrk0[NCel-1] = " << wrk0[NCel-1] <<"\n";
//  cout << " wrk0[NCB -1] = " << wrk0[NCB -1] <<"\n";

    #pragma omp parallel for default(none) \
         shared(gsd, NUFc, NCel, wrk0, wrk1, wrk3) \
         private(cnst1, cnst3 )
    for (int i=0; i<NUFc; i++) {
        sidtyp sidu = gsd.ISD[i];
        size_t L = static_cast<size_t>( sidu[4] ); 
        size_t M = static_cast<size_t>( sidu[5] ); 
        cnst3 = (float)sidu[2];

        // Only average with inner cells. 
        if( L < NCel && M < NCel ) {
            cnst1 = cnst3*(wrk0[M] + wrk0[L]);
            // Atomically update two cell variable.
            #pragma omp atomic update
            wrk1[L] += cnst1;
            #pragma omp atomic update
            wrk3[L] += cnst3;
            #pragma omp atomic update
            wrk1[M] += cnst1;
            #pragma omp atomic update
            wrk3[M] += cnst3;
        }  // End of if( L, M < NCel)

    }   // End of i-loop

    #pragma omp parallel for default(none) \
         shared(gsd, NVFc, NCel, wrk0, wrk2, wrk4) \
         private(cnst2, cnst4)
    for (int j=0; j<NVFc; j++) {
        sidtyp sidv = gsd.JSD[j];
        size_t L = static_cast<size_t>( sidv[4] ); 
        //ISD[3:6] already converted to C++ index from 0 to NCel+10;
        size_t M = static_cast<size_t>( sidv[5] ); 
        cnst4 = (float)sidv[2];

        // Only average with inner cells.
        if( L < NCel && M < NCel ) {
            cnst2 = cnst4*(wrk0[M] + wrk0[L]);
            // Atomically update two cell variable.
            #pragma omp atomic update
            wrk2[L] += cnst2;
            #pragma omp atomic update
            wrk4[L] += cnst4;
            #pragma omp atomic update
            wrk2[M] += cnst2;
            #pragma omp atomic update
            wrk4[M] += cnst4;
        }   // End of if( L0r1 etc ).

    }  // End of j-loop

    // Averaged cell-centre varb except for polar cells.
    #pragma omp parallel for default(none) \
         shared(NPol, NCel, wrk0, wrk1, wrk2, wrk3, wrk4) \
         private(cnst3, cnst4, cnst5, cnst6)
    for (int k=0; k<NCel-NPol; k++) {
        // Only update cell value when it is bounded with inner cells.
        cnst3 = 0.0;
        cnst4 = 0.0;
        cnst5 = 0.0;
        cnst6 = 0.0;
        if( wrk3[k] > 0.5f ) {
            cnst5 = wrk1[k]/wrk3[k];
            cnst3 = 2.0f;
        }
        if( wrk4[k] > 0.5f ) {
            cnst6 = wrk2[k]/wrk4[k];
            cnst4 = 2.0f;
        }
        if( abs(cnst5) + abs(cnst6) > 0.0 ) {
            wrk0[k] = (cnst5 + cnst6)/(cnst3 + cnst4);
        } 
    }  // End of k-loop

    // Copy averaged value back to varb[0:NCel]. 
    std::copy(wrk0.begin(), wrk0.begin() + NCel, varb.begin());

    return(0); 

    // End of SMCAverg function.
}

// End of module SMCAverg.cpp file.

