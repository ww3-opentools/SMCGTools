/* The SpctrIni module initializes the wave spectral field and stores it in 
 * the class WvSpectr storage space.                          JGLi09Jun2026 */

#include <omp.h>
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "SMCGrid.h"
#include "Working.h"

#include "WvSpectr.h"
#include "SMConsts.h"
using namespace SMConsts;

// Function to initialize SMC grid cell wave spectral field. 
int Specturn( std::vector<float>& spct, std::size_t NFrq, std::size_t NDrc, 
              float AnglD ) {
 
    if( NDrc < 1 || NFrq < 1 ) {
        cout <<" *** NDrc/NFrq < 1, NDrc = "<< NDrc <<"\n";
        exit(1);
    }

    // Rotation angle divided by directional bin width in degree. 
    cnst = AnglD*NDrc/360.0f; 
    int nsft = static_cast<int>(cnst);
    int ndrc = static_cast<int>(NDrc);
    float frac = cnst - nsft;
    // Restore the original wrapping line.  It caused the Arctic blocking!
    // This line seems no effect at all but it really did something.  Because 
    // NDrc is a size_t integer, nsft is first converted into size_t 64-bit 
    // long positive integer before the modulo operation.  This will be fine 
    // for positve integer nsft but for a negative nsft, the result will be a 
    // integer >= 0 as negative integer converted to size_t long integer is 
    // the complimentary value from the end of the maximum long integer. This 
    // will be wrong for this scheme.  JGLi18Aug2026
//  nsft %= NDrc;
    // This above line is replaced with the following line if required.
    nsft %= ndrc;
    // As normally abs(nsft) < NDrc, this line is not required unless the 
    // abs(AnglD) > 360.0f becomes true. 

    // Loop over all frequency range.
    #pragma omp parallel for default(none) \
         shared(NFrq, NDrc, spct, nsft, frac, ndrc) \
         private(wrk1, wrk2, wrk3, cnst1, cnst2)
    for(int n=0; n<NFrq; n++) {
        // Assign 0.0 to working variable of required dimensions.
        // Note assign first resize the vector if not present or overwrite it.
        wrk1.assign(NDrc, 0.0f); 
        wrk2.assign(NDrc, 0.0f); 
        wrk3.assign(NDrc, 0.0f); 

        // Extract one frequency spectrum to wrk1.
        std::copy(spct.begin()+NDrc*n, 
                  spct.begin()+NDrc*(n+1), wrk1.begin());

        // Shift nsft bins if abs(nsft) > 0.
        if( abs( nsft ) > 0 ) {
            for(int k=0; k<ndrc; k++) {
                int jj = k + nsft;
                if( jj < 0 )     jj += ndrc;
                if( jj >= ndrc ) jj -= ndrc;
                wrk2[k] = wrk1[jj];
            }
        }
        else { // Copy wrk1 to wrk2.
            wrk2 = wrk1;
//          std::copy(wrk1.begin(), wrk1.begin()+NDrc, wrk2.begin());
        }
        
        // Fractional partitioning clockwise if abs(frac) > 0.0. 
        if( frac > 0.0f ) {
            cnst1 = wrk2[0]*frac; 
            for(int k=NDrc-1; k>=0; k--) {
                cnst2 = wrk2[k]*frac;
                wrk3[k] = wrk2[k] + cnst1 - cnst2;
                cnst1 = cnst2;
            }
        }
        else {  // Anti-clockwise passing -frac part.
            cnst1 = wrk2[NDrc-1]*frac; 
            for(int k=0; k<NDrc; k++) {
                cnst2 = wrk2[k]*frac;
                wrk3[k] = wrk2[k] - cnst1 + cnst2;
                cnst1 = cnst2;
            }
        }
        // Copy wrk3 back to spct. Note spct is a reference input, so the 
        // update is automatically stored into spct in the calling scope.
        std::copy( wrk3.begin(), wrk3.begin()+NDrc, spct.begin()+NDrc*n );

    } // End of frequency n loop.

    return(0); 

}  // End of Specturn function.


