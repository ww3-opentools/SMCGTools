/* The ArctBndy module updates the linking boundary conditions between the  *
 * global and Arctic parts with full Spc class.              JGLi08Jul2026  */

#include <iostream>
using namespace std;
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <algorithm>

#include "Arctics.h"
#include "WvSpectr.h"
#include "SMConsts.h"
using namespace SMConsts;
#include "PrintVector.h"

// Function to initialize Great circle turning parameters.
int ArctBndy( WvSpectr& Spc, const Arctics& Arc, int NT, bool TstOut ) {

    // Fetch Arctic part related parameters. 
    auto NCel = Arc.NCel;
    auto NGlo = Arc.NGlo;
    auto NArc = Arc.NArc;
    auto NGlB = Arc.NGlB;
    auto NArB = Arc.NArB;
    auto NPol = Arc.NPol;

    auto NSea = Spc.NSeapt();
    auto NFrq = Spc.NFrqcy();
    auto NDrc = Spc.NDrctn();
    auto MSpc = Spc.MSpctr();

    if( NSea < NCel ) {
        cout <<" *** Spc NSea is smaller than NCel. *** \n";
        exit(1);
    }

    // Global boundary condition from Arctic part cells.
    for(int i=0; i<NGlB; i++) {
        int j= i + NGlo - NGlB;
        int k= Arc.MBGlo[i];
        if( k - NGlo < 0 ) {
            cout <<" *** Global link bdy i = "<< i <<" no Arctic match! \n";
            continue;
        }
        // Otherwise complete this boundary point spectral update.
        auto spcta = Spc[k];
        auto Angla = Arc.Angld[k-NGlo];

        if( TstOut and (i+1) % 12 == 0 and i < 40 ) {
            string Flnm = "Sa"+std::to_string(i)+"i";
            string Msgs = " Arctic Spc["+std::to_string(k)+
                  "] before turn for global i = "+ std::to_string(i); 
            PrintVector(spcta, MSpc, NT, Flnm, Msgs);
        } // End test output.
 
        // Turn the Arctic spectrum at k by Angla.
        Specturn( spcta, NFrq, NDrc, Angla );
        // Copy the rotated spectrum to global boundary cell.
        Spc[j] = spcta;

        if( TstOut and (i+1) % 12 == 0 and i < 40 ) {
            string Flnm = "Sa"+std::to_string(i)+"t";
            string Msgs = " Global Spc["+std::to_string(j)+
                  "] after turn Angla = "+ std::to_string(Angla);
            PrintVector(spcta, MSpc, NT, Flnm, Msgs);
        } // End test output.
 
    } // End of all global boundary point i loop.

    // Arctic boundary condition from global part cells.
    for(int l=0; l<NArB; l++) {
        int m= l + NGlo;
        int n= Arc.MBArc[l];
        if( n < 0 ) {
            cout <<" *** Arctic link bdy l = "<< l <<" no Glboal match! \n";
            continue;
        }
        // Otherwise complete this boundary point spectral update.
        auto spctg = Spc[n];
        auto Anglg = -Arc.Angld[l];

        if( TstOut and (l+1) % 12 == 0 and l < 40 ) {
            string Flnm = "Sg"+std::to_string(l)+"i";
            string Msgs = " Global Spc["+std::to_string(n)+
                  "] before turn for Arctic l = "+ std::to_string(l); 
            PrintVector(spctg, MSpc, NT, Flnm, Msgs);
        } // End test output.
 
        // Turn the Arctic spectrum at n by -Angld.
        Specturn( spctg, NFrq, NDrc, Anglg );
        // Copy the rotated spectrum to global boundary cell.
        Spc[m] = spctg;

        if( TstOut and (l+1) % 12 == 0 and l < 40 ) {
            string Flnm = "Sg"+std::to_string(l)+"t";
            string Msgs = " Arctic Spc["+std::to_string(m)+
                  "] after turn Anglg = "+ std::to_string(Anglg);
            PrintVector(spctg, MSpc, NT, Flnm, Msgs);
        } // End test output.

    } // End of all Arctic boundary point l loop.

    return(0);
  
} // End of function ArctBndy.cpp module.


