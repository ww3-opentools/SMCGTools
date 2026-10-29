/* The ArctBndy module updates the linking boundary conditions between the  *
 * global and Arctic parts with full Spc class.              JGLi08Jul2026  */

#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <omp.h>

#include "Arctics.h"
#include "SMConsts.h"
using namespace SMConsts;

// Function to update Arctic part linking boundary for scalar variable.
int ArcScalr( std::vector<float>& Sclr, const Arctics& Arc ) {

    // Fetch Arctic part related parameters. 
    auto NCel = Arc.NCel;
    auto NGlo = Arc.NGlo;
    auto NArc = Arc.NArc;
    auto NGlB = Arc.NGlB;
    auto NArB = Arc.NArB;
    auto NPol = Arc.NPol;

    // Check input Sclr size, not less than NCel.
    auto NSea = Sclr.size();

    if( NSea < NCel ) {
        std::cout <<" *** Sclr NSea is smaller than NCel. *** \n";
        exit(1);
    }

    // Global boundary condition from Arctic part cells.
    #pragma omp parallel for default(none) \
        shared(Sclr, Arc, NGlo, NGlB)
    for(int i=0; i<NGlB; i++) {
        int j= i + NGlo - NGlB;
        int k= Arc.MBGlo[i];
        Sclr[j] = Sclr[k];
    } // End of all global boundary point i loop.

    // Arctic boundary condition from global part cells.
    #pragma omp parallel for default(none) \
        shared(Sclr, Arc, NGlo, NArB)
    for(int l=0; l<NArB; l++) {
        int m= l + NGlo;
        int n= Arc.MBArc[l];
        Sclr[m] = Sclr[n];
    } // End of all Arctic boundary point l loop.

    return(0);
  
} // End of function ArcScalr. 

// Function to update Arctic part linking boundary for vector variable at cell centre.
int ArcVectr( std::vector<float>& Uctr, 
              std::vector<float>& Vctr, const Arctics& Arc ) {

    // Fetch Arctic part related parameters. 
    auto NCel = Arc.NCel;
    auto NGlo = Arc.NGlo;
    auto NArc = Arc.NArc;
    auto NGlB = Arc.NGlB;
    auto NArB = Arc.NArB;
    auto NPol = Arc.NPol;

    // Check input Sclr size, not less than NCel.
    auto NUct = Uctr.size();
    auto NVct = Uctr.size();

    if( NUct < NCel or NVct < NCel ) {
        std::cout <<" *** Uctr Vctr dimenions smaller than NCel. *** \n";
        exit(1);
    }

    // Global boundary condition from Arctic part cells.
    #pragma omp parallel for default(none) \
        shared(Uctr, Vctr, Arc, NGlo, NGlB)
    for(int i=0; i<NGlB; i++) {
        int j= i + NGlo - NGlB;
        int k= Arc.MBGlo[i];
        float Angd = D2Rad*Arc.Angld[k];
        Uctr[j] = Uctr[k]*cos(Angd) + Vctr[k]*sin(Angd);
        Vctr[j] = Vctr[k]*cos(Angd) - Uctr[k]*sin(Angd);
    } // End of all global boundary point i loop.

    // Arctic boundary condition from global part cells.
    #pragma omp parallel for default(none) \
        shared(Uctr, Vctr, Arc, NGlo, NArB)
    for(int l=0; l<NArB; l++) {
        int m= l + NGlo;
        int n= Arc.MBArc[l];
        float Angd = D2Rad*Arc.Angld[m];
        Uctr[m] = Uctr[n]*cos(Angd) - Vctr[n]*sin(Angd);
        Vctr[m] = Vctr[n]*cos(Angd) + Uctr[n]*sin(Angd);
    } // End of all Arctic boundary point l loop.

    return(0);
  
} // End of function ArcVectr. 

// End of ArcVSclr.cpp module.

