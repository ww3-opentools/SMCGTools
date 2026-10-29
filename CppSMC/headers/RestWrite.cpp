/***************************************************************************
 *  The RestWrite.cpp program saves the full spectral class data Spc into a  
 *  text file with the prefix Flnm plus NTS stamp.          JGLi19Aug2026  
 ***************************************************************************/
#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <cstdio>    // for C-style prontf format.
#include <string>
#include <vector>

#include "Restrfile.h"

void RestWrite( WvSpectr& Spc, int NTS, string Flnm, string& Msgs) {
    const int perRow = 10;
    std::string outfile = Flnm+std::to_string(100000+NTS)+".dat";
    ofstream fwrt(outfile);

    auto NSea = Spc.NSeapt();
    auto MSpc = Spc.MSpctr();

    fwrt << setw(8) << NTS << setw(9) << NSea << setw(6) << MSpc <<" "+Msgs << endl;
    for(int i = 0; i < NSea; ++i)
    {
        for(int k = 0; k < MSpc; ++k) {
            // Print value in scientific notation, width 12, precision 3
            fwrt << std::setw(12) << std::scientific << std::setprecision(3)
                 << Spc(i, k);

            // New line after every 10 elements or last element
            if ((k + 1) % perRow == 0 || k == MSpc - 1) fwrt << '\n';
        }
    }
    fwrt.close();
}
// End of RestWrite function. 

