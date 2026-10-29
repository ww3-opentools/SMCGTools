/***************************************************************************
 *  The PrintVector.cpp program save a vector of NCel elements in a file of 
 *  name FN100NTS.dat and a count line NTS NCel.           JGLi16Jun2026  
 ***************************************************************************/
#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <cstdio>    // for C-style prontf format.
#include <string>
#include <vector>

#include "PrintVector.h"

void PrintVector(const std::vector<float>& Vctr, int NCel, int NTS, string FN, const string& Msgs ) {
    const int perRow = 10;
    std::string outfile = FN+std::to_string(100000+NTS)+".dat";
    ofstream fwrt(outfile);

    fwrt << setw(8) << NTS << setw(9) << NCel <<" "+Msgs << endl;
    for (int i = 0; i < NCel; ++i)
    {
        // Print value in scientific notation, width 12, precision 3
        fwrt << std::setw(12) << std::scientific << std::setprecision(3)
             << Vctr[i];

        // New line after every 10 elements or last element
        if ((i + 1) % perRow == 0 || i == NCel - 1)
            fwrt << '\n';
    }
    fwrt.close();
}
// End of PrintVector function. 

