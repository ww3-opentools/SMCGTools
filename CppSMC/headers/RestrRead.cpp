/***************************************************************************
 *  The RestrRead.cpp program read the full Spc from a restart file written 
 *  by RestWrite.cpp program from previous run.        JGLi19Aug2026  
 ***************************************************************************/
#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <cstdio>    // for C-style prontf format.
#include <string>
#include <vector>

#include "Restrfile.h"

void RestrRead( WvSpectr& Spc, int& NSt, string RstrFile) {
    ifstream finp(RstrFile);
    if( !finp ) {
        cout <<RstrFile+" restart file not opened! \n";
        return;
    }

    int NSe, MSp;
    std::string line;
    auto NSea = Spc.NSeapt();
    auto MSpc = Spc.MSpctr();

    finp >> NSt;
    cout <<" *** Restart time NSt = "<< NSt <<endl;
    finp >> NSe;
    if( NSe != NSea ) {
        cout <<" *** Restart NSe not matching NSea: "<< NSe <<" "<< NSea<<endl;
        return;
    }
    finp >> MSp;
    if( MSp != MSpc ) {
        cout <<" *** Restart MSp not matching MSpc: "<< MSp <<" "<< MSpc<<endl;
        return;
    }
    std::getline(finp, line);  // Get the rest of count line, if any.

    // Start reading save full spectral file.
    for(int i = 0; i < NSea; ++i)
    {
        for(int k = 0; k < MSpc; ++k) finp >> Spc(i, k);
    }
    finp.close();

    cout <<" *** Restart file read from "<< RstrFile <<" done. \n";
    return;
}
// End of RestrRead function. 

