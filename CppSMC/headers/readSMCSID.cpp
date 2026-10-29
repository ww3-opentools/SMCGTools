/****************************************************************************
 *  This readSMCSID.cpp program reads face arrays from face array text files.
 *  Use vector<celtyp> to define cell arrays.      JGLi26May2026  
 ****************************************************************************/

#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <stdexcept>
#include <cmath>

#include "SMConsts.h"
using namespace SMConsts;

#include "SMCGrid.h"

/****************************************************************************
 *  This readside.cpp program reads side arrays from side array text files.
 *  Use vector<sidtyp> to define side arrays.      JGLi26Apr2026  
 ****************************************************************************/

namespace {
    struct SMCSID gsd;
    bool initialside = false;
}

/* SMC grid side array has 7/8 elements but the J-side last element is only 
 * used for sorting the face array and it is not required in propagation subs. 
 * So sidtyp definition will only use 7 elements for both I and J-side arrays. 
 * J-side last element will be skipped at file reading stage. JGLi26Mar2026 */

/* Define a format to input sidtyp variable.  */
istream &operator >>(istream &stream, sidtyp &one) {
    for( int i=0; i<7; i++) stream >> one[i];
    return stream;
}

/* Define a format to output sidtyp variable. */
ostream &operator <<(ostream &stream, sidtyp &one) {
    stream << setw(6) << one[0];    // I-index S/W
    stream << setw(6) << one[1];    // J-index S/W
    stream << setw(5) << one[2];    // Size of side
    stream << setw(8) << one[3];    // mm cell index
    stream << setw(8) << one[4];    // m  cell index
    stream << setw(8) << one[5];    // n  cell index
    stream << setw(8) << one[6];    // nn cell index
    stream << '\n';
    return stream;
}

/* A function to read face arrays and fill in SMCSID structure. */
void readSMCSID(const std::string& inpmodel, int nlv1, int NFM, int NPol) { 
    if( initialside ) return;

    /* Use SMCELL grd data by loading it.  */
    const auto& grd = getSMCELL();
    auto NCel = grd.ncel;
    auto NGlo = grd.nglo;
    cout << " Total cell number grd.ncel = " << NCel << endl;
    cout << " Global part cell number NGlo = " << NGlo << endl;
    cout << " Arctic switch of grd.Arctic = " << grd.Arctic << endl;

    /* Resize i-side variables in sid struct.   */
    gsd.isdlv.resize(nlv1);
    gsd.ISD.resize(NFM);

    /* Declare an input file stream from the infile. */
    std::string infile, line;
    infile = inpmodel + "ISid.dat";
    ifstream sidin(infile);
    if( !sidin ) cout << " *** "  << infile << " is not opened!" << endl;
    else cout << "Read I-side " << infile << endl;

    for( int n=0; n<nlv1; n++ ) sidin >> gsd.isdlv[n];  
    std::getline(sidin, line);  // Get the rest of count line.
    cout << "nufc = " << gsd.isdlv[0] << endl;

    // Assign the side array from input file.
    int nig, m;  // To mark the end of global part cell array.
    for( nig=0; nig<gsd.isdlv[0]; nig++ ) {
        sidin >> gsd.ISD[nig];
    }  

    cout << "First i-side is "; 
    cout << gsd.ISD[0];
    cout << " Last i-side is "; 
    cout << gsd.ISD[nig-1];
    cout << " nig = " << nig <<" "<< gsd.isdlv[0] << endl;

    sidin.close();   // Close the input file.

    if( grd.Arctic ) {   // Append Arctic side array to global one.
        /* Declare an input file stream from the infile. */
        infile = inpmodel + "AISd.dat";
        ifstream sidin(infile);
        if( !sidin ) cout << " *** "  << infile << " is not opened!" << endl;
        else cout << "Append Arctic " << infile << endl;

        int nia;    // Total Arctic I-side number.
        sidin >> nia;   // Only read the total Acrtic I-side number.  
        std::getline(sidin, line);  // Get the rest of count line.
     	cout << "nia  = " << nia  << endl;

        // Append Arctic side array to global one.
	for( int n=0; n<nia; n++ ) {
            sidin >> gsd.ISD[nig+n];
        }  
  
        // Add Arctic side counts to global ones.
        gsd.nuar = nia;
        gsd.isdlv[0] += nia; 
        gsd.isdlv[nlv1-1] += nia; 

        cout << "First Arctic i-side is "; 
        cout << gsd.ISD[nig];
        cout << " Last Arctic i-side is "; 
        cout << gsd.ISD[gsd.isdlv[0]-1];

        // Add NGlo to Arctic face array mm m n nn values.
	for( int n=0; n<nia; n++ ) {
            for( int i=3; i<7; i++ ) {
                if( gsd.ISD[nig+n][i] > 0 )  gsd.ISD[nig+n][i] += NGlo; 
            }
        }  
        cout << " Adjusted last i-side is "; 
        cout << gsd.ISD[gsd.isdlv[0]-1];

        sidin.close();   // Close the input file.
    }
    gsd.nufc = gsd.isdlv[0];     // Updated total ISD rows.
    cout << "Updated nufc = " << gsd.nufc << endl;
       
    // Change isdlv into level loop counts.  JGLi16Jun2026
    gsd.isdlv[0] = 0;
    for( int n=1; n<nlv1; n++ ) {
         gsd.isdlv[n] += gsd.isdlv[n-1];
    }
    cout <<" I-side loop counts set for NLvl + 1 = " << nlv1 << endl;
    for( int n=0; n<nlv1; n++ ) cout <<" " << gsd.isdlv[n];
    cout << endl;

    // Deduct -1 from ISD[3:6] to C++ index and change k=-9:0 index to NCel-k.
    sidtyp sidu;
    for(int n=0; n<gsd.nufc; n++ ) {
        sidu = gsd.ISD[n];
        for(int i=3; i<7; i++ ) {
            if( sidu[i] <= 0 ) gsd.ISD[n][i] = NCel - sidu[i]; 
            else               gsd.ISD[n][i] -= 1;
        }
    }  
    cout << " ISD[3:6] converted to C++ index. \n";

    /* Resize j-side variables in sid struct.   */
    gsd.jsdlv.resize(nlv1);
    gsd.JSD.resize(NFM);

    /* Declare J-side file stream from the infile. */
    infile = inpmodel + "JSid.dat";
    ifstream sjdin(infile);
    if( !sjdin ) cout << " *** "  << infile << " is not opened!" << endl;
    else cout << "Read J-side " << infile << endl;

    for( int n=0; n<nlv1; n++ ) sjdin >> gsd.jsdlv[n];  
    std::getline(sjdin, line);  // Get the rest of count line.
    cout << "nvfc = " << gsd.jsdlv[0] << endl;

    // Assign the side array from input file.
    int njg;  // To mark the end of global part cell array.
    for( njg=0; njg<gsd.jsdlv[0]; njg++ ) {
        sjdin >> gsd.JSD[njg];
        sjdin >> m;  // Skip last element of J-side array.
    }  

    cout << "First J-side is "; 
    cout << gsd.JSD[0];
    cout << " Last J-side is "; 
    cout << gsd.JSD[gsd.jsdlv[0]-1];
    cout << " njg = " << njg <<" "<< gsd.jsdlv[0] << endl;

    sjdin.close();   // Close the input file.

    if( grd.Arctic ) {   // Append Arctic side array to global one.
        /* Declare Arctic J-side file stream from the infile. */
        infile = inpmodel + "AJSd.dat";
        ifstream sjdin(infile);
        if( !sjdin ) cout << " *** "  << infile << " is not opened!" << endl;
        else cout << "Append Arctic " << infile << endl;

	int nja;     // Arctic J-side number.
        sjdin >> nja;    // Read only the total Arctic J-side number.  
        std::getline(sjdin, line);  // Get the rest of count line.
     	cout << "nja  = " << nja  << endl;

        // Append Arctic side array to global one.
	for(int n=0; n<nja; n++ ) {
            sjdin >> gsd.JSD[njg+n];
            sjdin >> m;  // Skip last element of J-side array.
        }  
  
        // Add Arctic side counts to global ones.
        gsd.nvar = nja; 
        gsd.jsdlv[     0] += nja; 
        gsd.jsdlv[nlv1-1] += nja;

        cout << "First Arctic J-side is "; 
        cout << gsd.JSD[njg];
        cout << " Last Arctic J-side is "; 
        cout << gsd.JSD[gsd.jsdlv[0]-1];

        // Add NGlo to Arctic face array mm m n nn values.
	for( int n=0; n<nja; n++ ) {
            for( int j=3; j<7; j++ ) {
                if( gsd.JSD[njg+n][j] > 0 )  gsd.JSD[njg+n][j] += NGlo; 
            }
        }  
        cout << " Adjusted last J-side is "; 
        cout << gsd.JSD[gsd.jsdlv[0]-1];
        sjdin.close();   // Close the input file.
    }
       
    gsd.nvfc = gsd.jsdlv[0];     // Updated total JSD rows.
    cout << "Updated nvfc = " << gsd.nvfc << endl;

    // Change jsdlv into level loop counts.  JGLi16Jun2026
    gsd.jsdlv[0] = 0;
    for( int n=1; n<nlv1; n++ ) {
         gsd.jsdlv[n] += gsd.jsdlv[n-1];
    }
    cout <<" J-side loop counts set for NLvl + 1 = " << nlv1 << endl;
    for(int n=0; n<nlv1; n++ ) cout <<" "<< gsd.jsdlv[n];
    cout << endl;

    // Deduct -1 from JSD[3:6] to C++ index and change k=-9:0 index to NCel-k.
    // Also calculate cosine of JSD latitude, CLatv. 
    sidtyp sidv;
    gsd.CLatv.resize(NFM);
    for( int n=0; n<gsd.nvfc; n++ ) {
         sidv = gsd.JSD[n];
         for( int i=3; i<7; i++ ) {
             if( sidv[i] <= 0 ) gsd.JSD[n][i] = NCel - sidv[i];
             else               gsd.JSD[n][i] -= 1;
         }
         gsd.CLatv[n] = cos( (grd.zlat + sidv[1]*grd.dlat)*D2Rad );
    }
    cout << " JSD[3:6] converted to C++ index. \n";

    initialside = true;

}

// End of readSMCSID function.
  
/* A function to fetch struct SMCSID for shared use. */
const SMCSID& getSMCSID() {
    if (!initialside) {
        throw std::runtime_error("*** SMCSID data not initialized.");
    }
    // Return the SMCSID gsd data, otherwise.
    return gsd;
}

// End of readSMCSID.cpp module.
 
