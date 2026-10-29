/****************************************************************************
 *  This readcell.cpp program reads cell arrays from a cell array text file.
 *  Use vector<celtyp> to define cell arrays.      JGLi26Apr2026  
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

namespace {     // Anonymous namespace only available to this file.
    struct SMCELL grd;
    bool initialized = false;
}

/* Define a format to input celtyp variable.  */
istream &operator >>(istream &stream, celtyp &one) {
    for( int i=0; i<5; i++) stream >> one[i];
    return stream;
}

/* Define a format to output celtyp variable. */
ostream &operator <<(ostream &stream, celtyp &one) {
    stream << setw(6) << one[0];
    stream << setw(6) << one[1];
    stream << setw(5) << one[2];
    stream << setw(4) << one[3];
    stream << setw(6) << one[4] << '\n';
    return stream;
}

/* A function to read cell array and fill in SMCELL structure. */
void readSMCELL(const std::string& inpmodel, int nlv1, int NCM, int NPol, 
                float zdlnlt[4], bool Arctic) {  
    if( initialized ) return;

    /* Store zdlnlt and NPol parameters to struct grd.  */
    grd.nlvl = nlv1 - 1;
    grd.npol = NPol;
    grd.zlon = zdlnlt[0];
    grd.zlat = zdlnlt[1];
    grd.dlon = zdlnlt[2];
    grd.dlat = zdlnlt[3];

    /* Reset Arctic = false if NPol = 0  */
    if( NPol == 0 ) {
        Arctic = false;
        cout << " No Arctic part as NPol = " << NPol << endl;
    }
    grd.Arctic = Arctic;

    /* Define a cell structure and resize some arrays.  */
    celtyp cell;
    grd.ncelv.resize(nlv1);
    grd.CEL.resize(NCM);

    /* Declare an input file stream from the infile. */
    std::string infile, line;
    infile = inpmodel + "Cels.dat";
    ifstream celin(infile);
    if( !celin ) cout << " *** "  << infile << " is not opened!" << endl;
    else cout << "Reading cell array from " << infile << endl;

    for( int n=0; n<nlv1; n++ ) {
        celin >> grd.ncelv[n];   // Read one count line integer.  
        cout << setw(8) << grd.ncelv[n];   // Print it out for check.
    }
    std::getline(celin, line);  // Get the rest of count line.
    grd.ncel = grd.ncelv[0];     // Total number as first element of vector ncelv.
    cout << "\n ncel = " << grd.ncel << endl;

    // Assign the cell array from input file.
    int ng;  // To mark the end of global part cell array.
    for( ng=0; ng<grd.ncel; ng++ ) {
        celin >> grd.CEL[ng];
    }
    grd.nglo = ng;

    if( grd.nglo != grd.ncel ) {
        cout << " *** Warning ncel not matching nglo! \n";
        cout << " nglo = " << grd.nglo << '\n';
    }
   
    cout << "First cell is "; 
    cout << grd.CEL[0];
    cout << " Last cell is "; 
    cout << grd.CEL[grd.ncel-1];

    celin.close();   // Close the input file.

    if( grd.Arctic ) {   // Append Arctic cell array to global one.
        /* Declare an input file stream from the infile. */
        grd.narcs.resize(3);
        infile = inpmodel + "BArc.dat";
        ifstream celin(infile);
        if( !celin ) cout << " *** "  << infile << " is not opened!" << endl;
        else cout << "Appending Arctic part from " << infile << endl;

        for ( int n=0; n<3; n++ ) {
            celin >> grd.narcs[n];  
            cout << setw(8) << grd.narcs[n];   // Print it out for check.
        }
        cout << endl;
        std::getline(celin, line);  // Get the rest of count line.

        // Add Arctic cell numbers to total number and base level number.
        grd.ncel += grd.narcs[0];	
        grd.ncelv[nlv1 -1] += grd.narcs[0]; 
     	cout << " Merged ncel = " << grd.ncel << endl;

        // Append Arctic cell array to global one.
	for( int n=0; n<grd.narcs[0]; n++ ) {
            celin >> grd.CEL[grd.nglo + n];
        } 
        cout << "First cell is "; 
        cout << grd.CEL[grd.nglo];
        cout << " Last cell is "; 
        cout << grd.CEL[grd.ncel-1];

        celin.close();   // Close the input file.
    }

    /* Append boundary cells to end of CEL.  JGLi04Jun2026 */
    int NCel = grd.ncel;  // Boundary cell starting index.
    int mm = 1;           // Start with size-1 cell.
    int MJ = 1 << (nlv1-2);   // bitshift of 2**(nlv1-2); 
   
    cell[0] = 0;
    cell[1] = 0;
    cell[4] = 10; 
    for( int n=0; n<10; n++ ) {
        cell[2] = mm; 
        cell[3] = min(mm, MJ); 
        grd.CEL[NCel+n] = cell;
        mm *= 2;
    }
    cout << "Maximum boundary cell j-side = " << MJ << endl;

    /* Calculate cosine of cell centre latitude and cell area. */
    grd.CLat.resize(NCM);
    grd.CAra.resize(NCM);
    float cnst, cnst0, dlathf = 0.5*grd.dlat;
    for( int n=0; n<NCel-NPol; n++ ) {
        cell = grd.CEL[n];
        cnst = grd.zlat + (2*cell[1] + cell[3])*dlathf;
        grd.CLat[n] = cos( cnst*D2Rad );
        grd.CAra[n] = 1.0/(float)(cell[2]*cell[3]);
    }
    /* Polar cell CLat is at cell edge and CAra is a circle one. */
    int n = NCel-NPol;
    if( NPol > 0 ) {
        cell = grd.CEL[n];
//  Polar cell is at 90.0 and Clat uses positive cosine of side latitude. 
        cnst = grd.zlat + cell[1]*grd.dlat;
        grd.CLat[n] = abs( cos( cnst*D2Rad ) );
        if( cnst > 0.0 ) cnst0 = 90.0 - cnst;
        else  cnst0 = -90.0 - cnst;
        grd.PCRd = cnst0;
        grd.CAra[n] = grd.dlon*grd.dlat*grd.CLat[n]/(Pie*cnst0*cnst0); 
        cout << " Polar cell radius in lat deg = " << grd.PCRd << endl;
        cout << " Polar cell area in radian**2 = " << grd.CAra[n] << endl;
        cout << " Polar cell side cos(latside) = " << grd.CLat[n] << endl;
    }
    if( NPol > 1 ) {  // Suppose two polar cells are of equal area.
        grd.CLat[n+1] = grd.CLat[n];
        grd.CAra[n+1] = grd.CAra[n];
    }        
    
    // Change ncelv into level loop counts.  JGLi16Jun2026
    grd.ncelv[0] = 0;
    for( int n=1; n<nlv1; n++ ) {
         grd.ncelv[n] += grd.ncelv[n-1];
    }
    cout <<" Level loop counts set for NLvl + 1 = " << nlv1 << endl; 
    for(int n=0; n<nlv1; n++ ) cout <<" "<<grd.ncelv[n];
    cout << endl;
 
    initialized = true;

}

// End of readSMCELL function.

/* A function to fetch struct SMCELL for shared use. */
const SMCELL& getSMCELL() {
    if (!initialized) {
        throw std::runtime_error("*** SMCELL data not initialized.");
    }
    // Return the SMCELL grd data, otherwise.
    return grd;
}

// End of readSMCELL.cpp module.

