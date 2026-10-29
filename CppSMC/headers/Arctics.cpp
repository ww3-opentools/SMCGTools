/* The ArctAngd module contains functions to set up Arctics Arc sturct data 
 * and calculate the map-east angles to local east.         JGLi03Jun2026    */

#include <iostream>
using namespace std;
#include <fstream>
#include <string>
#include <vector>
#include <cmath>

#include "SMCGrid.h"
#include "Arctics.h"
#include "Rotated.h"
#include "GtCrTurn.h"

#include "SMConsts.h"
using namespace SMConsts;

namespace {     // Anonymous namespace only available to this file.
    struct Arctics Arc;
    bool initialized = false;
}

// Function to set up polar part links with global part. 
int ArctAngd( bool TstOut ) {

    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();

// Open a file to store Arctic settings for test check.
    ofstream ztest;
    if( TstOut ) {
        ztest.open("ArcSetting.txt");
        if( !ztest.is_open() ) cout << " *** Arctic output not opened!" << endl;
        ztest << " Arctic setting for Arc. " << '\n';
    }

//  Fetch a few grd parameters.
    auto NCel = grd.ncel;
    auto NGlo = grd.nglo;
    auto NPol = grd.npol;
    auto cone = static_cast<celtyp>(grd.CEL[NGlo]);
    auto ctwo = static_cast<celtyp>(grd.CEL[NCel-1]);
    cout << " Total cell number = " << NCel << endl;
    cout << " Total global cell = " << NGlo << endl;
    cout << " Number Polar cell = " << NPol << endl;
    cout << " First and last Arctic cell arrays = \n";
    cout << cone;  // Use celtyp its own output format.
    cout << ctwo;  // Use celtyp its own output format.

//  Fetch grd Arctic related parameters for Arc.
    auto NArc = grd.narcs[0];
    auto NArB = grd.narcs[1];
    auto NGlB = grd.narcs[2];

//  Resize Arc vectors and initialize them.
    Arc.Angld.assign(NArc, 0.0f);
    Arc.ELats.assign(NArc, 0.0f);
    Arc.MBArc.assign(NArB, -1);
    Arc.MBGlo.assign(NGlB, -1);
 
    Arc.NCel = NCel; 
    Arc.NGlo = NGlo; 
    Arc.NArc = NArc; 
    Arc.NArB = NArB; 
    Arc.NGlB = NGlB; 
    Arc.NPol = NPol;
     
//  if( TstOut ) {  // Same condition for opening the file stream.
    if( ztest ) {   // Equivalent to check the file stream is open.
        ztest << " Global linking boundary matchi: i, j, k, cell. \n";
    }

//  Matching global boundary cells with Arctic inner cells.
    int i, j, k, m, n=0; 
    for(i=0; i<NGlB; i++) {
        j = i + NGlo - NGlB;
        m = 0;
        cone = grd.CEL[j];
        for(k=NGlo+NArB; k<NCel-NPol; k++) {
            ctwo = grd.CEL[k];
            if( cone[0] == ctwo[0] and cone[1] == ctwo[1] ) {
                Arc.MBGlo[i] = k;
                m += 1;
                n += 1;
                if( TstOut ) ztest << i <<" "<< j <<" "<< k <<" "<< ctwo; 
                break;
            }
        }
    }
    cout << " Matched global boundary cell n/NGlB = " << n <<"/"<< NGlB << endl;

    if( ztest ) ztest << " Arctic linking boundary matchi: i, j, k, cell. \n";
//  Match Arctic boundary cells with global inner cells
    n = 0;
    for(i=0; i<NArB; i++) {
        j = i + NGlo;
        m = 0;
        cone = grd.CEL[j];
        for(k=NGlo-NGlB-2*NArB; k<NGlo-NGlB; k++) {
            ctwo = grd.CEL[k];
            if( cone[0] == ctwo[0] and cone[1] == ctwo[1] ) {
                Arc.MBArc[i] = k;
                m += 1;
                n += 1;
                if( TstOut ) ztest << i <<" "<< j <<" "<< k <<" "<< ctwo; 
                break;
            }
        }
    }
    cout << " Matched Arctic boundary cell n/NArB = " << n <<"/"<< NArB << endl;

//  Set up North and South polar cells if any, assuming polar cells are at end.
    cone = grd.CEL[NCel-1];
    if( NPol == 2 ) {
        ctwo = grd.CEL[NCel-2];
        if( cone[1] > ctwo[1] ) {
            Arc.North = NCel - 1;
            Arc.South = NCel - 2;
        }
        else {
            Arc.North = NCel - 2;
            Arc.South = NCel - 1;
        }
    }
    else if( NPol == 1 ) {
             Arc.North = NCel - 1;
             Arc.South = NCel - 1;   // Merge with North polar cell.
         }
    else {
        cout << "Inconsistent Arctic and NPol setting NPol= " << NPol << endl;
    }
    cout << " North South pole IDs as " << Arc.North <<" "<< Arc.South << endl;
    cout << " South pole is set as " << Arc.South << endl;
    if( ztest ) ztest << 
            " North South pole IDs as " << Arc.North <<" "<< Arc.South << endl;

//  Calculate Arctic part rotation angle from map-esat to local east.
    auto ZLon = grd.zlon, ZLat = grd.zlat;
    auto DLon = grd.dlon, DLat = grd.dlat;
    auto HLon = 0.5*DLon, HLat = 0.5*DLat;
    auto PLon = Arc.PLon, PLat = Arc.PLat; 

    RotateArrays Rta(NArc);
    // Use default rotated Pole for map-east definition. 
    Rta.polon = PLon, Rta.polat = PLat;

    for(k=0; k<NArc-NPol; k++) {
        i = k + NGlo; 
        cone = grd.CEL[i];
        Rta.slon[k] = cone[0]*DLon + cone[2]*HLon + ZLon;
        Rta.slat[k] = cone[1]*DLat + cone[3]*HLat + ZLat;
    }

    k = Arc.North - NGlo;
    cone = grd.CEL[Arc.North];
    // North Polar cell centre coincide with NP but roation angle is specified
    // for its edge to avoid the singularity at Pole.
    Rta.slon[k] = cone[0]*DLon + ZLon;
    Rta.slat[k] = cone[1]*DLat + ZLat;
    // Set polar cell radius in lat deg and alarm if it differs from grd one.
    Arc.PCRd = 90.0f - Rta.slat[k];
    if( abs(Arc.PCRd - grd.PCRd) > 1.0E-6f ) {
        cout << " *** Arc.PCRDY = " << Arc.PCRd << endl; 
        cout << " *** grd.PCRDY = " << grd.PCRd << endl; 
    }
    if( ztest ) ztest <<" Arc.PCRDY = "<< Arc.PCRd << endl;

    // For two polar parts, add south polar cell slon slat. 
    if( NPol == 2 ) {
        m = Arc.South - NGlo;
        cone = static_cast<celtyp>(grd.CEL[Arc.South]);
        Rta.slon[m] = cone[0]*DLon + ZLon;
        Rta.slat[m] = cone[1]*DLat + ZLat;
    }

    // Calculate the roated angles for all Polar regions.
    lltoeqangled(Rta, NArc);   // Rta passes the reference pointer not the data.

    // Also update GtCrTurning parameter Cfr.cthr with rotated latitude.
    // Save the Rta.angd to Arc.Angld.   JGLi10Jul2026
    for(k=0; k<NArc; k++) {
        Arc.Angld[k] = Rta.angd[k];
        Arc.ELats[k] = Rta.elat[k];
    }
    cout << " Arctic map-east Angld and ELats calculated. \n";

    if( ztest ) {
        // Print a few Arctic AnglD values for checking, if any.
        ztest <<" Arc.Angld[0:NArb/2], NArb = "<< NArB << endl;
        for(int i=0; i< NArB/2; i++ ) {
            ztest <<" "<< Arc.Angld[i];
            if( (i+1) % 8 == 0 ) ztest << endl;
        }
        ztest << endl;
 
    }
    ztest.close();   // Close the ztest output file stream.

    // Set up initialized flag to be true. 
    initialized = true;

    return(0); 
    // End of ArctAngd function.
}

/* Function to fetch Arctics Arc for local use.  */
const Arctics&  getArctics() {
    if (!initialized) {
        throw std::runtime_error("*** Arctics data not initialized.");
    }
    // Return the Arctics Arc data, otherwise.
    return Arc;
}

// End of Arctics.cpp module.

