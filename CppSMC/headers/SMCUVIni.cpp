/* The SMCUVIni module initializes Pro.cv and u-v CFL values on a SMC grid. 
 * It works on portions of u-faces for multi-resolution use.  JGLi09Jun2026 */

#include <omp.h>
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>
#include <algorithm>

#include "SMCGrid.h"
#include "Rotated.h"
#include "Propgtn.h"
#include "SMCUVIni.h"

#include "SMConsts.h"
using namespace SMConsts;

namespace {     // Anonymous namespace only available to this file.
    struct ECLonLat ECLL;
    bool initialized = false;
}

// Function to initialize SMC grid cell value and velocity for propagation.
int SMCUVIni( Propxy& Pro ) {
    
    cout << " Initializaton for solid rotation test. \n";
    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    const auto& gsd = getSMCSID();   // Get gsd as read only const reference.
    auto NCel = grd.ncel;
    auto NUfc = gsd.nufc;
    auto NVfc = gsd.nvfc;
    cout << " Total number of cells grd.ncel = " << NCel << endl;
    cout << " Total number of USide gsd.nufc = " << NUfc << endl;
    cout << " Total number of VSide gsd.nvfc = " << NVfc << endl;
    auto zlon = grd.zlon;
    auto zlat = grd.zlat;
    auto dlon = grd.dlon;
    auto dlat = grd.dlat;

    //  Check nlvl and npol from grd struct.
    auto NLvl = grd.nlvl; 
    auto NPol = grd.npol;
    cout << " Number of levels = " << NLvl << endl;
    cout << " Number of poles  = " << NPol << endl;

    // Factor of maxium level cell size.
    int MFct = 1 << (NLvl - 1);
    cout << " Maximum level MFC = " << MFct << endl;
 
    // Initial Pro.cv background values to be 1.0.
    Pro.cv.assign(NCel, 1.0f);

    // Insert a belt 4.0 above background for cells between abs(lat < 11.0).
    celtyp cell;
    for (int i=0; i<NCel; i++) {
        cell = grd.CEL[i];
        float cnst = (float)(zlat + (cell[1]+0.5*cell[3])*dlat);
        if( abs(cnst) < 11.0 ) Pro.cv[i] += 4.0f;
    }

    // Setup rotated struct for cell and face wind initialization.
    int MSid = max(NUfc, NVfc);
    RotateArrays ra(MSid);
    ra.polon = -180.0, ra.polat = 0.0;

    // Couriant mumber factor for size-1 radian distance.
    float cnst1 = Pro.DTFc/(dlon*D2Rad);
    float cnst2 = Pro.DTFc/(dlat*D2Rad);

    // Initialize cell centre Courant number for size-1 cell.
    for (int i=0; i<NCel; i++) {
        cell = grd.CEL[i];
        ra.slat[i] = (float)(zlat + (cell[1]+0.5*cell[3])*dlat);
        ra.slon[i] = (float)(zlon + (cell[0]+0.5*cell[2])*dlon);
        ra.slat[i] = (float)(zlat + (cell[1]+0.5*cell[3])*dlat);
//      if( i >= NCel-NPol ) ra.slat[i] = (float)(zlat + (cell[1])*dlat);
    }
    cout << " Calling lltoeqangled ... \n";
    lltoeqangled(ra, NCel);   // ra passes the reference pointer not the data.

    cout << " Calculating cell centre velocities Ang36H = " << Ang36H << "\n"; 
    for (int i=0; i<NCel; i++) {
        float omega = Ang36H*cos(ra.elat[i]*D2Rad);
        Pro.uc[i] =  omega*cos(ra.angd[i]*D2Rad)*cnst1/grd.CLat[i];
        Pro.vc[i] = -omega*sin(ra.angd[i]*D2Rad)*cnst2;
    }
    
    // Initialize U-side centre Courant number for size-1 distance.
    sidtyp sdu;
    for (int i=0; i<NUfc; i++) {
        sdu = gsd.ISD[i];
        ra.slon[i] = (float)(zlon + sdu[0]*dlon);
        ra.slat[i] = (float)(zlat + (sdu[1]+0.5*sdu[2])*dlat);
    }
    cout << " Calling lltoeqangled ... \n";
    lltoeqangled(ra, NUfc);   // ra passes the reference pointer not the data.

    cout << " Calculating U-side centre velocities Ang36H = " << Ang36H << "\n"; 
    for (int i=0; i<NUfc; i++) {
        float omega = Ang36H*cos(ra.elat[i]*D2Rad);
        Pro.uu[i] =  omega*cos(ra.angd[i]*D2Rad)*cnst1/cos(ra.slat[i]*D2Rad);
    }
    
    // Initialize V-side centre Courant number for size-1 distance.
    sidtyp sdv;
    for (int j=0; j<NVfc; j++) {
        sdv = gsd.JSD[j];
        ra.slon[j] = (float)(zlon + (sdv[0]+0.5*sdv[2])*dlon);
        ra.slat[j] = (float)(zlat + (sdv[1])*dlat);
    }
    cout << " Calling lltoeqangled ... \n";
    lltoeqangled(ra, NVfc);   // ra passes the reference pointer not the data.

    cout << " Calculating V-side centre velocities Ang36H = " << Ang36H << "\n"; 
    for (int j=0; j<NVfc; j++) {
        float omega = Ang36H*cos(ra.elat[j]*D2Rad);
        Pro.vv[j] = -omega*sin(ra.angd[j]*D2Rad)*cnst2;
    }
    
    return(0); 

    // End of SMCUVIni function.
}

// Function to initialize cell value and velocity for deformation test.
int SMCUVDfm( Propxy& Pro ) {
    
    cout << " Initializaton for deformation test. \n";
    // Fetch the cell array and parameters. 
    const auto& grd = getSMCELL();   // Get grd as read only const reference.
    const auto& gsd = getSMCSID();   // Get gsd as read only const reference.
    auto NCel = grd.ncel;
    auto NUfc = gsd.nufc;
    auto NVfc = gsd.nvfc;
    cout << " Total number of cells grd.ncel = " << NCel << endl;
    cout << " Total number of USide gsd.nufc = " << NUfc << endl;
    cout << " Total number of VSide gsd.nvfc = " << NVfc << endl;
    auto zlon = grd.zlon;
    auto zlat = grd.zlat;
    auto dlon = grd.dlon;
    auto dlat = grd.dlat;

    //  Check nlvl and npol from grd struct.
    auto NLvl = grd.nlvl; 
    auto NPol = grd.npol;
    cout << " Number of levels = " << NLvl << endl;
    cout << " Number of poles  = " << NPol << endl;

    // Factor of maxium level cell size.
    int MFct = 1 << (NLvl - 1);
    cout << " Maximum level MFC = " << MFct << endl;
 
    // Initial Pro.cv background values to be 0.0.
    Pro.cv.assign(NCel, 0.0f);

    // Set ECLL dimension nsea and delcare ECLon and ECLat.
    ECLL.nsea = NCel;
    ECLL.ECLon.assign(NCel, 0.0f);
    ECLL.ECLat.assign(NCel, 0.0f);

    // Setup rotated struct for cell and face wind initialization.
    int MSid = max(NUfc, NVfc);
    RotateArrays ra(MSid);
    ra.polon = ECLL.PLonDf, ra.polat = ECLL.PLatDf;

    // Couriant mumber factor for size-1 radian distance.
    float cnst1 = Pro.DTFc/(dlon*D2Rad);
    float cnst2 = Pro.DTFc/(dlat*D2Rad);

    // Initialize cell centre Courant number for size-1 cell.
    celtyp cell;
    for (int i=0; i<NCel; i++) {
        cell = grd.CEL[i];
        ra.slon[i] = (float)(zlon + (cell[0]+0.5*cell[2])*dlon);
        ra.slat[i] = (float)(zlat + (cell[1]+0.5*cell[3])*dlat);
//      if( i >= NCel-NPol ) ra.slat[i] = (float)(zlat + (cell[1])*dlat);
    }
    cout << " Calling lltoeqangled ... \n";
    lltoeqangled(ra, NCel);   // ra passes the reference pointer not the data.

    float Alpha, Omega, DfSpd, AguDfm = ECLL.AguDfm; 
    cout << " Calculating cell centre velocities AguDfm = " << AguDfm << "\n"; 
    for (int i=0; i<NCel; i++) {
        Alpha = 3.0f*cos(ra.elat[i]*D2Rad);
        if( Alpha == 0.0f ) {
            Omega = AguDfm;
            DfSpd = 0.0; 
        }
        else {
            Omega = AguDfm*sinh(Alpha)/( Alpha*pow(cosh(Alpha), 3) );
            DfSpd = Omega*Alpha/3.0f; 
        }    
        Pro.uc[i] =  DfSpd*cos(ra.angd[i]*D2Rad)*cnst1/grd.CLat[i];
        Pro.vc[i] = -DfSpd*sin(ra.angd[i]*D2Rad)*cnst2;
        Pro.cv[i] = 1.0f - tanh( 0.2f*Alpha*sin(ra.elon[i]*D2Rad) );
        // Save a copy of elon and elat to ECLL.
        ECLL.ECLon[i] = ra.elon[i];
        ECLL.ECLat[i] = ra.elat[i];
    }
    
    // Initialize U-side centre Courant number for size-1 distance.
    sidtyp sdu;
    for (int i=0; i<NUfc; i++) {
        sdu = gsd.ISD[i];
        ra.slon[i] = (float)(zlon + sdu[0]*dlon);
        ra.slat[i] = (float)(zlat + (sdu[1]+0.5*sdu[2])*dlat);
    }
    cout << " Calling lltoeqangled ... \n";
    lltoeqangled(ra, NUfc);   // ra passes the reference pointer not the data.

    cout << " Calculating U-side centre velocities AguDfm = " << AguDfm << "\n"; 
    for (int i=0; i<NUfc; i++) {
        Alpha = 3.0f*cos(ra.elat[i]*D2Rad);
        if( Alpha == 0.0f ) {
            Omega = AguDfm;
            DfSpd = 0.0;
        }
        else {
            Omega = AguDfm*sinh(Alpha)/( Alpha*pow(cosh(Alpha), 3) );
            DfSpd = Omega*Alpha/3.0f;
        }
        Pro.uu[i] =  DfSpd*cos(ra.angd[i]*D2Rad)*cnst1/cos(ra.slat[i]*D2Rad);
    }
    
    // Initialize V-side centre Courant number for size-1 distance.
    sidtyp sdv;
    for (int j=0; j<NVfc; j++) {
        sdv = gsd.JSD[j];
        ra.slon[j] = (float)(zlon + (sdv[0]+0.5*sdv[2])*dlon);
        ra.slat[j] = (float)(zlat + (sdv[1])*dlat);
    }
    cout << " Calling lltoeqangled ... \n";
    lltoeqangled(ra, NVfc);   // ra passes the reference pointer not the data.

    cout << " Calculating V-side centre velocities AguDfm = " << AguDfm << "\n"; 
    for (int j=0; j<NVfc; j++) {
        Alpha = 3.0f*cos(ra.elat[j]*D2Rad);
        if( Alpha == 0.0f ) {
            Omega = AguDfm;
            DfSpd = 0.0;
        }
        else {
            Omega = AguDfm*sinh(Alpha)/( Alpha*pow(cosh(Alpha), 3) );
            DfSpd = Omega*Alpha/3.0f;
        }
        Pro.vv[j] = -DfSpd*sin(ra.angd[j]*D2Rad)*cnst2;
    }
    
    // Set up initialized flag to be true for Struct ECLL. 
    initialized = true;

    return(0); 

    // End of SMCUVDfm function.
}

// Function to calculate the exact deformation solution at a given timestep. 
int SMCNTDfm( std::vector<float>& CNT, float Time ) {
    
    // Set ECLL dimension nsea and delcare ECLon and ECLat.
    float NSea = ECLL.nsea;
    // ECLL.ECLon.assign(NCel, 0.0f);
    // ECLL.ECLat.assign(NCel, 0.0f);
    if( NSea > CNT.size() ) {
        cout << CNT.size() << " is smaller than NSea: "<< NSea <<'\n';
        return(1);
    }
        
    float Alpha, Omega, DfSpd, AguDfm = ECLL.AguDfm; 
    cout << " Calculating analytical solution with AguDfm = " << AguDfm 
         << " at time = "<< Time <<'\n'; 
    for (int i=0; i<NSea; i++) {
        Alpha = 3.0f*cos(ECLL.ECLat[i]*D2Rad);
        if( Alpha == 0.0f ) {
            Omega = AguDfm;
        }
        else {
            Omega = AguDfm*sinh(Alpha)/( Alpha*pow(cosh(Alpha), 3) );
        }    
        CNT[i]=1.0f - tanh(0.2f*Alpha*sin(ECLL.ECLon[i]*D2Rad - Omega*Time));
    }

    return(0); 

    // End of SMCNTDfm function.
}

/* Function to fetch ECLonLat ECLL for local use.  */
const ECLonLat&  getECLonLat() {
    if (!initialized) {
        throw std::runtime_error("*** ECLonLat data not initialized.");
    }
    // Return the ECLonLat ECLL data, otherwise.
    return ECLL;
}


