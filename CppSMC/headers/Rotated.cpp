/* The Rotated module contains functions for conversion between standard 
 * latitude-longitude grid and rotated or equitorial grid.   JGLi02Jun2026   */
// Note struct wk has float type and most calculations are double precision.
// So before store results into wr, it is better converted to float first.
  
#include <iostream>
using namespace std;
#include <string>
#include <vector>
#include <cmath>

#include "Rotated.h"
#include "SMConsts.h"
using namespace SMConsts;

double small = 1.0E-8;

// Function to convert standard lat-lon to equatorial ones. 
int lltoeqangled(RotateArrays& wk, int Nmbr ) {
    
//  Check input work array demension less than Nmbr.
    if( wk.slat.size() < Nmbr ) {
        cout << " wk.slat.size() = " << wk.slat.size() << endl;
        cout << " less than Nmbr = " << Nmbr << endl;
        return(1);
    }

//  Simply copy over for Polat=90.0 and Polon=0.0 case. 
    if( wk.polat == 90.0f and wk.polon == 0.0f ) {
        for (int i=0; i<Nmbr; i++) {
            wk.elon[i] = wk.slon[i];
            wk.elat[i] = wk.slat[i];
        }
        return(0);
    }

    float plonmd = wk.polon;
    
//  1. Initialise local constants
    if( plonmd > 180.0f ) plonmd -= 360.0f; 
    float zerolon = plonmd + 180.0f;

//  Sine and cosine of latitude of eq pole
    double polatrad = wk.polat*D2Rad;
    double sinplat = sin( polatrad );
    double cosplat = cos( polatrad );
    if( wk.polat < 0.0f ) {
        sinplat *= -1.0;
        cosplat *= -1.0;
    } 

//  2. Transform from standard to equatorial latitude-longitude
    for (int i=0; i<Nmbr; i++) {

//  Scale longitude to range -180 to +180 degs
        float slod = wk.slon[i] - zerolon;
        if( slod >   180.0f ) slod -= 360.0f;
        if( slod <= -180.0f ) slod += 360.0f;

//  Convert latitude & longitude to radians
        double slon = D2Rad*slod;
        double slat = D2Rad*wk.slat[i];

//  Compute eq latitude using equation (4.4)
        double arg = -cosplat*cos(slat)*cos(slon) + sinplat*sin(slat);
        arg = min(arg, 1.0);
        arg = max(arg,-1.0);
        double EPHI = asin(arg);
        wk.elat[i] = static_cast<float>(Rad2D*EPHI);

//  Compute eq longitude using equation (4.6)
        double TERM1 = sinplat*cos(slat)*cos(slon) + cosplat*sin(slat);
        double TERM2 = cos(EPHI);
        double ELMD = 0.0;
        if (abs(TERM2) > 0.0) {
            arg = TERM1/TERM2;
            arg = min(arg, 1.0);
            arg = max(arg,-1.0);
            ELMD = Rad2D * acos(arg); 
            ELMD = copysign(ELMD, slon);
        }

//  Scale longitude to range 0 to 360 degs
        if (ELMD >= 360.0)  ELMD -= -360.0;
        if (ELMD <    0.0)  ELMD +=  360.0;
        wk.elon[i] = static_cast<float>(ELMD);

//  Calculate turning angle for standard wind velocity
        double elmdrad = D2Rad*ELMD;

//  Formulae used are from eqs (4.19) and (4.21)
        TERM2 = sin(elmdrad);
        arg = sin(slon)*TERM2*sinplat + cos(slon)*cos(elmdrad);
        arg = min(arg, 1.0);
        arg = max(arg,-1.0);
        TERM1 = Rad2D * acos(arg);
        wk.angd[i] = static_cast<float>( copysign(TERM1, TERM2) );

    }  // End of i-loop

    return(0); 

    // End of lltoeqangled function.
}


// Function eqtollangled to reverse equatorial lat-lon deg to standard ones. 
// Note struct wk has float type and most calculations are double precision.
// So before store results into wr, it is better converted to float first.
  
int eqtollangled(RotateArrays& wk, int Nmbr ) {
    
//  Check input work array demension less than Nmbr.
    if( wk.elat.size() < Nmbr ) {
        cout << " wk.elat.size() = " << wk.elat.size() << endl;
        cout << " less than Nmbr = " << Nmbr << endl;
        return(1);
    }

//  Simply copy over if Polat=90.0 and Polon=0.0 
    if( wk.polat == 90.0f and wk.polon == 0.0f ) {
        for (int i=0; i<Nmbr; i++) {
            wk.slon[i] = wk.elon[i]; 
            wk.slat[i] = wk.elat[i]; 
        }
        return(0);
    }

    float plonmd = wk.polon;
    if( plonmd > 180.0f ) plonmd -= 360.0f;
    
//  1. Initialise local constants
    float zerolon = plonmd + 180.0f;

//  Sine and cosine of latitude of eq pole
    double polatrad = wk.polat*D2Rad;
    double sinplat = sin( polatrad );
    double cosplat = cos( polatrad );
    if( wk.polat < 0.0f ) {
        sinplat *= -1.0;
        cosplat *= -1.0;
    } 

//  2. Transform from equatorial to standard latitude-longitude
    for (int i=0; i<Nmbr; i++) {

//  Scale longitude to range -180 to +180 degs
        float elod = wk.elon[i];
        if( elod >=  180.0f ) elod -= 360.0f;
        if( elod <  -180.0f ) elod += 360.0f;

//  Convert latitude & longitude to radians
        double elon = D2Rad*elod;
        double elat = D2Rad*wk.elat[i];

//  Compute eq latitude using equation (4.7)
        double arg = cosplat*cos(elat)*cos(elon) + sinplat*sin(elat);
        arg = min(arg, 1.0);
        arg = max(arg,-1.0);
        double EPHI = asin(arg);
        wk.slat[i] = static_cast<float>(Rad2D*EPHI);

//  Compute eq longitude using equation (4.8)
        double TERM1 = sinplat*cos(elat)*cos(elon) - cosplat*sin(elat);
        double TERM2 = cos(EPHI);
        double SLMD = 0.0;
        if (abs(TERM2) > 0.0) {
            arg = TERM1/TERM2;
            arg = min(arg, 1.0);
            arg = max(arg,-1.0);
            SLMD = Rad2D * acos(arg); 
            SLMD = copysign(SLMD, elon);
            SLMD += zerolon;
        }

//  Scale longitude to range 0 to 360 degs
        if (SLMD >= 360.0)  SLMD -= 360.0;
        if (SLMD <    0.0)  SLMD += 360.0;
        wk.slon[i] = static_cast<float>(SLMD);

//  Calculate turning angle for standard wind velocity
        double slmdrad = D2Rad*(SLMD - zerolon);

//  Formulae used are from eqs (4.19) and (4.21)
        TERM2 = sin(elon);
        arg = sin(slmdrad)*TERM2*sinplat + cos(slmdrad)*cos(elon);
        arg = min(arg, 1.0);
        arg = max(arg,-1.0);
        TERM1 = Rad2D * acos(arg);
        wk.angd[i] = static_cast<float>( copysign(TERM1, TERM2) );

    }  // End of i-loop

    return(0); 

    // End of eqtollangled function.
}

