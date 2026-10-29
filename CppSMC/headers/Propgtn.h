/* This Propgtn.h header file declares a group of parameters used for the  
 * propagation modules, including SMCUNOx/y and SMCUVIni.    JGLi06Jun2026  */

#pragma once
#include <vector>

struct Propxy {
    float ADfu, ADfv, DTFc;
    std::vector<float> cv, cd, fx, fy, uc, vc, fxad, uu, vv;
    Propxy(int ncm, int nfm): 
        cv(ncm), cd(ncm), fx(ncm), fy(ncm), fxad(ncm), 
        uc(ncm), vc(ncm), uu(nfm), vv(nfm) {}
};

int SMCUNO2x(int Lvl, Propxy& pro, bool L0r1 = true );
int SMCUNO2y(int Lvl, Propxy& pro, bool L0r1 = true );
int SMCUNO3x(int Lvl, Propxy& pro, bool L0r1 = true );
int SMCUNO3y(int Lvl, Propxy& pro, bool L0r1 = true );

// End of Propgtn.h.

