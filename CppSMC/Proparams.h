#ifndef CPP_PROPARAMS_H
#define CPP_PROPARAMS_H

#include <string>
#include <vector>

/***************************************************************************
 *  The header declares parameters used for SMC grid propagation test.
 *  To be used by #include "./Proparams.h".        JGLi16Oct2026  
 ***************************************************************************/

// Set spectral dimensions as global size_t constants.
    const std::size_t NFrq = 30;
    const std::size_t NDrc = 24;

// Output file name and default restart file saved at run end.
    std::string Flnm{}, Msgs{};
    bool SaveRestart = true;

// Restart file name and default no restart file.
    std::string RstrFile{};
    bool Restarts = false;
    int  NSt = 1;

// Default no test output files and use UNO2 scheme.
    bool TestOutf = false;
    bool UNO3 = false;
 
// Timestep and number of steps for the run and write out.
    float DT = 3600.0f;   // Global timestep in second.
    int  NTS = 24;        // Default run for a day.
    int  NWP = 6; 

// Local timestep as frequency variables to be specified.
    std::vector<int>   NTFrq(NFrq);
    std::vector<float> DTFrq(NFrq);

/* End of Proparams.h  */

#endif

