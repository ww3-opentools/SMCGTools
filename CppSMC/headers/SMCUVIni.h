/* The SMCUVIni.h file declares a struct to hold deformation parameters and
 * some function interfaces for scalar test initialization.   JGLi16Oct2026 */
#pragma once
#include <vector>
#include <vector>
#include "Propgtn.h"

struct ECLonLat {
    std::size_t nsea {0};
    float AguDfm {2.59808f}; 
    float PLonDf {-178.568f}, PLatDf {81.818f}; 
    std::vector<float> ECLon, ECLat;
};

// Function to initial U V and CV values for scalar tests. 
const ECLonLat& getECLonLat();
int SMCUVIni( Propxy& Pro ); 
int SMCUVDfm( Propxy& Pro ); 
int SMCNTDfm( std::vector<float>& CNT, float TimeNT ); 

// End of SMCUVIni.h file.

