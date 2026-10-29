/***************************************************************************
 *  The SMCSclrProp.cpp program reads cell and face arrays from SMC1d text 
 *  files and does a solid rotation or deformation test.    JGLi16Oct2026  
 ***************************************************************************/
#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <cstdio>    // for C-style prontf format.
#include <string>
#include <cstring>
#include <cmath>
#include <ctime>
#include <vector>
#include <array>
#include <cassert>
#include <memory>
#include <omp.h>

#ifdef SMC61250
#include "SMC61250.h"
#elif SMC1R3
#include "SMC1R3.h"
#else  
#include "SMC1d.h"
#endif /* SMC61250 */

#include "SMCGrid.h"
#include "Arctics.h"
#include "Working.h"
#include "Rotated.h"
#include "Propgtn.h"

#include "SMCGradn.h"
#include "SMCUVIni.h"      // function interface only.
#include "CTimeDate.h"     // Time stamp print functions.
#include "PrintVector.h"   // Print a vector variable.

#include "SMConsts.h"
using namespace SMConsts;

int k=100, m=0;     // Test global and local k value.

int main(int argc, char *argv[]) {

   int k=0;
   cout << "Local k = " << k <<" Global ::k = "<< ::k << '\n';

   // Default output Msgs for PrintVector to be empty. 
   string Msgs{};
   bool UNO3 = true;

   // Open a message file to store information of this run.
   ofstream fmsg;
   fmsg.open("Msgs"+model+".txt");
   fmsg << "\n    Screen messages for "+model+" run \n\n";

   /* Function prototypes */
   double average(const std::vector<float>& x, const int nmb);  // Average function. 

   // Save ctime for time stamps, accurate to one second.
   std::time_t runstart = std::time(nullptr);
   printDateTime("Run started: ");
   std::cout << '\n';
   printDateTime("Run started: ", fmsg);
   fmsg << '\n';

   // Use input argument to select a propagation test.
   cout << "Full list of argv are: " << endl;
   for( int i=0; i<argc; ++i ) cout << i << "  " << argv[i] << endl;

   // Use first argv[1] for deformation test.
   bool Deform = false;
   if( argc > 1 ) {
       string Argv1 = argv[1];
       if( Argv1.size() > 3 ) {
           Deform = true;
           cout <<" * Deformaton test run. \n";
           fmsg <<" * Deformaton test run. \n";
       }
   }
   else {
       cout <<" * Solid rotation test. \n";
       fmsg <<" * Solid rotation test. \n";
   } 

   // Read cell array from input file, including Arctic if any.
   string inpmodel = indir + model;
   float zdlnlt[4] = {zlon, zlat, dlon, dlat};
   for( int i=0; i<4; ++i ) cout << "  " << zdlnlt[i]; 
   cout << endl;

    // Checking if OpenMP is enabled
    #ifdef _OPENMP
        cout << " * OpenMP enabled\n";
        fmsg << " * OpenMP enabled\n";
    #else
        cout << " * OpenMP NOT enabled\n";
        fmsg << " * OpenMP NOT enabled\n";
    #endif

    // Start a parallel session and output number of threads.
    #pragma omp parallel 
    {
        #pragma omp master 
        {
            #pragma omp critical
            cout << " * Number of Threads = "
                      << omp_get_num_threads() << std::endl;
            fmsg << " * Number of Threads = "
                      << omp_get_num_threads() << std::endl;
        }
        // Report from all threads.
        int tid = omp_get_thread_num();
        #pragma omp critical
            cout << "Hello from thread " << tid << std::endl;
  
    } // End of omp parallel thread check block.

   // Check whether Arctic part is present.
   if( Arctic ) cout <<" Arctic part will be included. \n";
   else      cout <<" No Arctic part will be included. \n";

   readSMCELL( inpmodel, nlv1, NCM, NPol, zdlnlt, Arctic );

   const auto& grd = getSMCELL();   // Read only grd as const reference.

   auto NCel = grd.ncel;
   auto cell = static_cast<celtyp>(grd.CEL[NCel-1]);
   cout << " * Total cell number = " << NCel << endl;
   fmsg << " * Total cell number = " << NCel << endl;
   cout << " * Last  cell arrays = ";
   cout << cell;  // Use celtyp its own output format.
   cout << " * Boundary cells appended. \n";
   for (int i=0; i<10; i++) {
       cell = static_cast<celtyp>(grd.CEL[NCel+i]);
       cout << cell; 
   }

   auto NLvl = grd.nlvl;
   int  MLFc = 1 << (NLvl - 1);
   cout << " * Number of Levels NLvl and MLFC = " << NLvl <<" "<< MLFc << endl;
   fmsg << " * Number of Levels NLvl and MLFC = " << NLvl <<" "<< MLFc << endl;

// Test using wrk1 array in Working.h as external vector.   JGLi08Jun2026 
   wrk1.resize(NCM);
   cout << " * Extracting cell depth from CEL ... \n";
   for( int i=0; i<NCel; i++ ) {
       wrk1[i] = grd.CEL[i][4];    //Defalt integer conversion to float.
       if (i == 0 || i == NCel -1) {
           cout <<"  "<< i <<": " << wrk1[i];
       }
   }
   cout << endl;

   auto zbase = average(wrk1, NCel);   // Only process NCel elements of DEPTH.
   cout << " * Mean depth = " << zbase << endl;
   fmsg << " * Mean depth = " << zbase << endl;

   auto NGlo = grd.nglo;
   cout << " * Global part number NGlo = " << NGlo << endl;

   // Read cell side/face array, including Arctic part if any.
   readSMCSID( inpmodel, nlv1, NFM, Arctic );

   // Get a gsd data.
   const auto& gsd = getSMCSID();

   auto NUfc = gsd.nufc;
   auto NVfc = gsd.nvfc;

   cout << " * Total i/j-face numbers = " << NUfc <<" "<< NVfc << endl;
   fmsg << " * Total i/j-face numbers = " << NUfc <<" "<< NVfc << endl;

   // Check i/j-side level loops.
   for( int i=0; i<nlv1; i++) {
        cout << "gsd.isdlv[" << i <<"] = "<<gsd.isdlv[i] <<endl;
        cout << "gsd.jsdlv[" << i <<"] = "<<gsd.jsdlv[i] <<endl;
   }

   auto isid = static_cast<sidtyp>(gsd.ISD[NUfc-1]);
   auto jsid = static_cast<sidtyp>(gsd.JSD[NVfc-1]);
   cout << " * Last i-side arrays = ";
   cout << isid;  // Use sidtype its own output format.
   cout << " * Last j-side arrays = ";
   cout << jsid;

   // Set up Arctics part, if any.
   Arctics Arc;
   if( grd.Arctic ) {
       cout << " * Call ArctAngd with NPol = " << NPol << endl;
       fmsg << " * Call ArctAngd with NPol = " << NPol << endl;
       bool TstOut = false;
       ArctAngd( TstOut );
       Arc = getArctics();
       cout << " * NPol passed to Arc.NPol = " << Arc.NPol << endl;
       cout << " * NArc and NGlo = " << Arc.NArc <<" "<< Arc.NGlo << endl;
       cout << " * NGlB and NArB = " << Arc.NGlB <<" "<< Arc.NArB << endl;
       fmsg << " * NArc and NGlo = " << Arc.NArc <<" "<< Arc.NGlo << endl;
       fmsg << " * NGlB and NArB = " << Arc.NGlB <<" "<< Arc.NArB << endl;
       if( Arc.NGlo != NGlo ) cout <<" *** Warning: Arc.NGlo not matching grd.nglo ! \n";
   }

   /* Solid sphere rotation test on SMC1d grid.  JGLi16Jun2026  */

   // Declare Propxy struct for propagation test.
   Propxy Pro(NCM, NFM);  // Allocate Propxy struct.

   // Check which advection scheme to be used.
   if( UNO3 ) {  // Default UNO2 scheme will be used.
       cout << " * UNO3 scheme will be used. \n";
       fmsg << " * UNO3 scheme will be used. \n";
   }

   // Set up timestep, number of steps, and writeup interval.
   float DT = 120.0;
   if( Deform ) DT = 0.02;
   int  NTS = 1080;  // For one round 36 hr rotation.
   if( Deform ) NTS = 1000;
   int  NWP = 90;   // Every 3 hrs.
   if( Deform ) NWP = 50;
   cout <<" * Starting NT = 1 until NTS = "<< NTS <<" NWP = "<< NWP << endl;
   fmsg <<" * Starting NT = 1 until NTS = "<< NTS <<" NWP = "<< NWP << endl;

   // Work out the maximum diffusity for size-1 cell.
   // Note that the diffusivity is measured in radian distance.
   float dxh = grd.dlon*D2Rad*0.5;   // Minimum dx half the Equotor one.
   float dy1 = grd.dlat*D2Rad;       // dy for size-1 is a constant.
   float dsmn = min(dxh, dy1);
// float Adfmx = 0.1/(float)MLFc;   // Maximum Fourier number = 0.1.
   float Adfmx = 0.0/(float)MLFc;   // Swithc off diffusion for SMC1d.
   cout << " * Maximum size-1 diffusion Fourier number = " << Adfmx << "\n";
   fmsg << " * Maximum size-1 diffusion Fourier number = " << Adfmx << "\n";
   Pro.ADfu = Adfmx;
   Pro.ADfv = Adfmx;
   Pro.DTFc = DT/MLFc;   // Fractional timestep for levl 1 propagation.

   // Call SMCUVIni to initialize Prx.cv, uc and Pry.vc.
   if( Deform ) SMCUVDfm( Pro );
   else         SMCUVIni( Pro );
   cout << " * Pro.cv[0] initialised as " << Pro.cv[0] << endl;
   cout << " * Pro.uc[0] initialised as " << Pro.uc[0] << endl;
   cout << " * Pro.vc[0] initialised as " << Pro.vc[0] << endl;

   // Write out the wind field.
   PrintVector(Pro.uc, NCel, 1, "UC", Msgs);
   PrintVector(Pro.vc, NCel, 2, "VC", Msgs);
   PrintVector(Pro.uu, NUfc, 3, "UU", Msgs);
   PrintVector(Pro.vv, NVfc, 4, "VV", Msgs);

   // Write out the initial state.
   cout << " * Write out at NT = 0 ";
   cout << " Pro.cv[0] = " << Pro.cv[0] << endl;
   PrintVector(Pro.cv, NCel, 0, "Cn", Msgs);

   printDateTime(" Loop started at ");
   std::cout << '\n';
   printDateTime(" Loop started at ", fmsg);
   fmsg << '\n';

   int LvF;
   std::vector<float> CNT(NCel);

   // Main timestep and level loops, all starting from 1.
   for( int NT=1; NT<=NTS; NT++ ) {
       for( int MLF=1; MLF<=MLFc; MLF++ ) {
           for( int Lvl=1; Lvl<=NLvl; Lvl++ ) {
               LvF = 1 << (Lvl - 1);
               if( MLF % LvF == 0 ) {
                   if( UNO3 ) { // Use UNO3 scheme.
                       SMCUNO3x(Lvl, Pro);
                       SMCUNO3y(Lvl, Pro);
                   }
                   else { // Default UNO2 scheme.
                       SMCUNO2x(Lvl, Pro);
                       SMCUNO2y(Lvl, Pro);
                   }
               }
           }
       }  // End of MLF loop for NT step.

       // Update Arctic linking boundary condition if any.
       if( Arctic ) ArcScalr( Pro.cv, Arc );

       if( NT % NWP == 0 ) {
           cout << " * Write out at NT = " << NT; 
           cout << " Pro.cv[0] = " << Pro.cv[0] << endl;
           PrintVector(Pro.cv, NCel, NT, "Cn", Msgs);
           if( Deform ) {
               SMCNTDfm( CNT, DT*NT );
               PrintVector(CNT, NCel, NT, "Cx", Msgs);
           }
       }
   }

   // Save the run end data and time.
   printElapsed(runstart);
   printElapsed(runstart, fmsg);

   printDateTime("Run finished: ");
   std::cout << '\n';
   printDateTime("Run finished: ", fmsg);
   fmsg << '\n';

   fmsg.close();
   cout << " * All done. \n";

   return(0);
}

/*   End of main function   */

/*  This function calculates the average  */

double average(const std::vector<float>& x, const int nmb) 
{
   int k;
   double avrg=0.0;
   for (k=0; k<nmb; k++)
       avrg += x[k];

   if (nmb > 0)
       avrg /= (float)nmb;

   return (avrg);
}
/*  End of function average  */

