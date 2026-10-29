/***************************************************************************
 *  The SMCGridProp.cpp program reads cell and face arrays from a given grid 
 *  files and does a spectral progagation test.           JGLi16Jun2026  
 *  Last modified for SMCGTools propagation test model.   JGLi16Oct2026
 ***************************************************************************/
#include <iostream>
#include <iomanip>
using namespace std;
#include <fstream>
#include <cstdio>    // for C-style printf format.
#include <string>
#include <cstring>
#include <vector>
#include <array>
#include <cmath>
#include <ctime>
#include <cassert>
#include <memory>
#include <chrono>
#include <omp.h>

#ifdef SMC36125
#include "SMC36125.h"   // SMC36125 grid 4-level grid.
#elif SMC125120
#include "SMC125120.h"  // SMC125120 5-level grid.
#elif SMC61250
#include "SMC61250.h"   // SMC61250 4-level grid.
#elif SWEs1d
#include "SWEs1d.h"     // Full global 1 deg grid.
#elif SMC1R3
#include "SMC1R3.h"     // Full global 1-deg 3-level grid.
#else  
#include "SMC61250.h"   // Default SMC61250 grid parameters.
#endif /* SMC36125 */

#include "SMCGrid.h"    // SMC grid parameter header.
#include "Arctics.h"    // SMC Arctic part parameters.
#include "Working.h"    // Temporary working vectors.
#include "Rotated.h"    // Rotated grid functions.
#include "Propgtn.h"    // Scalar propagation variables.

#include "SMCGradn.h"   // Gradient for SMC grid cells.
#include "Spectrum.h"   // Const spectral parameters. 
#include "WvCgrpKs.h"   // Wave number and group speed.
#include "WvSpectr.h"   // Full domain wave spectra.
#include "GtCrTurn.h"   // Gt circle turning header.
#include "SMCAverg.h"   // 1-2-1 Average on SMC grid.
#include "SMCUVIni.h"   // Solid rotation velocity field.
#include "Restrfile.h"  // Restart file write/read.
#include "CTimeDate.h"     // Time stamp print functions.
#include "PrintVector.h"   // Print a vector variable.

#include "SMConsts.h"      // Some constants for share.
using namespace SMConsts;

#include "./Proparams.h"   // Configuration parameters.

int k=100, m=0;     // Test global and local k value.

int main(int argc, char *argv[]) {

   int k=0;    // Local k value.
   cout << "Local k = " << k <<" Global ::k = "<< ::k << '\n';

   // Open a message file to store information of this run.
   ofstream fmsg;
   fmsg.open("Msgs"+model+".txt");
   fmsg << "\n    Screen messages for "+model+" run \n\n";
 
   /* Function prototypes */
   double average(const std::vector<float>& x, const int nmb);  // Average function. 
   void printarow(int i, float a, float b, float c, float d, 
                  std::ostream& ostrm = std::cout);

   // Use chrono to time the program, precision with fractional second.
   std::chrono::time_point<std::chrono::steady_clock> Startm, StLoop, Endstm;
   Startm = std::chrono::steady_clock::now();

   // Save ctime for time stamps, accurate to one second.
   std::time_t runstart = std::time(nullptr);
   printDateTime("Run started: ");
   std::cout << '\n';
   printDateTime("Run started: ", fmsg);
   fmsg << '\n';

   // Use input argument to select a propagation test.
   cout << "Full list of argv are: " << endl;
   for( int i=0; i<argc; ++i ) cout << i << "  " << argv[i] << endl;

   // Use first argv[1] for restart flag, default no restart file.
   if( argc > 1 ) {
       RstrFile = argv[1];
       ifstream Rstf(RstrFile);
       if( Rstf ) {
           Restarts = true;
           cout <<" * Run will restart with file "+RstrFile << endl;
           fmsg <<" * Run will restart with file "+RstrFile << endl;
           Rstf.close();
       }
   }

   // Read cell array from input file, including Arctic if any.
   string inpmodel = indir + model;
   float zdlnlt[4] = {zlon, zlat, dlon, dlat};
   printarow(NLV, zlon, zlat, dlon, dlat);
   printarow(NLV, zlon, zlat, dlon, dlat, fmsg);

   // Set number of threads for OpenMP parallelization.
// omp_set_num_threads(8);

    // Checking if OpenMP is enabled
    #ifdef _OPENMP
        cout << " * OpenMP enabled\n";
        fmsg << " * OpenMP enabled\n";
    #else
        cout << " * OpenMP NOT enabled\n";
        fmsg << " * OpenMP NOT enabled\n";
    #endif

//  // Check number of threads in OpenMP parallel mode.
    #pragma omp parallel 
    {
        #pragma omp master 
        {
            #pragma omp critical
            cout << " Number of Threads = "
                      << omp_get_num_threads() << std::endl;
            fmsg << " Number of Threads = "
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
           cout << i << "  " << wrk1[i] << endl;
       }
   }

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
        cout << " * gsd.isdlv[" << i <<"] = "<<gsd.isdlv[i] << 
                "   gsd.jsdlv[" << i <<"] = "<<gsd.jsdlv[i] <<endl;
   }

   auto isid = static_cast<sidtyp>(gsd.ISD[NUfc-1]);
   auto jsid = static_cast<sidtyp>(gsd.JSD[NVfc-1]);
// for (int i=2; i<6; i++) isid[i] -= NGlo;   
   cout << " * Last i-side arrays = ";
   cout << isid;  // Use sidtype its own output format.
// for (int j=2; j<6; j++) jsid[j] -= NGlo;   
   cout << " * Last j-side arrays = ";
   cout << jsid;

   // Set up wave spectral parameters.
   setSpectrum( NFrq,  NDrc );
   const auto& spm = getSpectrum();   // Read only spectral parameters.
   float Cgmx = 0.6*Gravty/spm.sigma[0];
   cout << " * Maximum group speed at gravity depth Cgmx  = " << Cgmx << endl;
   cout << " * Spectral parameters set up spm.DThta = " << spm.DThta << endl;
   if( TestOutf ) {
        Msgs = " Spectral parameters spm.theta[0:NDrc]";
        PrintVector(spm.theta, NDrc, 0, "Ta", Msgs);
        Msgs = " Spectral parameters spm.cthta[0:NDrc]";
        PrintVector(spm.cthta, NDrc, 1, "Tc", Msgs);
        Msgs = " Spectral parameters spm.sthta[0:NDrc]";
        PrintVector(spm.sthta, NDrc, 2, "Ts", Msgs);
   }

   // Set up Arctics part and fetch Arc, if any.
   Arctics Arc; 
   if( grd.Arctic ) {
       cout << " * Call ArctAngd with NPol = " << NPol << endl; 
       fmsg << " * Call ArctAngd with NPol = " << NPol << endl; 
       ArctAngd( TestOutf );
       Arc = getArctics();
       cout << " * NArc and NGlo = " << Arc.NArc <<" "<< Arc.NGlo << endl;
       cout << " * NGlB and NArB = " << Arc.NGlB <<" "<< Arc.NArB << endl;
       fmsg << " * NArc and NGlo = " << Arc.NArc <<" "<< Arc.NGlo << endl;
       fmsg << " * NGlB and NArB = " << Arc.NGlB <<" "<< Arc.NArB << endl;
       if( Arc.NGlo != NGlo ) cout <<" *** Warning: Arc.NGlo not matching grd.nglo ! \n";
   }

   // Test function SMCAverg module, re-using wrk1 vector.
   for( int i=0; i<NCel; i++ ) {
       wrk1[i] = grd.CEL[i][4];    //Defalt integer conversion to float.
   }
   cout << " wrk1[10] = " << wrk1[10] << endl;
   cout << " * Calling SMCAverg module ... \n";
   SMCAverg( wrk1 );
   cout << " wrk1[10] = " << wrk1[10] << endl;

   // Declare the GtCrTn structure and initialize Cfr.
   GtCrTn Cfr(NCel, NDrc);
   GtCrTIni( Cfr );

   // Declare and initialize wave number and group speed struct.
   WvCgrpKs CgK(NCM, NFrq); 
   WvCGrpKIni( CgK ); 
   cout << " * Wave group speed CgK.Cgrp[0][0] = " << CgK.Cgrp[0][0] << endl;
   // Save a copy of wave group speed and number K for first frequency.
   if( TestOutf ) {
       wrk2.resize(NCM);   // Reuse wrk1 as well.
       for(int i=0; i<NCel; i++) {
           wrk1[i] = CgK.Cgrp[i][0];
           wrk2[i] = CgK.WnmK[i][0];
       }
       Msgs = " Group speed at all sea points.";
       PrintVector(wrk1, NCel, 0, "Cg", Msgs);
       Msgs = " Wave number at all sea points.";
       PrintVector(wrk2, NCel, 1, "WK", Msgs);
       cout <<" * Wave group speed and number K saved. \n";
   }

   // Create full domain spectral storage class.
   WvSpectr Spc(NCel, NFrq, NDrc);

   auto MSpc = Spc.MSpctr();
   auto NSea = Spc.NSeapt();

   cout <<" NFrq NDrc MSpc = "<< NFrq <<" "<< NDrc <<" "<< MSpc << endl;
   fmsg <<" NFrq NDrc MSpc = "<< NFrq <<" "<< NDrc <<" "<< MSpc << endl;
 
   // Read restart file if run restarts.
   if( Restarts ) {
       RestrRead( Spc, NSt, RstrFile );
   }
   else { // Initialize Spc.
       // Initialize wave spectral field.
       Spectini( Spc );
       cout << " Initialise spectral Spc(0, 0) = "<< Spc(0, 0) << endl;
       cout << " And equivalently as Spc[0][0] = "<< Spc[0][0] << endl;
   } 

   // Spectral propagation test on SMC grid with Arctic part.  JGLi08Jul2026  

   // Declare Propxy struct for propagation test.
   Propxy Pro(NCM, NFM);  // Allocate Propxy struct.

   cout <<" Local timesteps with DT = "<< DT << " s, and CFL = "
        << CFL_Advctn <<'\n';
   fmsg <<" Local timesteps with DT = "<< DT << " s, and CFL = "
        << CFL_Advctn <<'\n';

   // Work out maximum CFL number and decide local timestep loop count. 
   float dxh = grd.dlon*D2Rad*0.5;   // Minimum dx half the Equotor one.
   float dy1 = grd.dlat*D2Rad;       // dy for size-1 is a constant.
   float dsmn = min(dxh, dy1);       // Minimum cell side length in radian.
   int   NTLc;
   float DTLc, CFLcx, CFLcy; 

   // Define local timestep for each frequency.
   cnst = 0.5f/Pie;    // Factor to convert Sigma to f in Hz.
   cout << " M NTLc DTLc  DTC06  Cgmx m/s  f Hz \n";
   fmsg << " M NTLc DTLc  DTC06  Cgmx m/s  f Hz \n";
   for(int M=0; M<NFrq; M++) {
       float Cgmx = 0.6*Gravty/spm.sigma[M];
       float cflmx = Cgmx*DT/(dsmn*MLFc*REarth); // Max Courant number with DT.
       float DTC06 = CFL_Advctn*DT/cflmx;   // Local timestep limited by CFL.
       NTLc = static_cast<int>( DT/DTC06 ) + 1;   // Rounded local timestep.
       DTLc = DT/(float)NTLc;    // Local timestep as integer fraction of DT. 
       NTFrq[M] = NTLc;
       DTFrq[M] = DTLc;    // Timestep for base resolution cells.

       cout << setw(3) << M;
       printarow(NTLc, DTLc, DTC06, Cgmx, spm.sigma[M]*cnst);
       fmsg << setw(3) << M;
       printarow(NTLc, DTLc, DTC06, Cgmx, spm.sigma[M]*cnst, fmsg);
   }

   // Work out the maximum diffusity Fourier number for size-1 cell.
   float Adfmx = Fourier_Df/(float)MLFc;  // Limited by maximum Fourier number.
   Pro.ADfu = Adfmx;
   Pro.ADfv = Adfmx;
   cout << " Size-1 diffusion Fourier number = " << Adfmx <<
           " and maximum Fourier number = "<< Fourier_Df <<'\n';
   fmsg << " Size-1 diffusion Fourier number = " << Adfmx <<
           " and maximum Fourier number = "<< Fourier_Df <<'\n';

   // Some parameters for timestep loop.
   float cflx = 1.0/(grd.dlon*D2Rad*REarth);
   float cfly = 1.0/(grd.dlat*D2Rad*REarth);
   bool L0r1 = false;  // Use cell centre CFL for propagation.
 
   // Test Arctic part linking boundary conditions.
   if( grd.Arctic ) {
       ArctBndy( Spc, Arc, 0, TestOutf );
   }

   // Declare SWHt and write out the initial SWHt if NSt = 1.
   std::vector<float> SWHt(NSea);
   if( NSt == 1 ) {
       SWHt = SWHeight( Spc );
       Msgs = " Initial SWHt all sea points."; 
       PrintVector(SWHt, NSea, NSt-1, "wh", Msgs);
       cout << NSt-1 << " SWHt[0] = " << SWHt[0] << endl;
   }

   cout <<" * Starting NT = "<< NSt <<" until NSt+NTS = " << NSt+NTS << endl;
   fmsg <<" * Starting NT = "<< NSt <<" until NSt+NTS = " << NSt+NTS << endl;

   int LvF, kfd;   // Sub-level factor and spectral component index.
   vector<float> Cgrpkf(NCM);

   // Mark the loop start time.
   StLoop = std::chrono::steady_clock::now();

   printDateTime(" Loop started: ");
   std::cout << '\n';

   // Main timestep and level loops, all starting from 1.
   for(int NT=NSt; NT<NSt+NTS; NT++ ) {
       for(int MF=0; MF<NFrq; MF++) {
           // Setup fractional timestep and local loop count for MF frequency.
           NTLc = NTFrq[MF];
           DTLc = DTFrq[MF];
           for(int i=0; i<NCM; i++) {
               Cgrpkf[i] = CgK.Cgrp[i][MF];
           }
           // Fractional timestep for levl 1 propagation.
           Pro.DTFc = DTLc/(float)MLFc;
           CFLcx = Pro.DTFc*cflx;
           CFLcy = Pro.DTFc*cfly;

           if( NT == NSt ) cout <<" MF = "<< MF <<" NTLc = "<< NTLc << endl;
           for(int NTL=0; NTL<NTLc; NTL++ ) {
           for(int LD=0; LD<NDrc; LD++ ) {

               // Reset cell variables to zero.
               Pro.cv.assign(NCM, 0.0f);
               Pro.cd.assign(NCM, 0.0f);
               Pro.uc.assign(NCM, 0.0f);
               Pro.vc.assign(NCM, 0.0f);

               // Use WvSpectr::Compnt(kfd) to pass component to Pro.cv.
               kfd = MF*NDrc + LD;
               auto Spckfd = Spc.Compnt( kfd );

               for(int i=0; i<NGlo; i++) {
                   Pro.cv[i] = Spckfd[i]; 
                   Pro.uc[i] = Cgrpkf[i]*spm.cthta[LD]*CFLcx/grd.CLat[i];
                   Pro.vc[i] = Cgrpkf[i]*spm.sthta[LD]*CFLcy;
               } // End of global part assignment.
               if( grd.Arctic ) {
                   for(int i=NGlo; i<NCel; i++) {
                       Pro.cv[i] = Spckfd[i]; 
                       // Exclude Polar cells as no velocity is defined.
                       if( i < NCel-NPol ) {
                           cnst0 = spm.theta[LD] - Arc.Angld[i-NGlo]*D2Rad;
                           Pro.uc[i] = Cgrpkf[i]*cos(cnst0)*CFLcx/grd.CLat[i];
                           Pro.vc[i] = Cgrpkf[i]*sin(cnst0)*CFLcy;
                       }
                   }
               } // End of Arctic assignment.

               // Reset face variables to zero.
               Pro.fxad.assign(NCM, 0.0f);
               Pro.fx.assign(NCM, 0.0f);
               Pro.fy.assign(NCM, 0.0f);

               // Calculate the propagation of this spectral component.     
               for(int MLF=1; MLF<=MLFc; MLF++ ) {
                   if( (NT == NSt and NTL == 0) and (MF == 0 and LD == 0) ) cout <<" "<< MLF <<" :";
                   for(int Lvl=1; Lvl<=NLvl; Lvl++ ) {
                       LvF = 1 << (Lvl - 1);
                       if( MLF % LvF == 0 ) {
                           if( UNO3 ) { // Use UNO3 scheme.
                               SMCUNO3x(Lvl, Pro, L0r1);
                               SMCUNO3y(Lvl, Pro, L0r1);
                           }
                           else { // Default UNO2 scheme.
                               SMCUNO2x(Lvl, Pro, L0r1);
                               SMCUNO2y(Lvl, Pro, L0r1);
                           }
                           if( (NT == NSt and NTL == 0) and (MF == 0 and LD == 0) ) cout <<" "<< LvF;
                       }
                   } // End of Lvl loop.
                   if( (NT == NSt and NTL == 0) and (MF == 0 and LD == 0) ) cout << endl;
               } // End of MLF loop.

               // Call 1-2-1 average if at NWP timestep.  JGLi22Sep2026
               if( NT % NWP == 0 )  SMCAverg( Pro.cv );

               // Pass propagated component back to Spc.
               Spc.setCompnt( kfd, Pro.cv );

           } // End of directional LD loop.
           } // End of local NTLc loop.

       } // End of NFrq loop.

       // GtCrTurning and refraction if any.
       GtCrTurn( Spc, Cfr, CgK, DTLc );

       // Update Arctic part linking boundary conditions at NT step.
       if( grd.Arctic ) {
           ArctBndy( Spc, Arc, NT );
       }

//     Note the % operation is done with immediate number so NWP/2 needs to 
//     be put within a pair of paratheses if NT is modulo with NWP/2. 
       if( NT == NSt || NT % (NWP/2) == 0 ) {
           printTime();
           cout << " NT = "<< NT <<" completed. \n";
       }

       if( NT % NWP == 0 ) {
           cout << " * Write out at NT = " << NT << endl;
           SWHt = SWHeight( Spc );
           Msgs = " SWHt at all sea points.";
           PrintVector(SWHt, NSea, NT, "wh", Msgs);
           cout << NT << " SWHt[0] = " << SWHt[0] << endl;
       }
   } // End of timestep NT loop.

   // Write a restart file for next step at NSt+NTS.
   Msgs = " Restart file for "+model;
   if( SaveRestart ) RestWrite( Spc, NSt+NTS, "Rst", Msgs ); 

   // Timing at end of program.
   Endstm = std::chrono::steady_clock::now();
   std::chrono::duration<double> elapsed_seconds = Endstm - Startm;
   double t = elapsed_seconds.count(); 
   cout << " Total elpased time in s = " << t << '\n';
   fmsg << " Total elpased time in s = " << t << '\n';
   elapsed_seconds = Endstm - StLoop;
   t = elapsed_seconds.count(); 
   cout << " Loops elpased time in s = " << t << '\n';
   fmsg << " Loops elpased time in s = " << t << '\n';

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

void printarow(int i, float a, float b, float c, float d, 
               std::ostream& ostrm = std::cout) {
    ostrm <<"  "<< i <<"  "<< a <<"  "<< b 
          <<"  "<< c <<"  "<< d << '\n';
}

