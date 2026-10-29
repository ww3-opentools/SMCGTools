#ifndef CPP_SMCGRID_H
#define CPP_SMCGRID_H

/***************************************************************************
 *  The SMCGrid.h declares the cell and side varaible type for SMC grid and 
 *  also declares the SMCELL and SMCSID structures.          JGLi26May2026  
 ***************************************************************************/
#include <array>
#include <vector>
#include <string>
#include <iosfwd>

// typedef int celtyp[5];    //Define a new data type for a single cell.
using celtyp = std::array<int, 5>;   // C++17 type definition.

// Declare celtyp i/o format functions to be shared.
std::istream &operator >>(istream &stream, celtyp &one);
std::ostream &operator <<(ostream &stream, celtyp &one);

/* SMC grid side array has 7/8 elements but the J-side last element is only 
 * used for sorting the face array and it is not required in propagation subs. 
 * So sidtyp definition will only use 7 elements for both I and J-side arrays. 
 * J-side last element will be skipped at file reading stage. JGLi26Mar2026 */

// typedef long int sidtyp[7];   //Define side type for a single side.
using sidtyp = std::array<int, 7>;

// Declare sidtyp i/o format functions to be shared.
std::istream &operator >>(istream &stream, sidtyp &one);
std::ostream &operator <<(ostream &stream, sidtyp &one);

/* Declare cell and side structures.  */
struct SMCELL {
    int   nlvl, nglo, ncel, npol;
    float zlon, zlat, dlon, dlat, PCRd;
    bool Arctic;
    std::vector<int>  ncelv, narcs;  
    std::vector<float>  CLat, CAra;  
    std::vector<celtyp>  CEL;     
};

struct SMCSID {
    int  nufc, nuar, nvfc, nvar;
    std::vector<int>  isdlv, jsdlv; 
    std::vector<float>  CLatv; 
    std::vector<sidtyp>  ISD, JSD;     //Declare I/JSD arrays. 
};

// Initial program for struct SMCELL and SMCSID.
void readSMCELL(const std::string& inpmodel, int nlv1, int NCM, int NPol, 
                float zdlnlt[4], bool Arctic = true);
void readSMCSID(const std::string& inpmodel, int nlv1, int NFM, int NPol);
 
// Get read-only SMCELL and SMCSID variables.
const SMCELL& getSMCELL();
const SMCSID& getSMCSID();

/* End of SMCGrid.h  */

#endif

