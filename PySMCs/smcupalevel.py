"""
;;  Adapted from idl script for PVWave  12 Sept 2005
;;  It reads MFT model output files and plots contours of the
;;  tracer concentration in a ps file.
;;
;; First created: For 4 views     J G Li   26 Nov 2008
;; Last modified: Color grids     J G Li   26 Aug 2009
;; Adapted for 25km SMC grids     J G Li    3 Feb 2010
;; Adapted for 6-25km SMC grids   J G Li   18 Feb 2010
;; Add Arctic boundary rings.     J G Li   10 Jun 2011
;; S4R sterographic projection.   J G Li   30 Jun 2011
;; Modified for SMC625 grid.      J G Li   15 Dec 2011
;; Add 33 buoy locations.         J G Li   11 Apr 2012
;; Adapted for Atlantic model.    J G Li   25 Jul 2012
;; Encapsulated file for report.  J G Li   21 Nov 2012
;; Regional plot around UK.       J G Li   21 Nov 2012
;; Regional plot around Med.      J G Li    4 Jul 2013
;; Regional plot around Med.      J G Li    4 Jul 2013
;; Refined for Euro12 comparison. J G Li   16 Aug 2013
;; SMC6125 for European region.   J G Li   26 Nov 2013
;; SMC6125 for European Arctic.   J G Li   16 Dec 2013
;; Adapted for 4-level SMC36125.  J G Li   25 Feb 2014
;; Use readcell and steromap.     J G Li   19 May 2014
;; Updated /Onecl projection.     J G Li   11 Nov 2014
## Converted into Python.         JGLi02Nov2020
"""

import sys

def smcupalevel(cell_file, skiprows=0 ):
    """
    Read cells from cell_file (if given) or the default 'tempcell.dat'
    Assume no header to be skipped or skip first skiprows (>0) lines.
    Draw the cells as rectangular cells within its ranges.
    Wrap around if nxwrap > 0 is provided, usually NX for global grid.
    """

    import numpy as np

##  Use default cell file if not defined
    if( len(cell_file) < 2 ):
        cellfile = 'tempcels.dat'
    else:
        cellfile = cell_file

    print(" Read cells from ", cellfile)

##  Read cells with s by a procedure
    cels = np.genfromtxt(cellfile, skip_header=skiprows)
    cel  = np.array( cels ).astype(np.int32)
    ncs  = cel.shape[0]
    print(" Total number of cells =", ncs)

##  Cell array output format for each cell.
    fmtcel='%6d %5d %4d %3d %5d'

##  Double all cell i,j,di,dj to up a level.
    for i in range(ncs):
        for j in range(4):
            cel[i,j] = cel[i,j]*2

#           if( cel[i,j] % 2 == 0 ):
#               cel[i,j] = int(cel[i,j]/2)
#           else: 
#               print("Cell[i,j] not converted!", i, j, cel[i,j])

##  Save converted cell array.
    Cell_New = cellfile+"_New"
    np.savetxt(Cell_New, cel, fmt=fmtcel, comments='')

##  End of smc325t650 conversion function.
    return (0)


def main():
    print(" Input argvs =", sys.argv[1:])
    infile=' '
    skprow=0
    nagv = len(sys.argv) 
    if( nagv > 1 ):
        if( len(sys.argv[1]) > 2 ): infile=sys.argv[1] 
    if( nagv > 2 ):
        if( int(sys.argv[2]) > 0 ): skprow=int(sys.argv[2])

    smcupalevel(infile, skiprows=skprow )

if __name__ == '__main__':
    main()

