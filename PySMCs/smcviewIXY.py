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

def smcviewrap(cell_file, nxwrap=0, skiprows=0 ):
    """
    Read cells from cell_file (if given) or the default 'tempcell.dat'
    Assume no header to be skipped or skip first skiprows (>0) lines.
    Draw the cells as rectangular cells within its ranges.
    Wrap around if nxwrap > 0 is provided, usually NX for global grid.
    """

    import numpy as np
    import matplotlib.pyplot as plt
    from rgbcolor import rgbcolor
    from scale_depth import scale_depth

##  Use default cell file if not defined
    if( len(cell_file) < 2 ):
        cellfile = 'tempcels.dat'
    else:
        cellfile = cell_file

    print(" Read cells from ", cellfile)

##  Read cells with s by a procedure
    cels = np.genfromtxt(cellfile, skip_header=skiprows)
    cel  = np.array( cels ).astype(np.int)
    ncs  = cel.shape[0]
    print(" Total number of cells =", ncs)

##  Wrap cel[:,0] values if nxwrap > 0 
    if( nxwrap > 0 ):
        cel[:,0] = np.where(cel[:,0]>nxwrap/2, cel[:,0]-nxwrap, cel[:,0])

##  Range of selected cells
    ixmin=min(cel[:,0])-2
    ixmax=max(cel[:,0])+8
    jymin=min(cel[:,1])-2
    jymax=max(cel[:,1])+8 
    print( ' Cel i, j range =', ixmin, ixmax, jymin, jymax )
    figx=(ixmax-ixmin)/10
    figy=(jymax-jymin)/10
    print( 'Figure size =',figx, figy)

##  Use color to show depth
    colrfile = 'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )
##  Use only first 131 colors in colrs(0:255). 
    depth, factr, cstar, marks, ncstr, nclrm = scale_depth(nclrm=136)
    ndeps = np.zeros( (ncs), dtype=np.int )
    ndep0 = ncstr + np.rint( (cstar-np.log10(11))*factr ).astype(np.int)

    fig=plt.figure(figsize=[figx,figy])
    ax=fig.add_subplot(1,1,1)
    xprop={'xlim':[ixmin,ixmax]} 
    yprop={'ylim':[jymin,jymax]}
    ax.set(**xprop)
    ax.set(**yprop)
    ax.set_aspect('equal')
    ax.set_autoscale_on(False)
    ax.set_axis_off()

    for i in range(ncs - 1):
        xc=[cel[i,0],cel[i,0]+cel[i,2],cel[i,0]+cel[i,2],cel[i,0],cel[i,0]]
        yc=[cel[i,1],cel[i,1],cel[i,3]+cel[i,1],cel[i,3]+cel[i,1],cel[i,1]]
        if( cel[i,4] > -11 ): 
            ndeps[i] = ncstr + np.rint( (cstar-np.log10(cel[i,4]+11))*factr ).astype(np.int)
        if( ndeps[i] == ndep0 ):
            ax.plot(xc, yc, linestyle='-', linewidth=1.0, color='k')
        else:
            ax.plot(xc, yc, linestyle='-', linewidth=1.0, color=colrs(ndeps[i]))
            
        plt.text( (xc[0]+xc[2])*0.5, (yc[0]+yc[2])*0.5, '.', 
            horizontalalignment='center', fontsize=11, color='r' )

##  Readjust axis range in case not exactly as specified.
    ax.set_ylim(bottom=jymin, top=jymax)
    ax.set_xlim(left=ixmin, right=ixmax)
    plt.subplots_adjust(wspace=0, hspace=0)

##  Put i, j index at central cross 
    istar = int(ixmin/8)*8
    jcntr = int((jymin + jymax)/8)*4
    for i in range(istar, ixmax, 8): 
        ii = i
        if( nxwrap > 0 and i < 0): ii = i + nxwrap
        iii=int(ii/8) + 1
        plt.text( i, jcntr, str(iii), rotation=90, 
            horizontalalignment='right', fontsize=11, color='r' )

    jstar = int(jymin/8)*8
    icntr = int((ixmin + ixmax)/8)*4
    for j in range(jstar, jymax, 8):
        jj = int((j+2816)/8) + 1
        plt.text( icntr, j+0.2, str(jj), 
            horizontalalignment='right', fontsize=11, color='r' )

    plt.savefig('./tmpsmcgrd.png', dpi=None,facecolor='w',edgecolor='w',
                 orientation='portrait')
                 

    plt.show()

def main():
    print(sys.argv[1:])
    infile=' '
    nxwrap=0
    skprow=0
    nagv = len(sys.argv) 
    if( nagv > 1 ):
        if( len(sys.argv[1]) > 2 ): infile=sys.argv[1] 
    if( nagv > 2 ):
        if( int(sys.argv[2]) > 0 ): nxwrap=int(sys.argv[2])
    if( nagv > 3 ):
        if( int(sys.argv[3]) > 0 ): skprow=int(sys.argv[3])

    smcviewrap(infile, nxwrap=nxwrap, skiprows=skprow )

if __name__ == '__main__':
    main()

