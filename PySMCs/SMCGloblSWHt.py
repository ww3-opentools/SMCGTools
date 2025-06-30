"""
## Program to draw SWH field on global SMC grid.
##
## First created:    JGLi30Jun2011
## Converted into Python.    JGLi20Dec2018
## Last modified:    JGLi06Jun2025
##
"""

def main(): 

## Import relevant modules and functions
    import sys
    import numpy as np
    import pandas as pd
    import matplotlib.pyplot as plt

    from datetime import datetime, timedelta

    from readtext import readtext
    from readcell import readcell
    from rgbcolor import rgbcolor
    from smcswhcv import smcswhcv 
    from smcfield import smcfield 
    from addtexts import addtexts 

## Check input information file name if provided.
    print(sys.argv)
    if( len(sys.argv) > 1 ):
        if( len(sys.argv[1]) > 3 ):
            gridfile = sys.argv[1]
    else:
        gridfile = 'GridInfo61250.txt'

## Read global grid information file. 
    with open( gridfile, 'r' ) as flhdl:
## First line contains grid name and number of resolution levels.
        nxlne = flhdl.readline().split()
        Gname = nxlne[0]
        Level = int(nxlne[1])
        print(" Input grid name and number of levl= ", Gname, Level)
## Second line contains zlon zlat dlon dlat of size-1 cell parameters.
        nxlne = flhdl.readline().split()
        zdlnlt = np.array(nxlne, dtype=float)
        print(" Input grid zlon zlat dlon dlat = \n", zdlnlt) 
## Third line is the working directory and cell array subdirectory.
        nxlne = flhdl.readline().split()
        Wrkdir=nxlne[0]
        DatGMC=nxlne[1]
        print(" Working directory and DatGMC = \n", nxlne) 
## Fourth line starts with the number of polar cells.
        nxlne = flhdl.readline().split()
        npl = int(nxlne[0])
        print(" Number of polar cells = ", npl)
## Fifty line is the SWH files and propagation test output directories.
        nxlne = flhdl.readline().split()
        SWHdir=nxlne[0]
        OutDat=nxlne[1]
        print(" SWH and Prop OutDat = \n", nxlne)

## Read the cell array and Arctic part if any.
    Cel_file = [ DatGMC+Gname+'Cels.dat' ]
    if( npl > 0 ):
        Arc_file = DatGMC+Gname+'BArc.dat'
        Cel_file.append( Arc_file )

    headrs, cel = readcell( Cel_file ) 
    ng = int( headrs[0].split()[0] )
    nc = ng
    if( npl > 0 ):
        NArB = headrs[1].split()
        na = int( NArB[0] )
        nb = int( NArB[1] )
        nbg= int( NArB[2] )
        nc = ng + na 
    print (' Merged total cel number = %d' % nc )

## Default Python library path of MyCodes.
    MyCodes='./'

## Use own color map and defined depth colors 
    colrfile = MyCodes+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

## Read start and end datetime from Wrkdir/strendat.
    datefl = open( MyCodes+'strendat', 'r')
    strend = datefl.read().split()
    datefl.close()

## Convert into datetime variables
    start = datetime.strptime(strend[0], '%y%m%d%H')
    endat = datetime.strptime(strend[1], '%y%m%d%H')
    timdx = pd.date_range(start=start, end=endat, freq=strend[2])

## Possible selection of your plot types. 
    gorloc={0:'Global',1:'EuroArc',2:'Pacific',3:'Regular'}

## Prompt selection choices and ask for one input
    print (" \n ", gorloc)
    instr = input("\n *** Please enter your selected number here > ")
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print ("Invalid selection, program terminated.")
        exit()

    print (" Draw SWH plots "+pltype)

## Choose global or local verts from different files.
    vrfile = DatGMC+Gname+'Vrts'+pltype[0:4]+'.npz'
    vrtcls = np.load( vrfile )

    if( pltype == 'Global' ):
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']
        svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
        config = vrtcls['cnfg']
        print (' n/svrts/cels config read ') 

    else:
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel'] 
        config = vrtcls['cnfg']
        print (' nvrts, ncels and config read ') 

    if( pltype == 'Regular' ):
        from smcfrglr import smcfrglr

## Selected plot configuration parameters.
    sztpxy = config[1]
    rngsxy = config[2]

## Alternative font sizes.
    fntsz=12.0
    fntsa=1.20*fntsz 
    fntsb=1.50*fntsz

## Loop over datetime of swh files
    print (" SWH file loop started at ", datetime.now())

## Use ijk to count how many times to draw.
    ijk=0
    for dt in timdx:
        swhfl = SWHdir + dt.strftime('%y%m%d%H') + '.hs'

        hdlist, swh2d = readtext(swhfl)
        mc = int(hdlist[4])
        swhs = swh2d.flatten()[0:mc]

## Skip Arctic polar cell if nc = nga
        if( mc != nc ):
            print (' Unmatching mc/nc = %d %d' % (mc, nc)) 
            exit()
        else:
            print (' Plotting cell number mc = %d' % mc )

## Convert time step for output file
        datms = dt.strftime('%Y%m%d%H')

## Convert swh field into color indexes.
        nswh, swhmnx, swhscl = smcswhcv( swhs )

        txtary=[ [Gname+' SWH',    'k', fntsb],
                 ['SWHmn='+swhmnx[0], 'b', fntsa],
                 ['SWHmx='+swhmnx[1]+' m', 'r', fntsa],
                 [datms,     'k', fntsb] ] 

## Call function to draw the swh plot.
        epsfl = Wrkdir + 'swh' + pltype[0:4]+dt.strftime('%y%m%d%H') + '.eps'
        if( pltype == 'Global' ):
## Set up first subplot and axis for northern hemisphere
            fig=plt.figure(figsize=sztpxy[0:2])
            ax1=fig.add_subplot(1,2,1)

## Draw SWH field on northern heimisphere panel.
            smcfield(ax1, nswh, nvrts, ncels, colrs, config,
                     vscle=swhscl, vunit='SWH m')

## Draw field on southern hemisphere subplot panel.
            ax2=fig.add_subplot(1,2,2)
            smcfield(ax2, nswh, svrts, scels, colrs, config,
                     vscle=swhscl, vunit=' ')

## Put statistic information inside subplot ax2
            ax2.text(sztpxy[2], 9.0, txtary[0][0], color='r', 
                horizontalalignment='center', fontsize=fntsb)
            xydxdy=[sztpxy[2], sztpxy[3], 0.0, 0.6]
            addtexts(ax2, xydxdy, txtary[1:])

            plt.subplots_adjust(left=0.01, bottom=0.0, right=0.99, 
                                 top=1.0, wspace=0.01, hspace=0.0)

        elif( pltype == 'Regular' ):
## Draw field as regular grid box format.
            fig, ax = plt.subplots(figsize=sztpxy[0:2])

            smcfrglr(ax, nswh, nvrts, ncels, colrs, config,
                     vscle=swhscl, vunit='SWH m') 

            tpx=(rngsxy[1] - rngsxy[0])*0.5
            tpy= rngsxy[2]+1.0
            ax.text(tpx, tpy, Gname+' SWHmn='+swhmnx[0]+' SWHmx='+swhmnx[1]+' at '+datms,
                horizontalalignment='center', fontsize=fntsa, color='k')
            plt.subplots_adjust(left=0.03, bottom=0.04, \
                           right=0.99, top=0.98)

        else:
## Draw other regional plots.
            fig, ax = plt.subplots(figsize=sztpxy[0:2])

            smcfield(ax, nswh, nvrts, ncels, colrs, config,
                     vscle=swhscl, vunit='SWH m')

## Put statistic information inside plot ax.
            xydxdy=[sztpxy[2], sztpxy[3]-0.5, 0.0, -0.6]
            addtexts(ax, xydxdy, txtary)
            plt.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)

## Save plot and close figure plt.
        print(" ... Saving plot as", epsfl )
        plt.savefig(epsfl, dpi=None,facecolor='w',edgecolor='w', 
            orientation='portrait')
        plt.close()

## Increase ijk for next plot
        ijk += 1
        print (" Finish plot No.", ijk," at ", datetime.now())

## End of date loop

## End of main() function.

if __name__ == '__main__':
    main()

## End of SMCGloblSWHt.py program. 

