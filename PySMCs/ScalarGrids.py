"""
##
## Program to draw scalar tests SMC grids, inclding SMC1d, SMCHfd and SMC1R3. 
##
## First created:    JGLi26Nov2008
## Updated to Python 3.    JGLi02Apr2019
## Last modified:    JGLi11Nov2026
##
"""

def main():

## Import relevant modules and functions
    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    from datetime import datetime
    from readcell import readcell   
    from steromap import steromap
    from rgbcolor import rgbcolor
    from smcgrids import smcgrids 
    from addtexts import addtexts 
    from smcelvrts import smcelvrts

    print( " Program started at %s " % datetime.now().strftime('%F %H:%M:%S') )

## Find the $DATADIR directory on the system.
    datadir = os.popen('echo $DATADIR').read()[0:-1]
    homedir = os.popen('echo $HOME').read()[0:-1]

## Check input information file name if provided.
    print(sys.argv)
    if( len(sys.argv) > 1 ):
        if( len(sys.argv[1]) > 3 ):
            gridfile = sys.argv[1]
    else:
        gridfile = 'GridInfoSMC1R3.txt'

## Grid information from file.
    with open( gridfile, 'r' ) as flhdl:
        nxline = flhdl.readline().split()
        Model = nxline[0]
        NLevl = int(nxline[1])
        print(" Input grid name and NLevl= ", Model, NLevl)
        nxline = flhdl.readline().split()
        zdlnlt = np.array( nxline, dtype=float )
        print(" Input grid zlon zlat dlon dlat = \n", zdlnlt)
        nxline = flhdl.readline().split()
        Wrkdir=nxline[0]
        DatSWE=nxline[1]
        nxline = flhdl.readline().split()
        NPol = int(nxline[0])
        Arctic = int(nxline[1])
        print(" Numbers of polar cells NPol & Arctic = ", NPol, Arctic)

## Read global and Arctic part cells. 
    Cel_file = [ DatSWE+Model+'Cels.dat' ]
    if( Arctic > 0 ):
        Cel_file.append( DatSWE+Model+'BArc.dat' )

    headrs, cel = readcell( Cel_file ) 
    ng = int( headrs[0].split()[0] )
    if( Arctic > 0 ):
        NArB = headrs[1].split()
        na  = int( NArB[0] )
        nab = int( NArB[1] )
        ngb = int( NArB[2] )
    else:
        na = 0
        nab = 0
        ngb = 0
    nc = ng + na
    print (' Merged total cel number = %d' % nc )
    nabgpl = [na, nab, ngb, NPol]

## Initial txtary list as [txt, colr, fontsize].
    fntsz = 10.0
    fntsa = 1.25*fntsz
    fntsb = 1.50*fntsz

    if(   Model == 'SMC1d' ): Gname='SMC$\mathregular{1^o}$'
    elif( Model == 'SMCHfd'): Gname='SMC$\mathregular{0.5^o}$'
    elif( Model == 'SMC1R3'): Gname='SMC$\mathregular{1^o}$R3'
    else:                     Gname=Model.strip()

    txtary = [ [Gname+' Grid',   'k', fntsb],
               ['NC='+str(nc),   'r', fntsa], 
               ['NPl='+str(NPol),'r', fntsa] ] 

## Reset cell depth to 4km for dark color grid plot.
    cel[:,4]=4000

## Mark all cells on the Equator and 2 polar cells
    cel1=list(cel[:,1])
    ind0=[]
    for j in range(ng):
        if(cel1[j] == 0): ind0.append(j)
    excids= ind0 + [nc-2,nc-1]
## But remove the first cell at zero meridian.
    excids.remove(ind0[0])

## Use own color map and defined depth colors 
    colrfile = './rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

## Maximum mapping radius.
    radius=10.0

## Possible selection of your plot types. 
    gorloc={0:'Global', 1:'Single', 2:'Sglcky', \
            3:'Deform', 4:'Regulr', 5:'RegNth',6:'SRogle'}

## Prompt selection choices and ask for one input
    print (" \n ", gorloc)
    instr = input(' *** Please enter your selected number here > ')
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print ("Invalid selection, program terminated.")
        exit()

    print (" Draw SMC grid "+Model+pltype)

    if( pltype == 'Global'):
## Whole global projection angle from N Pole to be 90.0
        Pnrds=10.0
        pangle=90.0
        plon=  0.0 
        plat= 40.0 
        clrbxy=[ -9.6,-12.0, 19.0,  1.0]
        sztpxy=[ 16.0, 9.0, -10.2, -9.8]
        rngsxy=[-10.0, 10.0,-12.1, 10.0]

    if( pltype == 'Single'):
##  Singl hemispherical plot viewed at zero-meridion.
        Pnrds=4.0
        pangle=90.0
        plon=  0.0 
        plat= 40.0 
        clrbxy=[ 10.1, -9.6, 1.0, 19.0]
        sztpxy=[  9.0,  8.0, -9.6, -9.9]
        rngsxy=[-10.1, 12.1,-10.1, 10.1]

    if( pltype == 'Sglcky'):
##  Same as Single but color key at lower side.
        Pnrds=4.0
        pangle=90.0
        plon=  0.0 
        plat= 40.0 
        clrbxy=[ -9.6,-12.0, 19.0, 1.0]
        sztpxy=[  8.0,  9.0, -9.6, -9.9]
        rngsxy=[-10.1, 10.1,-12.1, 10.1]

    if( pltype == 'SRogle'):
##  Singl hemispherical plot for solid rotation test.
        Pnrds=10.0
        pangle=90.0
        plon=-60.0 
        plat= 45.0 
        clrbxy=[ -9.6,-12.0, 19.0,  1.0]
        sztpxy=[  8.0,  8.8, -9.6, -9.9]
        rngsxy=[-10.1, 10.1,-12.1, 10.1]

    if( pltype == 'Deform'):
## Same as 'Single' hemisphere but color key at bottom.
        Pnrds=100.0
        pangle=90.0 
        plon=  1.5
        plat= 75.0
        clrbxy=[ -9.6,-12.0, 19.0,  1.0]
        sztpxy=[  8.0,  8.8, -9.6, -9.9]
        rngsxy=[-10.1, 10.1,-12.1, 10.1]

    if( pltype == 'Regulr'):
## Global rectangular regular lat-lon view.
        pangle=90.0 
        plon=-180.0
        plat=  90.0
        clrbxy=[366.0,-86.0, 10.0, 170.0]
        sztpxy=[ 16.0,  8.0, 60.0, 30.0]
        rngsxy=[ 0.0, 390.0,-90.0, 90.0]

    if( pltype == 'RegNth'):
## Rectangular regular lat-lon view for North hemisphere.
        pangle=90.0
        plon=-180.0 
        plat=  90.0
        clrbxy=[366.0,  3.0, 10.0, 85.0]
        sztpxy=[ 16.0,  6.0, 60.0, 30.0]
        rngsxy=[ 0.0, 390.0,  0.0, 90.0]

    print( " Start loop over cells at %s " % datetime.now().strftime('%F %H:%M:%S') )

    papror='portrait'

## Projection pole lon-lat and angle
    rdpols=[radius, pangle, plon, plat]

## Create cell vertices
    if( pltype == 'Global' or pltype == 'Single' or pltype == 'Sglcky' 
                           or pltype == 'Deform' or pltype == 'SRogle'):
## Create cell vertices
        from smcelvrts import smcelvrts
        nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel, zdlnlt, rdpols, 
               rngsxy, NArB=nabgpl, excids=excids, Pnrds=Pnrds )

        nabmrk = [na, nab, nsmrk[0], nsmrk[1]]
        config=np.array([rdpols, sztpxy, rngsxy, clrbxy, nabmrk])

    if( pltype == 'Regulr' or pltype == 'RegNth'):
        from smcelvrglr import smcelvrglr
        nvrts, ncels, nmark = smcelvrglr( cel, zdlnlt, rngsxy, 
                  excids=excids, NArB=nabgpl )

## For regular plot there is no need to store the projection rdpols but zdlnlt.
        config=np.array([zdlnlt, sztpxy, rngsxy, clrbxy, nabgpl])

## Store selected north and south verts and cell numbers for swh plots.
## Use the np.savez to save 3/5 variables in one file.  JGLi22Feb2019 
    pzfile=DatSWE+Model+'Vrts'+pltype[0:4]+'.npz'
    if( pltype == 'Global' ):
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config, 
                          svrt=svrts, scel=scels)
    else:
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 

## These variables could be loaded back by
#   vrtcls = np.load(DatSWE+Model+'Vrts'+pltype[0:4]+'.npz')
#   nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']; config=vrtcls['cnfg']; 
#   svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
##

## Draw your selected grid plot.
    epsfile=Wrkdir+Model+pltype[0:4]+'grd.eps' 

    if( pltype == 'Global'):
        fig=plt.figure(figsize=sztpxy[0:2])
        ax1=fig.add_subplot(1,2,1)

        smcgrids(ax1, cel, nvrts,ncels,colrs, config, Arctic=Arctic, 
                 fntsz=fntsz, nmark=nsmrk[0]) 

        ax2=fig.add_subplot(1,2,2)
        smcgrids(ax2, cel, svrts,scels,colrs, config, Arctic=Arctic, 
                 fntsz=fntsz, nmark=nsmrk[1]) 

        xydxdy=[sztpxy[2], sztpxy[3], 0.0, 0.6]
        addtexts(ax2, xydxdy, txtary)

        plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, 
                             top=0.99, wspace=0.01, hspace=0.01)
        plt.savefig(epsfile, dpi=None,facecolor='w',edgecolor='w', \
                    orientation=papror)
        plt.close()

    elif( pltype == 'Single' or pltype == 'Sglcky' or pltype == 'Deform'
                             or pltype == 'SRogle'):
        fig, ax = plt.subplots(figsize=sztpxy[0:2])
        plt.subplots_adjust(left=0.0, bottom=0.0, right=1.0, top=1.0)

        smcgrids(ax, cel, nvrts,ncels,colrs,config, Arctic=Arctic, 
                 fntsz=fntsz, nmark=nsmrk[0], Pnrds=Pnrds)

        xydxdy=[-9.6, -9.6, 0.0, 0.6]
        addtexts(ax, xydxdy, txtary, hrzlgn='left')
        xydxdy=[ 9.6, -9.6, 0.0, 0.6]
        proary = [ ['Projection Pole',   'k', fntsb],
                   [f'PLon={plon:8.2f}', 'r', fntsa], 
                   [f'PLat={plat:7.2f}', 'r', fntsa] ] 
        addtexts(ax, xydxdy, proary, hrzlgn='right')

        print (" Save the grid plot ... " )
        plt.savefig(epsfile, dpi=None, facecolor='w', \
            edgecolor='w', orientation=papror)
        plt.close()

    else:
        fig, ax = plt.subplots(figsize=sztpxy[0:2])
        plt.subplots_adjust(left=0.03, bottom=0.04, \
                           right=0.99, top=0.98)
        from smcgrglr import smcgrglr

        smcgrglr(ax, cel, nvrts, ncels, colrs, config, 
                 fntsz=fntsz, nmark=nmark)

        tpx=(rngsxy[1] - rngsxy[0])*0.5
        tpy=rngsxy[2]+1.0
        ax.text(tpx, tpy, Gname+' NC ='+str(nc),  
            horizontalalignment='center', fontsize=fntsb, color='k')

        print (" Save the grid plot ... " )
        plt.savefig(epsfile, dpi=None, facecolor='w', \
            edgecolor='w', orientation=papror)

    print(" Program finished at %s " % \
           datetime.now().strftime('%F %H:%M:%S') )

## End of main:() function.

if __name__ == '__main__':
    main()

## End of ScalarGrids program. 

