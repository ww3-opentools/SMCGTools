"""
## Program to draw scalar tests on SMC1MR3 or SMC1d grids.
## 
## First created:    JGLi30Jun2011
## Converted into Python.    JGLi20Dec2018
## Last modified:    JGLi06Nov2026
##
"""

def main():

## Import relevant modules and functions
    import os
    import sys
    import numpy as np
    import matplotlib.pyplot as plt

    from datetime import datetime
    from readtext import readtext
    from readcell import readcell
    from rgbcolor import rgbcolor
    from whtliner import whtliner 
    from smcfield import smcfield 
    from addtexts import addtexts 

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
        zdlonlat = np.array( nxline, dtype=float )
        print(" Input grid zlon zlat dlon dlat = \n", zdlonlat)
        nxline = flhdl.readline().split()
        Wrkdir=nxline[0]
        DatSWE=nxline[1]
        nxline = flhdl.readline().split()
        NPol = int(nxline[0])
        Arctic = int(nxline[1])
        print(" Number of polar cells NPol = ", NPol)

## Model and grid names.
    if(   Model == 'SMC1d' ): Gname='SMC$\mathregular{1^o}$'
    elif( Model == 'SMCHfd'): Gname='SMC$\mathregular{0.5^o}$'
    elif( Model == 'SMC1R3'): Gname='SMC$\mathregular{1^o}$R3'
    else:                     Gname=Model.strip()

## Read global and Arctic part cells. 
    Cel_file = [ DatSWE+Model+'Cels.dat' ]
    if( Arctic > 0 ):
        Cel_file.append( DatSWE+Model+'BArc.dat' )

    headrs, cel = readcell( Cel_file ) 
    ng = int( headrs[0].split()[0] )
    if( Arctic > 0 ):
        NArB = headrs[1].split()
        na = int( NArB[0] )
        nab = int( NArB[1] )
        ngb = int( NArB[2] )
    else:
        na = 0
        nab = 0
        ngb = 0
    nc = ng + na
    print (' Merged total cel number = %d' % nc )

## Height scale used for plot.
    height0=[0,1,2,3,4,5,6]
    height1=[0.0,0.5,1.0,1.5,2.0]
    height2=[1000,1500,2000,2500,3000,3500]

## Possible selection of your depth scale.
    gorloc={0:'SRotatn 0-6',1:'Deform',2:'SWEs 1-3k'}
    print (" \n ", gorloc)
    indps = input("\n *** Please enter your selected scale number > ")
    k = int(indps)
    dpnmbr=gorloc.get(k, 'Invalid_selection')
    if( dpnmbr == 'Invalid_selection' ): 
        exit()
    if( k == 0 ): heights = height0
    if( k == 1 ): heights = height1
    if( k == 2 ): heights = height2
 
## Font sizes
    fntsz=10.0
    fntsa=1.20*fntsz 
    fntsb=1.50*fntsz

## Use own color map and defined depth colors 
    colrfile = './rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

## Possible selection of your plot types. 
    gorloc={0:'Global',1:'Single',2:'Sglcky',3:'Deform',4:'SRogle'}

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
    vrfile = DatSWE+Model+'Vrts'+pltype[0:4]+'.npz'
    vrtcls = np.load( vrfile, allow_pickle=True )

    if( pltype == 'Global' ):
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']
        svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
        config = vrtcls['cnfg']
        print (' n/svrts/cels config read ' )

    else:
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel'] 
        config = vrtcls['cnfg']
        print (' nvrts, ncels and config read ' )

    papror='portrait'

## Two bottom side text messages positions.
    sztpxy=config[1]
    tpx=sztpxy[2]; tpy=sztpxy[3]+0.5
    xydxdylft=[ tpx, tpy, 0.0, 0.5]
    xydxdyrgt=[-tpx, tpy, 0.0, 0.5]

## Use ijk to count how many times to draw.
    ijk=0

## Floor bottom height needs to be removed for SWEs tests as
## output is Hw + Btm
    Btm= -3000.0 

## Read in cell concentration data files from a list file
    hdr, cnfiles = readtext(Wrkdir+'cfiles.txt')
    cfiles = cnfiles.astype(str).reshape(len(cnfiles))

    print (' last file to draw is '+cfiles[-1] )
## Solid rotation test requires half or full circle timesteps for RMSE.
    if( k == 0 ):    
        instr = input("\n *** Please enter SR test RMSE timestep > ")
        nrmse = int( instr )

## loop over available files 
    for nn in range(0,len(cnfiles)):
        dfile=Wrkdir+cfiles[nn] 

        hdlist, swh2d = readtext(dfile)
        mt = int(hdlist[0])
        mc = int(hdlist[1])
        swhs = swh2d.flatten()[0:mc] 
        if( k == 2 ): swhs = swhs - Btm

## Skip Arctic polar cell if nc = nga
        if( mc != nc ):
            print(' Unmatching mc/nc = %d %d' % (mc, nc)) 
            exit()
        else:
            print(' Plotting cell number mc = %d' % mc )

## Convert time step for output file
        ntsp='NTS = %5d' % (mt)

## Convert water height into linear color indexes.
        nwhts, whtmnx, whtscl = whtliner( swhs, heights )

        txtary=[ ['WHmn='+whtmnx[0], 'b', fntsa],
                 ['WHmx='+whtmnx[1], 'r', fntsa],
                 [ntsp,      'k', fntsb] ] 

        epsfl = Wrkdir+'Hw'+pltype[0:4]+cfiles[nn][2:8]+'.eps'

## Call function to draw the swh plot.
        if( pltype == 'Global' ):

## Draw WH field on northern heimisphere panel.
            fig = plt.figure( figsize=(sztpxy[0:2]) )
            ax1=fig.add_subplot(1,2,1)
            smcfield(ax1, nwhts, nvrts, ncels, colrs, config,
                     vscle=whtscl, vunit='WHt m')

## Draw field on southern hemisphere subplot panel.
            ax2=fig.add_subplot(1,2,2)
            smcfield(ax2, nwhts, svrts, scels, colrs, config,
                     vscle=whtscl, vunit=' ')

## Put statistic information inside subplot ax2
            ax2.text(sztpxy[2], -sztpxy[3]-0.9, Model+' WHt m',
                horizontalalignment='center', fontsize=fntsb, color='r')
            xydxdy=[sztpxy[2], sztpxy[3], 0.0, 0.6]
            addtexts(ax2, xydxdy, txtary)
            plt.subplots_adjust(left=0.01, bottom=0.01, right=0.99, 
                                 top=0.99, wspace=0.01, hspace=0.01)

## All other single panel figures.
        else: 
            fig, ax = plt.subplots(figsize=sztpxy[0:2])

## Draw the water height field on given axis panel.
            smcfield(ax, nwhts, nvrts, ncels, colrs, config, 
                     vscle=whtscl, fntsz=fntsz, vunit='WHt m')

## Add texts on two sides for each panel.
            addtexts(ax, xydxdylft, txtary, hrzlgn='left')
            txtrgt=[ [Gname,   'k', fntsa], 
                     ['WHt m', 'r', fntsa], 
                     ['UNO3',  'b', fntsa] ]
            addtexts(ax, xydxdyrgt, txtrgt, hrzlgn='right')
            plt.subplots_adjust(left=0, bottom=0, right=1, top=1)

## Add rmse values for solid rotation or deformation tests.
            ## Solid roation test save initial field and rmse at half circles.
            if( k == 0 ):
                if( mt == 0 ):
                    swh0=swhs
                elif( mt % nrmse == 0 ):
                    swhdf = np.array(swhs) - np.array(swh0)
                    rmse=np.sqrt((swhdf**2).mean())
                    rmsetx=f'rmse={rmse:7.5f}'
                    ax.text(-tpx,-tpy-0.9, rmsetx,
                       horizontalalignment='right', fontsize=fntsa, color='k')
                    print(" Rotation "+rmsetx )

            ## Deformation test, read exact solution.
            if( k == 1 and mt > 0 ):  
                xfile=Wrkdir+'Cx'+cfiles[nn][2:] 
                hdlist, swh2d = readtext(xfile)
                mxt = int(hdlist[0])
                mxc = int(hdlist[1])
                exhs = swh2d.flatten()[0:mxc] 
                swhdf = np.array(swhs) - np.array(exhs)
                rmse=np.sqrt((swhdf**2).mean())
                rmsetx=f'rmse={rmse:7.5f}' 
                ax.text(-tpx,-tpy-0.9, rmsetx,
                   horizontalalignment='right', fontsize=fntsa, color='k')
                print(" Mod - Exc "+rmsetx )

## Save plot and close the fig object.
        print (" Save the plot "+epsfl )
        plt.savefig(epsfl, dpi=None,facecolor='w',edgecolor='w', \
                    orientation=papror)
        plt.close(fig)

## Increase ijk for next plot
        ijk += 1
        print (" Finish plot No.", ijk," at ", datetime.now())

## End of date loop

## End of main() function.
    
if __name__ == '__main__':
    main()

## End of program ScalarPropC.py.

