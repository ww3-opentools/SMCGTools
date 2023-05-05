"""
#;;
#;  Adapted from idl script for PVWave  12 Sept 2005
#;  It reads MFT model output files and plots contours of the
#;  tracer concentration in a ps file.
#;  S4R sterographic projection.  JGLi30Jun2011
#;  SMC625 global part only.      JGLi12Dec2011
#;  SMC625 SWH plot from WW3 text files.  JGLi16Feb2012
#;  SMC625 full grid and Arc/Atl.   JGLi20Jun2012
#;  SMC625 30 Frequency bin SWH plot from SMC625Tx.  JGLi02Oct2012
#;  SMC6125 full grid spectral test.  JGLi08Jan2013
#;  Input yymm from file.                      JGLi23Oct2013
#;  Adapted for G50SMC grid swh.               JGLi19Nov2013
#;  Extended to include Arctic.                JGLi11Feb2014
#;  Simplified to use readproj.                JGLi14Nov2014
#;  Adapted for CMA SMC6125 grid.              JGLi30Jul2018
##  Converted into Python.                     JGLi20Dec2018
##  Adapted for SMC36125 grid swh plot.        JGLi04Jan2019
##  Use PolyCollection to draw swh plot.       JGLi30Jan2019
##  Adapted for SMC61250 grid swh plot.        JGLi22Feb2019
##  Adapted for SMC61250 propagation test.     JGLi26Feb2019
##  Use swhglobl and swhlocal to plot.         JGLi14Mar2019
##  Modified for SMCGTolls PropOMP test.       JGLi19Oct2021
##
"""

def main( ):
##  Import relevant modules and functions

    import numpy as np

    from matplotlib.collections import PolyCollection
    from datetime import datetime, timedelta

    from readtext import readtext
    from steromap import steromap
    from rgbcolor import rgbcolor

##  Path of the cell projection files
    PySMCs='../PySMCs/'
    DatGMC='../DatGMC/'
    OutDat='../OutDat/'
    WrkDir='../tmpfls/'

##  Model name to be used for all plots.
    ModlName='SMC61250'

##  Use own color map and defined depth colors 
    colrfile = PySMCs+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Possible selection of your plot types. 
    gorloc={0:'Global',1:'EuroArc',2:'Pacific'}

##  Prompt selection choices and ask for one input
    print (" \n ", gorloc)
    instr = input("\n *** Please enter your selected number here > ")
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print ("Invalid selection, program terminated.")
        exit()

    print (" Draw SWH plots "+pltype)

##  Choose global or local verts from different files.
    vrfile = DatGMC+'S650Vrts'+pltype[0:4]+'.npz'
    vrtcls = np.load( vrfile )

    if( pltype == 'Global' ):
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']
        svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
        config = vrtcls['cnfg']
        print (' n/svrts/cels config read ' )
        from swhglobl import swhglobl

    else:
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel'] 
        config = vrtcls['cnfg']
        print (' nvrts, ncels and config read ' )
        from swhlocal import swhlocal

##  EuroArc and Pacific paper orientations
    if( pltype == 'Pacific' ):
        papror='landscape'
    else:
        papror='portrait'

#;; Define spectral direction
    ndir=36
    theta=np.arange(ndir)*np.pi*2.0/ndir

#;; Add a spectral array plots for the Northern stripe 
    x0= 3.0
    y0= 5.0
    t0=np.pi*0.25
    cs=np.cos(theta + t0)
    xn=theta*0.0+x0
    yn=theta*0.0+y0
    for i in range(ndir):
        if( cs[i] > 0.0 ): 
            spc=1.2*cs[i]*cs[i]
            xn[i]=x0+spc*np.cos(theta[i])
            yn[i]=y0+spc*np.sin(theta[i])

#;; Add another spectral array plots for the Southern stripe 
    x1=-2.0
    y1=-8.0
    cs=np.cos(theta - t0)
    xs=theta*0.0+x1
    ys=theta*0.0+y1
    for i in range(ndir):
        if( cs[i] > 0.0 ): 
            spc=1.2*cs[i]*cs[i]
            xs[i]=x1+spc*np.cos(theta[i])
            ys[i]=y1+spc*np.sin(theta[i])

#;;  Polar disk spectral array uses the Southern one but new location
    xp= 3.0
    yp= 8.0

#;;  Atlantic disk spectral array uses the Southern one but new location
    xt=-1.0
    yt=-2.0

##  Use ijk to count how many times to draw.
    ijk=0

#;; Specify number of steps per hour (DTG = 900 s)
    nhr=4

#;; Read in cell concentration data files from a list file
    hdr, cnfiles = readtext(OutDat+'cfiles.txt')
    cfiles = cnfiles.astype(np.str).reshape(len(cnfiles))

##  loop over available files 
#   for nn in range(0,len(cnfiles),2):
    for nn in range(0,len(cnfiles),1):
        dfile=OutDat+cfiles[nn] 

        hdlist, swh2d = readtext(dfile)
        mt = int(hdlist[0])
        mc = int(hdlist[1])
        swhs = swh2d.flatten()[0:mc]

        print (' Plotting cell number mc = %d' % mc )

#;; Convert time step for output file
        ntsp='NTS = %5d' % (mt)
        thrs='T = %6.2f hr' % (float(mt)/nhr) 

##  Call function to draw the swh plot.
        figfl = WrkDir + 'Hs' + cfiles[nn][2:7] + '.ps'
        if( pltype == 'Global' ):
            swhglobl(swhs,nvrts,ncels,svrts,scels,colrs,config,
                 mdlname= ModlName, datx=thrs,psfile=figfl)
        else:
            swhlocal(swhs,nvrts,ncels,colrs,config,
                 mdlname= ModlName, datx=thrs,psfile=figfl,
                 paprorn=papror )

##  Increase ijk for next plot
        ijk += 1
        print (" Finish plot No.", ijk," at ", datetime.now())

##  End of date loop

##  End of SMC61250Props.py main program ##


if __name__ == '__main__':
    main()

