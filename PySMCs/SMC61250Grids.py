"""
##
##  It reads cell files and uses own projection to draw 
##  the SMC grid. Projected polygons are collected into 
##  vert variables for grid and subsequent swh plots. 
##
#;  First created: For 4 views     J G Li   26 Nov 2008
#;  Modified for Color grids       J G Li   26 Aug 2009
#;  Adapted for 25km SMC grids     J G Li    3 Feb 2010
#;  Updated G25SMC-12c grids.      J G Li   25 May 2011
#;  Sterographic projection nR.    J G Li   30 Jun 2011
#;  Adapted for G50SMC grid.       J G Li   19 Aug 2011
#;  Extended to fill the Arctic.   J G Li    5 Oct 2011
#;  Rectify polar cell position.   J G Li   25 Oct 2011
#;  Simplify with readcell and steromap.  JGLi12Nov2014
##
##  Converted into a Python function.     JGLi05Dec2018
##  Save ELat/Lon and sx/yc in file.      JGLi11Dec2018
##  Add color map and draw color bar.     JGLi14Dec2018
##  Adapted for SMC36125 grid plot.       JGLi03Jan2019
##  Import resource and set stacksize.    JGLi07Jan2019
##  Use polycollections for two plots.    JGLi30Jan2019
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Updated to Python 3.6.8 GCC 7.3.0.    JGLi02Apr2019
##  Modified for updated SMC61250 grid.   JGLi18May2021
##  Modified for SMCGTools package.       JGLi06Oct2021
##
"""

def main():

##  Import relevant modules and functions
    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from datetime import datetime
    from readcell import readcell   
    from rgbcolor import rgbcolor
    from smcelvrts import smcelvrts

    print( " Program started at %s " % datetime.now().strftime('%F %H:%M:%S') )

##  Read global and Arctic part cells. 
    Wrkdir='../'
    buspct=Wrkdir+'Bathys/SPBuoys.dat'
    Cel_file = Wrkdir+'DatGMC/SMC61250Cels.dat'
    Arc_file = Wrkdir+'DatGMC/SMC61250BArc.dat'

    headrs, cel = readcell( [Cel_file, Arc_file] ) 
    ng = int( headrs[0].split()[0] )
    NArB = headrs[1].split()
    na = int( NArB[0] )
    nb = int( NArB[1] )
    nc = ng + na
    print (' Merged total cel number = %d' % nc )

##  Size-1 cell increments
#   dxlon=0.703125
#   dylat=0.468750
    dxlon=0.087890625
    dylat=0.058593750
#   dxlon=0.0439453125
#   dylat=0.029296875

##  Reference or i=j=0 point lon lat.
    zrlon=0.0
    zrlat=0.0

##  Maximum j row number in Global part
    jmxglb = np.max( cel[0:nc-na,1] )
    print (' Maximum j row =', jmxglb )

##  Grid related parameters
    zdlnlt = [zrlon, zrlat, dxlon, dylat]

##  Extra array in config for Arctic part
    Arctic = True
    ngabjm = [ng, na, nb, jmxglb]

##  Use own color map and defined depth colors 
    colrfile = Wrkdir+'PySMCs/rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Maximum mapping radius.
    radius=10.0

##  Possible selection of your plot types. 
    gorloc={0:'Global',1:'EuroArc',2:'Pacific'}

##  Prompt selection choices and ask for one input
    print (" \n ", gorloc)
    instr = input(' *** Please enter your selected number here > ')
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print ("Invalid selection, program terminated.")
        exit()

    print (" Draw SMC grid "+pltype)

    if( pltype == 'Global'):
##  Whole global projection angle from N Pole to be 90.0
        pangle=90.0
        plon= 0.0 
        plat= 23.5 
        clrbxy=[ -9.6,-12.0, 19.0,  1.0]
        sztpxy=[ 16.0,  9.0,-10.2,-10.3]
        rngsxy=[-10.0, 10.0,-12.1, 10.0]
        from smcglobl import smcglobl

    if( pltype == 'EuroArc'):
##  Euro-Arctic regional plot
        pangle=27.5 
        plon=  0.0
        plat= 69.0
        clrbxy=[  8.0, -4.5,   0.8,  7.0]
        sztpxy=[ 11.0, 15.0,   4.0, -5.4]
        rngsxy=[-11.0, 11.0, -15.0, 15.0]
        papror='portrait'

    if( pltype == 'Pacific'):
##  West Pacific regional plot
        plon= 138.0
        plat=  18.0
        pangle=33.0 
        clrbxy=[-13.6,  7.5,   7.5,  0.8]
        sztpxy=[ 15.0, 10.0, -10.0,  7.0]
        rngsxy=[-15.0, 15.0, -10.0, 10.0]
        papror='landscape'

    print( " Start loop over cells at %s " % datetime.now().strftime('%F %H:%M:%S') )

##  Projection pole lon-lat and angle
    rdpols=[radius, pangle, plon, plat]

##  Create cell vertices
    nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel, zdlnlt, rdpols, rngsxy, NArB=NArB )

##  Set plot size and limits and message out anchor point.
    config=np.array([rdpols, sztpxy, rngsxy, clrbxy, ngabjm])
    pzfile=Wrkdir+'DatGMC/S650Vrts'+pltype[0:4]+'.npz'

##  Store selected north and south verts and cell numbers for swh plots.
##  Use the np.savez to save 3/5 variables in one file.  JGLi22Feb2019 
    if( pltype == 'Global' ):
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config, 
                          svrt=svrts, scel=scels)
    else:
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 

##  These variables could be loaded back by
#   vrtcls = np.load(Wrkdir+'DatGMC/S650Vrts'+pltype[0:4]+'.npz')
#   nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']; config=vrtcls['cnfg']; 
#   svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
##

##  Draw your selected grid plot.
    psfile=Wrkdir+'tmpfls/smc61250grd'+pltype[0:4]+'.ps' 
    bufile=Wrkdir+'Bathys/ECBuoys.dat'

    if( pltype == 'Global'):
        from smcglobl import smcglobl
        smcglobl( cel, nvrts,ncels,svrts,scels,colrs, config,
             mdlname='SMC61250', buoys=buspct, psfile=psfile)

    else:
        from smclocal import smclocal
        smclocal( cel, nvrts,ncels,colrs,config, Arctic=Arctic, 
              mdlname='SMC61250', buoys=bufile, psfile=psfile,
              paprorn=papror)

    print( " Program finished at %s " % datetime.now().strftime('%F %H:%M:%S') )

## End of SMC61250Grids main program ##


if __name__ == '__main__':
    main()

