"""
#;  Adapted from idl script for PVWave  12 Sept 2005
#;  It reads MFT model output files and plots contours of the
#;  tracer concentration in a ps file.
#;  S4R sterographic projection.  JGLi30Jun2011
#;  SMC625 global part only.      JGLi12Dec2011
#;  SMC625 SWH plot from WW3 text files.  JGLi16Feb2012
#;  SMC625 30 Frequency bin SWH plot.     JGLi02Oct2012
#;  SMC6125 grid swh plot from Vn4 SMC6125Tx.  JGLi10Jan2013
#;  Input yymm from file.                      JGLi23Oct2013
#;  Adapted for G50SMC grid swh.               JGLi19Nov2013
#;  Extended to include Arctic.                JGLi11Feb2014
#;  Converted into Python.                     JGLi20Dec2018
##  Adapted for SMC36125 grid swh plot.        JGLi04Jan2019
##  Use PolyCollection to draw swh plot.       JGLi30Jan2019
##  Adapted for SMC61250 grid swh plot.        JGLi22Feb2019
##  Use swhglobl function to draw plot.        JGLi28Feb2019
##  Optional global or local function.         JGLi01Mar2019
##  Add paper orientation and format options.  JGLi06Mar2019
##  Modified for updated SMC61250 grid.        JGLi18May2021
##
"""


def main(): 
##  Import relevant modules and functions

    import numpy as np
    import pandas as pd

    from datetime import datetime, timedelta

    from readtext import readtext
    from readcell import readcell
    from rgbcolor import rgbcolor

##  Path of the cell projection files
    DatGMC='../DatGMC/'
    Wrkdir='../tmpfls/' 
    MCodes='../PySMCs/'
    swhdir='../../S61250Tx/ww3.'

##  Read global and Arctic part cells. 
    Cel_file = DatGMC+'SMC61250Cels.dat'
    Arc_file = DatGMC+'SMC61250BArc.dat'

    headrs, cel = readcell( [Cel_file, Arc_file] )
    ng = int( headrs[0].split()[0] )
    na = int( headrs[1].split()[0] )
    nb = int( headrs[1].split()[1] )
    nc = ng + na 
    print (' Merged total cel number = %d' % nc )

##  Maximum j row number in Global part
    jmxglb = cel[ng-1,1] 
    print (' Maximum j row = %d' % jmxglb )

##  Use own color map and defined depth colors 
    colrfile = MCodes+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Read start and end datetime from fdate
    datefl = open( MCodes+'strendat', 'r')
    strend = datefl.read().split()
    datefl.close()

##  Convert into datetime variables
    start = datetime.strptime(strend[0], '%y%m%d%H')
    endat = datetime.strptime(strend[1], '%y%m%d%H')
    timdx = pd.date_range(start=start, end=endat, freq=strend[2])

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
        print (' n/svrts/cels config read ') 
        from swhglobl import swhglobl

    else:
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel'] 
        config = vrtcls['cnfg']
        print (' nvrts, ncels and config read ') 
        from swhlocal import swhlocal

##  EuroArc and Pacific paper orientations
    if( pltype == 'Pacific' ):
        papror='landscape'
    else:
        papror='portrait'
    
##  Loop over datetime of swh files
    swhdir = '../../S61250Tx/ww3.'

    print (" SWH file loop started at ", datetime.now())

##  Use ijk to count how many times to draw.
    ijk=0

#   for i in range( ndays*4 ):
    for dt in timdx:
        swhfl = swhdir + dt.strftime('%y%m%d%H') + '.hs'

        hdlist, swh2d = readtext(swhfl)
        mc = int(hdlist[4])
        swhs = swh2d.flatten()[0:mc]

##  Skip Arctic polar cell if nc = nga
        if( mc != nc ):
            print ( ' Unmatching mc/nc = %d %d' % (mc, nc) ) 
            exit()
        else:
            print (' Plotting cell number mc = %d' % mc )

##  Convert time step for output file
        datms = dt.strftime('%Y%m%d%H')

##  Call function to draw the swh plot.
        figfl = Wrkdir + 'swh' + pltype[0:4]+dt.strftime('%y%m%d%H') + '.ps'
        if( pltype == 'Global' ):
            swhglobl(swhs,nvrts,ncels,svrts,scels,colrs,config,
                 mdlname='SMC61250',datx=datms,psfile=figfl)
        else:
            swhlocal(swhs,nvrts,ncels,colrs,config,
                 mdlname='SMC61250',datx=datms,psfile=figfl,
                 paprorn=papror )

##  Increase ijk for next plot
        ijk += 1
        print (" Finish plot No.", ijk," at ", datetime.now())

##  End of date loop

##  End of SMC61250SWHts.py main function.

if __name__ == '__main__':
    main()

