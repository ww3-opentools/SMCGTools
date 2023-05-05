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
##  Updated for Python 3.                      JGLi28Mar2019
##  Adapted for Pacf36125 sub-grid model.      JGLi06Oct2020
##  Suspend Arctic part from Atln36125.        JGLi18Nov2020
##  Restore Arctic part from Atln36125.        JGLi08Jan2021
##  Adapted for 3 sub-grids wht plots.         JGLi25Jan2021
##  Adapted for Sub61250 grids swht plots.     JGLi11Oct2021
##
"""

if( 1 > 0 ):
##  Import relevant modules and functions

    import numpy as np
    import pandas as pd

    from datetime import datetime, timedelta

    from readtext import readtext
    from readcell import readcell
    from rgbcolor import rgbcolor

##  Path of the cell projection files
    DatSub='../DatSub/'
    WrkDir='../tmpfls/' 

##  Use own color map and defined depth colors 
    colrfile = 'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Read start and end datetime from fdate
    datefl = open( 'strendat', 'r')
    strend = datefl.read().split()
    datefl.close()

##  Convert into datetime variables
    start = datetime.strptime(strend[0], '%y%m%d%H')
    endat = datetime.strptime(strend[1], '%y%m%d%H')
    timdx = pd.date_range(start=start, end=endat, freq=strend[2])

##  Possible selection of your plot types. 
#   gorloc={0:'SubG',1:'Pacf',2:'Atln'}
    gorloc={0:'SubG',1:'Soth',2:'Pacf',3:'Atln'}

##  Prompt selection choices and ask for one input
    print (" \n ", gorloc)
    instr = input("\n *** Please enter your selected number here > ")
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print ("Invalid selection, program terminated.")
        exit()

##  Model name to be used for all plots.
    ModlName=pltype+'61250'

    print (" Draw SWH plots for "+ModlName)

##  Choose global or local verts from different files.
    vrfile = DatSub+ModlName+'Vrts.npz'
    vrtcls = np.load( vrfile )

    if( pltype == 'SubG' ): 
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']
        svrts = vrtcls['svrt'] ; scels = vrtcls['scel']
        config = vrtcls['cnfg']
        print (' n/svrts/cels config read ') 
        from swhglobl import swhglobl
        papror='landscape'

    else:
        nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel'] 
        config = vrtcls['cnfg']
        print (' nvrts, ncels and config read ') 
        from swhlocal import swhlocal
        papror='portrait'

    
##  Loop over datetime of swh files
    swhSoth = '../../Sub650Soth/ww3.'
    swhPacf = '../../Sub650Pacf/ww3.'
    swhAtln = '../../Sub650Atln/ww3.'

    print (" SWH file loop started at ", datetime.now() )

##  Use ijk to count how many times to draw.
    ijk=0

#   for i in range( ndays*4 ):
    for dt in timdx:
##  Read Soth61250 swh.
        swhfl = swhSoth + dt.strftime('%y%m%d%H') + '.hs'

        hdlist, swh2d = readtext(swhfl)
        ms = int(hdlist[4])
        swh1 = swh2d.flatten()[0:ms]

##  Read Pacf61250 swh.
        swhfl = swhPacf + dt.strftime('%y%m%d%H') + '.hs'

        hdlist, swh2d = readtext(swhfl)
        mc = int(hdlist[4])
        swh2 = swh2d.flatten()[0:mc]

##  Read or append Atln51250 swh.
        swhfl = swhAtln + dt.strftime('%y%m%d%H') + '.hs'

        hdlist, swh2d = readtext(swhfl)
        mt = int(hdlist[4])
        swh3 = swh2d.flatten()[0:mt]

##  Convert time step for output file
        datms = dt.strftime('%Y%m%d%H')

##  Call function to draw the swh plot.
        figfl = WrkDir + 'swh' + pltype + dt.strftime('%y%m%d%H') + '.ps'

        if( pltype == 'SubG' ): 
##  Merge all sub-grid SWHs together.
            swhs = np.hstack( (swh1, swh2, swh3) )

            swhglobl(swhs,nvrts,ncels,svrts,scels,colrs,config,
                 mdlname= ModlName,datx=datms,psfile=figfl)
        else:
            if( m == 1 ): swhs = swh1
            if( m == 2 ): swhs = swh2
            if( m == 3 ): swhs = swh3

            swhlocal(swhs,nvrts,ncels,colrs,config,
                 mdlname= ModlName, datx=datms,psfile=figfl,
                 paprorn=papror )

##  Increase ijk for next plot
        ijk += 1
        print (" Finish plot No.", ijk," at ", datetime.now())

##  End of date loop

