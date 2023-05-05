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
##  Adapted for SMC61250 global grid.     JGLi18Feb2019
##  Updated to Python 3.6.8 GCC 7.3.0.    JGLi02Apr2019
##  Adapted for UK12H rotated SMC grid.   JGLi20Jun2019
##  Adapted for SMC Med36125 grid.        JGLi08Jul2019
##  Adapted for SMC Indian Ocean grid.    JGLi18Nov2019
##  Adapted for SMC36125 sub grids.       JGLi01Sep2020
##  Add boundary cell processing.         JGLi03Sep2020
##  Suspend Arctic part in Atln36125.     JGLi18Nov2020
##  Restore Arctic part in Atln36125.     JGLi08Jan2021
##  Draw 3 sub-grids split from SMC36125. JGLi12Jan2021
##  Draw 3 sub-grids split from SMC61250. JGLi08Oct2021
##  Modified to use smcelvrts function.   JGLi26Oct2021
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
    from smcglobl import smcglobl 
    from smclocal import smclocal 
    from smcelvrts import smcelvrts
    from smcellmap import smcell, smcmap, smcids

    print( " Program started at %s " % datetime.now().strftime('%F %H:%M:%S') )

##  Read global and Arctic part cells. 
    DatSub='../DatSub/'
    Wrkdir='../tmpfls/'
    Arcell='../DatGMC/SMC61250BArc.dat'

##  Size-1 cell increments
#   dxlon=0.703125
#   dylat=0.468750
    dxlon=0.087890625
    dylat=0.058593750
#   dxlon=0.0439453125
#   dylat=0.029296875
    zrlon=0.0
    zrlat=0.0
    zdlnlt = [zrlon, zrlat, dxlon, dylat]

##  Use own color map and defined depth colors 
    colrfile = 'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Maximum mapping radius.
    radius=10.0

##  Possible selection of your plot types. 
    gorloc={0:'SubG',1:'Soth',2:'Pacf',3:'Atln'}

##  Prompt selection choices and ask for one input
    print (" \n ", gorloc)
    instr = input(' *** Please enter your selected number here > ')
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print ("Invalid selection, program terminated.")
        exit()

    print( " Start draw "+pltype+" at %s " % datetime.now().strftime('%F %H:%M:%S') )

    if( pltype == 'Soth' ):
##  Whole global projection angle from N Pole to be 90.0
        pangle= 90.0
        plon=-136.0
        plat= 21.0
        clrbxy=[ -9.0,  3.6,   7.0,  0.7]
        sztpxy=[ 12.0, 11.0,  -7.0,  8.0]
        rngsxy=[-10.0, 10.0,  -9.5,  9.0]
        papror='portrait'

    if( pltype == 'Pacf' ):
##  Whole global projection angle from N Pole to be 90.0
        pangle=90.0
        plon=  10.0 
        plat= 22.6 
        clrbxy=[ -9.3, -9.0,  9.9,  0.9]
        sztpxy=[ 11.0, 10.5, -7.5, -2.5]
        rngsxy=[-10.1, 10.1, -9.6, 10.1]
        papror='portrait'
        rdpols=[radius, pangle, plon, plat]

##  Atlantic regional plot for Atln61250 sub-grid
    if( pltype == 'Atln'):
        pangle= 53.0 
        plon=-50.0
        plat= 50.0
        clrbxy=[ -9.0,  1.1,   0.7,  7.0]
        sztpxy=[ 10.0, 12.0,  -7.0, -7.5]
        rngsxy=[-10.0, 10.0, -12.0, 12.0]
        papror='portrait'

##  Full global plot of all sub-grids.
    if( pltype == 'SubG' ):
        pangle=90.0
        plon= 0.0
        plat= 23.5
        clrbxy=[ -9.6,-12.6, 19.0,  1.0]
        sztpxy=[ 16.0, 10.0,-10.2,-10.3]
        rngsxy=[-10.0, 10.0,-13.0, 10.0]

##  Radius, projection angle, and pole lon-lat.
    rdpols=[radius, pangle, plon, plat]

##  Model name for selected grid. 
    ModlName=pltype+'61250'

    if( m > 0 ):
##  Process sub-grid and its boundary cells
        Cellfile = [DatSub+ModlName+'Cels.dat']
        Bndyfile = [DatSub+ModlName+'Bdys.dat']
##  Append Arctic part for Atln sub-grid.
        if( m == 3 ):
            Cellfile.append(Arcell)
 
##  Read sub-grid cell
        headrs, cel = readcell( Cellfile )
        nc = int( headrs[0].split()[0] )
        NArB = []
        if( m == 3 ):
            NArB = headrs[1].split()
            na = int( NArB[0] )
            nb = int( NArB[1] )
            ng = nc 
            nc = ng + na 

##  Read sub-grid boudnary cell
        headrs, bcel = readcell( Bndyfile ) 
        nbdy = int( headrs[0].split()[0] )

##  Find boundary cell ids for the sub-grid.
        ncbdy = smcids(bcel, cel)

##  Find boundary cell central xlon ylat. 
        xlon, ylat = smcell(ncbdy, cel, zdlnlt)
        
##  Generate cell vertices for the sub-grid
        nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel, 
               zdlnlt, rdpols, rngsxy, excids=ncbdy, NArB=NArB )

##  Save boundary cell sequential number list for Pacf61250 grid.
        fmt = '%8d '
        Bndylist = Wrkdir+ModlName+'Blst.dat'
##  Boundary cell sequential numbers have to increase by 1 for WW3. JGLi06Oct2020
        np.savetxt(Bndylist, np.array(ncbdy)+1, fmt=fmt, header='',  comments='')

##  Save the boundary cell xlon, ylat to generate boundary condition in WW3.
        hdr = f'{nbdy:8d} \n'
        fms = '%s '
        Bndylnlt = Wrkdir+ModlName+'Blnlt.dat'
        with open(Bndylnlt, 'w') as flhd:
            flhd.writelines(hdr)
            for j in range(nbdy):
                if( ncbdy[j] >= 0 ): 
                    flhd.write(f'{xlon[j]:9.3f} {ylat[j]:8.3f}   0.0  0.0  1 \n' )
        print(" Boundary cells saved in "+Bndylnlt )

##  For Soth and Pacf sub-grids draw the southern hemisphere
        if( pltype == 'Soth' or pltype == 'Pacf' ):
            Arctic= False
            np000 = [nc, 0, nsmrk[0], nsmrk[1]]
            config=np.array([rdpols, sztpxy, rngsxy, clrbxy, np000])
            pzfile=DatSub+ModlName+'Vrts.npz'
            np.savez( pzfile, nvrt=svrts, ncel=scels, cnfg=config) 
            print(" Verts data saved in "+pzfile )

            psfile=Wrkdir+ModlName+'grd.ps'
            smclocal(cel,svrts,scels,colrs,config,Arctic=Arctic,
                     mdlname=ModlName, psfile=psfile,
                     paprorn=papror, nmark=nsmrk[1])

##  For Atln sub-grid draw the northern hemisphere
        if( pltype == 'Atln' ):
            Arctic = True 
            jmxglb = np.max( cel[0:ng,1] )
            print (' Maximum j row =', jmxglb)
            ngabjm = [ng, na, nb, jmxglb]
            config=np.array([rdpols, sztpxy, rngsxy, clrbxy, ngabjm])
            pzfile=DatSub+ModlName+'Vrts.npz'
            np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 
            print(" Verts data saved in "+pzfile )

            psfile=Wrkdir+ModlName+'grd.ps' 
            smclocal( cel,nvrts,ncels,colrs,config,Arctic=Arctic,
                      mdlname=ModlName, psfile=psfile, 
                      paprorn=papror, nmark=nsmrk[0] )


##  For m = 0 or full global grid.
    if( pltype == 'SubG' ):
##  Read all sub-grid cells plus Arctic part.
        Subs=['Soth','Pacf','Atln']

        Celfiles = [ DatSub+Subs[i]+'61250Cels.dat' for i in [0, 1, 2] ]
        Celfiles.append(Arcell)
        Arctic = True
 
        headrs, cel = readcell( Celfiles )

        ncls = [ int( headrs[i].split()[0] ) for i in [0, 1, 2, 3] ]
        ncs = np.array(ncls).sum()
        print (' Summed total cell number ncs =', ncs)

        nsh = int( headrs[0].split()[0] )
        npc = int( headrs[1].split()[0] )
        ntl = int( headrs[2].split()[0] )
        ng = nsh + npc + ntl
        print(" nsh npc ntl ng =", nsh, npc, ntl, ng )

        NArB = headrs[3].split()
        na = int( NArB[0] )
        nb = int( NArB[1] )
        nc = ng + na
        jmxglb = np.max( cel[0:nc-na,1] )
        print (' Maximum j row =', jmxglb)
        print(" nc ng na nb jmxglb=", nc, ng, na, nb, jmxglb)

##  Read sub-grid boundary cell lists.
        BCels = []
        NBdys = np.zeros((3), dtype=np.int)
        for k in range(3):
            Bndyfile = DatSub+Subs[k]+'61250Bdys.dat'
            headrs, bcel = readcell( [Bndyfile] ) 
            NBdys[k] = int( headrs[0].split()[0] )
            BCels.append(bcel)
            print (Bndyfile,'cel number = %d' % NBdys[k] )

##  Find boundary cell ids for each sub-grid and merge them together.
        ncbdS = smcids(BCels[0], cel[0:nsh,         :])
        ncbdP = smcids(BCels[1], cel[  nsh:nsh+npc, :]) + nsh
        ncbdA = smcids(BCels[2], cel[  nsh+npc:ng,  :]) + nsh + npc
        ncbdy = np.hstack( (ncbdS, ncbdP, ncbdA) )

##  Find boundary cell central xlon ylat.
        xlon, ylat = smcell(ncbdy, cel, zdlnlt)

##  Generate cell vertices for the sub-grid
        nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel,
               zdlnlt, rdpols, rngsxy, excids=ncbdy, NArB=NArB )

##  Save cell vertices without boundary cells for field plots.
        ngabjm = [ng, na, nb, jmxglb]
        config=np.array([rdpols, sztpxy, rngsxy, clrbxy, ngabjm])
        pzfile=DatSub+ModlName+'Vrts.npz'
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config, 
                          svrt=svrts, scel=scels)
        print(" Verts data saved in "+pzfile )

##  Draw the global plot for SubG61250 in global view. 
        psfile=Wrkdir+ModlName+'grd.ps' 
        smcglobl( cel, nvrts,ncels,svrts,scels,colrs, config,
            mdlname=ModlName, psfile=psfile, Arctic=Arctic,
            nmark=nsmrk[0], smark=nsmrk[1] )

##  All done and print timing line.
    print( " Program finished at %s " % datetime.now().strftime('%F %H:%M:%S') )

## End of Sub61250Grid3 main program ##

if __name__ == '__main__':
    main()

