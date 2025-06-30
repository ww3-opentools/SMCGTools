"""
## Program to draw the Sub61250 sub-grids and merged one.
##
## First created:        JGLi12Jan2021
## Last modified:        JGLi02May2025
##
"""

def main():

## Import relevant modules and functions
    import numpy as np
    import matplotlib.pyplot as plt

    from datetime import datetime
    from readcell import readcell   
    from rgbcolor import rgbcolor
    from smcgrids import smcgrids
    from addtexts import addtexts
    from smcelvrts import smcelvrts
    from smcellmap import smcell, smcmap, smcids

    print( " Program started at %s " % datetime.now().strftime('%F %H:%M:%S') )

## Read global and Arctic part cells. 
    DatSub='../DatSub/'
    Wrkdir='../tmpfls/'
    Arcell='../DatGMC/SMC61250BArc.dat'

## Size-1 cell increments
    dxlon=0.087890625
    dylat=0.058593750
#   dxlon=0.0439453125
#   dylat=0.029296875
    zrlon=0.0
    zrlat=0.0
    zdlnlt = [zrlon, zrlat, dxlon, dylat]

## Use own color map and defined depth colors 
    colrfile = 'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

## Maximum mapping radius.
    radius= 10.0
    fntsz = 10.0
    fntsa = 1.20*fntsz
    fntsb = 1.50*fntsz

## Possible selection of your plot types. 
    gorloc={0:'SubG',1:'Soth',2:'Pacf',3:'Atln'}

## Prompt selection choices and ask for one input
    print (" \n ", gorloc)
    instr = input(' *** Please enter your selected number here > ')
    m = int(instr)
    pltype=gorloc.get(m, 'Invalid_selection')
    if( pltype == 'Invalid_selection' ): 
        print ("Invalid selection, program terminated.")
        exit()

    print( " Start draw "+pltype+" at %s " % datetime.now().strftime('%F %H:%M:%S') )

    if( pltype == 'Soth' ):
## Whole global projection angle from N Pole to be 90.0
        pangle= 90.0
        plon=-136.0
        plat= 21.0
        clrbxy=[ -9.0,  3.6,   7.0,  0.7]
        sztpxy=[ 12.0, 11.0,  -7.0,  8.0]
        rngsxy=[-10.0, 10.0,  -9.5,  9.0]
        papror='portrait'

    if( pltype == 'Pacf' ):
## Whole global projection angle from N Pole to be 90.0
        pangle=90.0
        plon=  10.0 
        plat= 22.6 
        clrbxy=[ -9.3, -9.0,  9.9,  0.9]
        sztpxy=[ 11.0, 10.5, -7.5, -2.5]
        rngsxy=[-10.1, 10.1, -9.6, 10.1]
        papror='portrait'
        rdpols=[radius, pangle, plon, plat]

## Atlantic regional plot for Atln61250 sub-grid
    if( pltype == 'Atln'):
        pangle= 50.0 
        plon=-53.0
        plat= 51.0
        clrbxy=[ -9.0,  1.1,   0.7,  7.0]
        sztpxy=[ 10.0, 12.0,  -7.0, -7.6]
        rngsxy=[-10.0, 10.0, -12.0, 12.0]
        papror='portrait'

## Full global plot of all sub-grids.
    if( pltype == 'SubG' ):
        pangle=90.0
        plon= 0.0
        plat= 23.5
        clrbxy=[ -9.6,-12.0, 19.0,  1.0]
        sztpxy=[ 16.0,  9.0,-10.2,-10.3]
        rngsxy=[-10.0, 10.0,-12.1, 10.0]

## Radius, projection angle, and pole lon-lat.
    rdpols=[radius, pangle, plon, plat]

## Model name for selected grid. 
    ModlName=pltype+'61250'
    epsfile=Wrkdir+ModlName+'grd.eps' 
    spbuofl='../Bathys/SPBuoys.dat'
    paprorn='portrait'

    if( m > 0 ):
## Process sub-grid and its boundary cells
        Cellfile = [DatSub+ModlName+'Cels.dat']
        Bndyfile = [DatSub+ModlName+'Bdys.dat']
## Append Arctic part for Atln sub-grid.
        if( m == 3 ):
            Cellfile.append(Arcell)
 
## Read sub-grid cell
        headrs, cel = readcell( Cellfile )
        nc = int( headrs[0].split()[0] )
        NArB = []
        if( m == 3 ):
            NArB = headrs[1].split()
            na = int( NArB[0] )
            nb = int( NArB[1] )
            nbg= int( NArB[2] )
            npl= 1
            ng = nc 
            nc = ng + na 

## Read sub-grid boudnary cell
        headrs, bcel = readcell( Bndyfile ) 
        nbdy = int( headrs[0].split()[0] )

## Find boundary cell ids for the sub-grid.
        ncbdy = smcids(bcel, cel)

## Find boundary cell central xlon ylat. 
        xlon, ylat = smcell(ncbdy, cel, zdlnlt)
        
## Generate cell vertices for the sub-grid
        nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel, 
               zdlnlt, rdpols, rngsxy, excids=ncbdy, NArB=NArB )

## Save boundary cell sequential number list for Pacf61250 grid.
        fmt = '%8d '
        Bndylist = Wrkdir+ModlName+'Blst.dat'
## Boundary cell sequential numbers have to increase by 1 for WW3. JGLi06Oct2020
        np.savetxt(Bndylist, np.array(ncbdy)+1, fmt=fmt, header='',  comments='')

## Save the boundary cell xlon, ylat to generate boundary condition in WW3.
        hdr = f'{nbdy:8d} \n'
        fms = '%s '
        Bndylnlt = Wrkdir+ModlName+'Blnlt.dat'
        with open(Bndylnlt, 'w') as flhd:
            flhd.writelines(hdr)
            for j in range(nbdy):
                if( ncbdy[j] >= 0 ): 
                    flhd.write(f'{xlon[j]:9.3f} {ylat[j]:8.3f}   0.0  0.0  1 \n' )
        print(" Boundary cells saved in "+Bndylnlt )

## For Soth and Pacf sub-grids draw the southern hemisphere
        if( pltype == 'Soth' or pltype == 'Pacf' ):
            Arctic= False
            np000 = [nc, 0, nsmrk[0], nsmrk[1]]
            config=np.array([rdpols, sztpxy, rngsxy, clrbxy, np000])
            pzfile=DatSub+ModlName+'Vrts.npz'
            np.savez( pzfile, nvrt=svrts, ncel=scels, cnfg=config) 
            print(" Verts data saved in "+pzfile )

            fig, ax = plt.subplots(figsize=sztpxy[0:2])

            smcgrids(ax, cel, svrts,scels,colrs,config, 
                Arctic=Arctic, nmark=nsmrk[1]) 

            xydxdy=[sztpxy[2], sztpxy[3]-0.6, 0.0, -0.6]
            txtary= [ [ModlName+' Grid',  'k', fntsb],
                      ['NC='+str(nc),     'b', fntsa],
                      ['NBdy='+str(nbdy), 'r', fntsa] ]
            addtexts(ax, xydxdy, txtary)
            plt.subplots_adjust(left=0,bottom=0,right=1,top=1)

## For Atln sub-grid draw the northern hemisphere
        if( pltype == 'Atln' ):
            Arctic = True 
            nabgpl = [na, nb, nbg, npl]
            config=np.array([rdpols, sztpxy, rngsxy, clrbxy, nabgpl])
            pzfile=DatSub+ModlName+'Vrts.npz'
            np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 
            print(" Verts data saved in "+pzfile )

            fig, ax = plt.subplots(figsize=sztpxy[0:2])

            smcgrids(ax, cel, nvrts,ncels,colrs,config,
                Arctic=Arctic, nmark=nsmrk[0])

            xydxdy=[sztpxy[2], sztpxy[3]-0.6, 0.0, -0.6]
            txtary= [ [ModlName+' Grid',  'k', fntsb],
                      ['NC='+str(nc),     'b', fntsa],
                      ['NA='+str(na),     'b', fntsa],
                      ['NBdy='+str(nbdy), 'r', fntsa] ]
            addtexts(ax, xydxdy, txtary)
            plt.subplots_adjust(left=0,bottom=0,right=1,top=1)

## If m = 0 or full global grid.
    if( pltype == 'SubG' ):
## Read all sub-grid cells plus Arctic part.
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
        nbg= int( NArB[2] )
        npl= 1
        nc = ng + na
        jmxglb = np.max( cel[0:nc-na,1] )
        print (' Maximum j row =', jmxglb)
        print(" nc ng na nb jmxglb=", nc, ng, na, nb, jmxglb)

## Read sub-grid boundary cell lists.
        BCels = []
        NBdys = np.zeros((3), dtype=int)
        for k in range(3):
            Bndyfile = DatSub+Subs[k]+'61250Bdys.dat'
            headrs, bcel = readcell( [Bndyfile] ) 
            NBdys[k] = int( headrs[0].split()[0] )
            BCels.append(bcel)
            print (Bndyfile,'cel number = %d' % NBdys[k] )

## Find boundary cell ids for each sub-grid and merge them together.
        ncbdS = smcids(BCels[0], cel[0:nsh,         :])
        ncbdP = smcids(BCels[1], cel[  nsh:nsh+npc, :]) + nsh
        ncbdA = smcids(BCels[2], cel[  nsh+npc:ng,  :]) + nsh + npc
        ncbdy = np.hstack( (ncbdS, ncbdP, ncbdA) )

## Find boundary cell central xlon ylat.
        xlon, ylat = smcell(ncbdy, cel, zdlnlt)

## Generate cell vertices for the sub-grid
        nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel,
               zdlnlt, rdpols, rngsxy, excids=ncbdy, NArB=NArB )

## Save cell vertices without boundary cells for field plots.
        nabgpl = [na, nb, nbg, npl]
        config=np.array([rdpols, sztpxy, rngsxy, clrbxy, nabgpl])
        pzfile=DatSub+ModlName+'Vrts.npz'
        np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config, 
                          svrt=svrts, scel=scels)
        print(" Verts data saved in "+pzfile )

## Draw the global plot for SubG61250 in global view. 
        fig=plt.figure(figsize=sztpxy[0:2])
        ax1=fig.add_subplot(1,2,1)

        smcgrids(ax1, cel, nvrts,ncels,colrs, config, Arctic=True, 
                 nmark=nsmrk[0], buoys=spbuofl, hemis=1.0)

        ax2=fig.add_subplot(1,2,2)
        smcgrids(ax2, cel, svrts,scels,colrs, config, Arctic=False, 
                 nmark=nsmrk[1], buoys=spbuofl, hemis=-1.0)

        xydxdy=[-10.2, -9.8,  0.0,  0.6]
        txtary=[ [ModlName+' Grid',  'k', fntsb],
                 ['NC='+str(nc),     'b', fntsa],
                 ['NA='+str(na),     'b', fntsa],
                 ['NB='+str(nb),     'r', fntsa] ]               
        addtexts(ax2, xydxdy, txtary)
        plt.subplots_adjust(left=0.01, bottom=0.0, right=0.99, 
                             top=1.0, wspace=0.01, hspace=0.0)

## Save grid plot as eps file.
    print (" ... saving the smc grid plot as ", epsfile )
    plt.savefig(epsfile, dpi=None,facecolor='w',edgecolor='w', \
                orientation=paprorn)
    plt.close()

##  All done and print timing line.
    print( " Program finished at %s " % datetime.now().strftime('%F %H:%M:%S') )

## End of main() function. 

if __name__ == '__main__':
    main()

## End of Sub61250Grid3.py program.

