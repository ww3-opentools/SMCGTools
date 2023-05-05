"""
##  smcelvrts function generates the polygon vertices used to draw cells
##  on grid plots or output field plots when field colour is filled. 
##  First Created on 12 Jan 2021  by Jian-Guo Li
##  Last Modified on 28 Oct 2021  by Jian-Guo Li
#
# name:     smcelvrts
#
# purpose:  genearate cell vertices using sterographic map.
#
# usage:    nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel, 
#                  zdlnlt, rdpols, rngsxy, excids=[], NArB=[] )
#
# input:    cel --- cell array to be projected. 
#           zdlnlt = [zrlon, zrlat, dxlon, dylat] --- i=j=0 lon-lat and size-1 increments.
#           rdpols = [radius, pangle, plon, plat] --- projection radius, angle, and pole.
#           rngsxy = [-10.0, 10.0,-13.0, 10.0] --- plot ranges (assume radius = 10.0 ).
#           excids = [] --- excluded cell ids to be appended at the end for marks. 
#           NArB = [] --- defaul empty list. Provided unempty for the Arctic part.
# output:   nvrts, ncels --- cell vertices and cell ids on the northern hemisphere.
#           svrts, scles --- cell vertices and cell ids on the southern hemisphere. 
#
"""

def smcelvrts( cel, zdlnlt, rdpols, rngsxy, excids=[], NArB=[] ):

    import numpy as np
    from steromap import steromap

##  Process input parameters.
    zrlon=zdlnlt[0]; zrlat=zdlnlt[1]
    dxlon=zdlnlt[2]; dylat=zdlnlt[3]
    nc = cel.shape[0]
    ng = nc
    nba = nbg = 0
    Arctic = False

    if( len(NArB) > 2 ):
        Arctic = True
        na  = int( NArB[0] )
        nba = int( NArB[1] )
        nbg = int( NArB[2] )
        ng  = nc - na

    nexc = len(excids)
    if( nexc > 0 ): print( " Excluded cell numbers =", nexc )

##  Rotated pole lon-lat and projection angle
    radius = rdpols[0]
    pangle = rdpols[1]
    plon =   rdpols[2]
    plat =   rdpols[3]

    if( pangle <= 0.0 ): 
        print( ' Pangl has to be 0 < Pangl <= 90.0!' )
        print( ' Input Pangl is equal to', pangle )
        return 

    if( radius != 10.0 ):
        print( ' Projection radius is not the default 10.0 !' )
        print( ' Input radius is equal to', radius )

##  Initial verts and ncels variable for polycollections.
    nvrts = []
    ncels = []
    svrts = []
    scels = []
    nmark  = 0        
    smark  = 0
    if( nexc > 0 ):
        nbvrt = []
        nbcel = []
        sbvrt = []
        sbcel = []

##  Loop over all cells except Polar Cell, if any.
    for i in range(nc):

        if( (i < ng-nbg) or (i >= ng+nba) ):

##  Polar cell is projected as a square box
            if( Arctic and i == nc-1 ):
                slon=np.arange(4)*90.0
                slat=np.zeros(4) + cel[i,1]*dylat + zrlat

##  Other cells are rectangular boxes by their sizes.
            else:
                xc=np.array( [cel[i,0],cel[i,0]+cel[i,2],cel[i,0]+cel[i,2],cel[i,0]] )
                yc=np.array( [cel[i,1],cel[i,1],cel[i,3]+cel[i,1],cel[i,3]+cel[i,1]] )
                slon=xc*dxlon + zrlon
                slat=yc*dylat + zrlat

##  Convert slat slon to elat elon with given new pole
            elat,elon,sxc,syc = steromap(slat,slon,plat,plon,Pangl=pangle,radius=radius,Onecl=True)

##  Check projected cell position within range and separate two hemispheres.
            if( (rngsxy[0] < sxc[0] < rngsxy[1]) and 
                (rngsxy[2] < syc[0] < rngsxy[3]) ):

##  Check whether this cell is excluded from excids
                kexc = 0
                if( nexc > 0 ):
                    for k in range(nexc):
                        if( i == excids[k] ):
                            kexc += 1
                            if( elat[0] >= 0.0 ): 
                                nbvrt.append(list(zip(sxc,syc)))
                                nbcel.append(i)
                                nmark += 1
                            else:
                                sbvrt.append(list(zip(sxc,syc)))
                                sbcel.append(i)
                                smark += 1
                            break
                
                if( kexc < 1 ):
                    if ( elat[0] >= 0.0 ): 
                        nvrts.append(list(zip(sxc,syc)))
                        ncels.append(i)
                    else:
                        svrts.append(list(zip(sxc,syc)))
                        scels.append(i)


##  Append excluded cell verts if any. 
    if( nmark > 0 ):
        nvrts = np.concatenate( [nvrts, nbvrt], axis=0 )
        ncels = np.concatenate( [ncels, nbcel], axis=0 )
    if( smark > 0 ):
        svrts = np.concatenate( [svrts, sbvrt], axis=0 )
        scels = np.concatenate( [scels, sbcel], axis=0 )
    nsmrk = [nmark, smark]

## All done. Return projected cell vertices.

    return ( nvrts, ncels, svrts, scels, nsmrk )

##  print('... Finishing smcelvrts.py ...')


##
def main():

    import numpy as np
 
    from readcell import readcell   
    from rgbcolor import rgbcolor
    from smcglobl import smcglobl

##  Read global and Arctic part cells. 
    DatGMC='../DatGMC/'
    MyCode='./'
    Wrkdir='../tmpfls'

    Cel_file = DatGMC+'G50SMCels.dat'
    Arc_file = DatGMC+'G50SMCBAr.dat'

    hdrs, cel = readcell( [Cel_file, Arc_file] ) 
    ng = int( hdrs[0].split()[0] )
    NArB = hdrs[1].split()
    na = int( NArB[0] )
    nb = int( NArB[1] )
    nc = ng + na

##  Maximum j row number in Global part
    jmxglb = np.max( cel[0:nc-na,1] )
    print (' Maximum j row =', jmxglb )

##  Extra array in config for Arctic part
    Arctic = True
    ngabjm = [ng, na, nb, jmxglb]

##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

##  Size-1 cell increments
    dxlon=0.703125
    dylat=0.468750
    zrlon=0.0
    zrlat=0.0
    zdlnlt = [zrlon, zrlat, dxlon, dylat]

##  Mark a few cells by excids for demonstration. 
    nc2 = int(nc/2 + nc/4)
    excids = [nc2-2, nc2-1, nc2, nc2+1, nc2+2]
    print(cel[excids,:])

##  Maximum mapping radius.
    radius=10.0

##  New Pole as projection direction (central point)
    pangle = 90.0
    plon= 0.0 
    plat= 23.5 
    rdpols=[radius, pangle, plon, plat]

##  Use own color map and defined depth colors 
    colrfile = MyCode+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Other plotting parameters.
    clrbxy=[ -9.6,-12.6, 19.0,  1.0]
    sztpxy=[ 16.0, 10.0,-10.2,-10.3]
    rngsxy=[-10.0, 10.0,-13.0, 10.0]
    config=np.array([rdpols, sztpxy, rngsxy, clrbxy, ngabjm])

##  Use fucntion smcelvrts to calculate cell vertices in plot.
    nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel, 
           zdlnlt, rdpols, rngsxy, excids=excids, NArB=NArB)

    print(" nsmrk =", nsmrk[0], nsmrk[1])

##  Save plotting data, which may be used by other field drawing. 
#   pzfile=Wrkdir+'SMC50Vrts.npz'
#   np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config, 
#                     svrt=svrts, scel=scels)

##  Draw the SMC50km global grid with marks.
    psfile=Wrkdir+'SMC50gridmk.ps'
    smcglobl( cel, nvrts,ncels,svrts,scels,colrs, config,
              mdlname='SMC50km', psfile=psfile, 
              nmark=nsmrk[0], smark=nsmrk[1] )
 
## End of main function.


if __name__ == '__main__':
    main()


