"""
##  smcelvrglr function generates the polygon vertices used to draw cells
##  on grid plots or output field plots when field colour is filled. 
##  First Created on 12 Jan 2021  by Jian-Guo Li
##  Last Modified on 12 Sep 2022  by Jian-Guo Li
#
# name:     smcelvrglr
#
# purpose:  genearate cell vertices using lat-lon rectangular map.
#
# usage:    nvrts, ncels, nmark = smcelvrglr( cel, 
#                  zdlnlt, rngsxy, excids=[], NArB=[])
#
# input:    cel --- cell array to be projected. 
#           zdlnlt = [zrlon, zrlat, dxlon, dylat] --- i=j=0 lon-lat and size-1 increments.
#           rngsxy = [0.0, 360.0,-90.0, 90.0] --- plot ranges (assume lat-lon grid).
#           excids = [] --- excluded cell ids to be appended at the end for marks. 
#           NArB = [] --- default empty list. Provided unempty for the Arctic part.
# output:   nvrts, ncels --- cell vertices and cell ids on the northern hemisphere.
#           svrts, scles --- cell vertices and cell ids on the southern hemisphere. 
#
"""

def smcelvrglr( cel, zdlnlt, rngsxy, excids=[], NArB=[] ):

    import numpy as np

##  Process input parameters.
    zrlon=zdlnlt[0]; zrlat=zdlnlt[1]
    dxlon=zdlnlt[2]; dylat=zdlnlt[3]
    nc = cel.shape[0]
    ng = nc
    nba = nbg = npl = 0
    Arctic = False

    if( len(NArB) > 2 ):
        Arctic = True
        na  = int( NArB[0] )
        nba = int( NArB[1] )
        nbg = int( NArB[2] )
        ng  = nc - na

##  Default 1 polar cell.
        if( len(NArB) > 3 ):
            npl = int( NArB[3] )
        else:
            npl = 1

    nexc = len(excids)
    if( nexc > 0 ): print( " Marked cell numbers =", nexc )

##  Check lat-lon box range and return if not consistent.
    if( rngsxy[0] < rngsxy[1] ): 
        print( ' Longitude range:', rngsxy[0:2] )
    else:
        print( ' Inconsistent longitude range:', rngsxy[0:2] )
        return 

    if( rngsxy[2] < rngsxy[3] ):
        print( ' Latitude  range:', rngsxy[2:4] )
    else:
        print( ' Inconsistent latitude range:', rngsxy[2:4] )
        return 

##  Initial verts and ncels variable for polycollections.
    nvrts = []
    ncels = []
    nmark  = 0        
    if( nexc > 0 ):
        nbvrt = []
        nbcel = []

##  Loop over all cells except Polar Cell, if any.
    for i in range(nc):

        if( (i < ng-nbg) or (i >= ng+nba) ):

##  Polar cell is projected as a square box
            if( Arctic and i >= nc-npl ):
                sxc=np.array( [0.0, 360.0, 360.0, 0.0] )
                yp1=cel[i,1]*dylat + zrlat
                yp2=(cel[i,1]+cel[i,3])*dylat + zrlat
                syc=np.array( [yp1, yp1, yp2, yp2] ) 

##  Other cells are rectangular boxes by their sizes.
            else:
                xc=np.array( [cel[i,0],cel[i,0]+cel[i,2],cel[i,0]+cel[i,2],cel[i,0]] )
                yc=np.array( [cel[i,1],cel[i,1],cel[i,3]+cel[i,1],cel[i,3]+cel[i,1]] )
                sxc=xc*dxlon + zrlon
                syc=yc*dylat + zrlat

##  Check projected cell position within range. 
            if( (rngsxy[0] <= sxc[0] < rngsxy[1]) and 
                (rngsxy[2] <= syc[0] < rngsxy[3]) ):

##  Check whether this cell is excluded from excids
                kexc = 0
                if( nexc > 0 ):
                    for k in range(nexc):
                        if( i == excids[k] ):
                            kexc += 1
                            nbvrt.append(list(zip(sxc,syc)))
                            nbcel.append(i)
                            nmark += 1
                            break
                
                if( kexc < 1 ):
                    nvrts.append(list(zip(sxc,syc)))
                    ncels.append(i)

##  Append excluded cell verts if any. 
    if( nmark > 0 ):
        nvrts = np.concatenate( [nvrts, nbvrt], axis=0 )
        ncels = np.concatenate( [ncels, nbcel], axis=0 )

## All done. Return projected cell vertices.

    return ( nvrts, ncels, nmark )

##  print('... Finishing smcelvrglr.py ...')


##
def main():

    import numpy as np
    import matplotlib.pyplot as plt
 
    from readcell import readcell   
    from rgbcolor import rgbcolor
    from smcrglgrd import smcrglgrd

##  Read global and Arctic part cells. 
    DatGMC='../DatGMC/'
    MyCode='../PySMCs/'
    Wrkdir='./'

    Cel_file = DatGMC+'G50SMCels.dat'
    Arc_file = DatGMC+'G50SMCBAr.dat'

    hdrs, cel = readcell( [Cel_file, Arc_file] ) 
    ng = int( hdrs[0].split()[0] )
    NArB = hdrs[1].split()
    na = int( NArB[0] )
    nb = int( NArB[1] )
    nbg = int( NArB[2] )
    npl = 1
    NArB.append(npl)
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

##  Use own color map and defined depth colors 
    colrfile = MyCode+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

##  Other plotting parameters.
    clrbxy=[ 5.0, -89.0, 350.0, 6.0]
    sztpxy=[ 16.0, 11.0, 60.0, 30.0]
    rngsxy=[0.0, 360.0, -90.0, 90.0]
    config=np.array([sztpxy, rngsxy, clrbxy, ngabjm])

##  Use fucntion smcelvrts to calculate cell vertices in plot.
    nvrts, ncels, nmark = smcelvrglr( cel, 
           zdlnlt, rngsxy, excids=excids, NArB=NArB)

    print(" nmark =", nmark)

##  Save plotting data, which may be used by other field drawing. 
#   pzfile=Wrkdir+'SMC50rglVrts.npz'
#   np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 

##  Draw the SMC50km global grid with marks.
    psfile=Wrkdir+'SMC50gridrglr.ps'
    fig=plt.figure(figsize=sztpxy[0:2])
    ax=fig.add_subplot(1,1,1)

    smcrglgrd(ax, cel, nvrts, ncels, colrs, config, 
             grid='SMC50km', fontsz=11, nmark=nmark )
 
##  Save plot as ps file
    print (" Save the smc grid local plot ... " )
    plt.subplots_adjust(left=0.03,bottom=0.03,right=0.99,top=0.98)
    plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', \
                orientation='landscape',papertype='a3')

## End of main function.


if __name__ == '__main__':
    main()


