"""
## Program to draw the AMM 1.5-3 km SMC grid and mark boundary and dry cells.
##
## First created:    JGLi06Jul2023
## Last modified:    JGLi18Jun2025
##
"""

def main():

## Import relevant modules and functions
    import numpy as np
    import matplotlib.pyplot as plt

    from datetime import datetime
    from readcell import readcell   
    from rgbcolor import rgbcolor
    from eq2lldeg import eq2lldeg
    from smcgrids import smcgrids
    from addtexts import addtexts
    from smcelvrts import smcelvrts
    from smcellmap import smcell, smcmap, smcids

    print( " Program started at %s " % datetime.now().strftime('%F %H:%M:%S') )

## Read regional grid cell array.
    MyCode='/home/users/jianguo.li/Python/MyCodes/'
    Wrkdir='/data/users/jianguo.li/Bathy/AMM15Wave/'
    DatSMC=Wrkdir+'DatSMC/'
    GridNm='AMM153km'
    Cel_file=DatSMC+GridNm+'Cels.dat'
    Bndyfile=DatSMC+GridNm+'Bdys.dat'
    Dry_file=DatSMC+GridNm+'Drys.dat'

    headrs, cel = readcell( [Cel_file] ) 
    numbrs=headrs[0].split()
    ng = int( numbrs[0] )
    n1 = int( numbrs[1] )
    n2 = int( numbrs[2] )
    print( numbrs )

## No Arctic part so set na nb to be zero.
    Arctic = False
    na = nb = 0
    nc = ng + na
    print( GridNm+' total cell number = %d' % nc )
    print( GridNm+' N1 and N2 numbers = %d, %d' % (n1,n2) )

## Find out cell covered range.
    ic0=cel[:,0]; ic2=cel[:,2]
    istr = np.min(ic0)
    iend = np.max(ic0 + ic2) 
    jc1=cel[:,1]; jc3=cel[:,3]
    jstr = np.min(jc1)
    jend = np.max(jc1 + jc3)
    print (" Cell i range:", istr, iend)
    print (" Cell j range:", jstr, jend)

## Read grid boudnary cell array.
    headrs, bcel = readcell( [Bndyfile] )
    nbdy = int( headrs[0].split()[0] )
    print( GridNm+' bdy cell number = %d' % nbdy )

## Read dry (depth <= 0) cell array.
    headrs, dcel = readcell( [Dry_file] )
    ndry = int( headrs[0].split()[0] )
    print( GridNm+' dry cell number = %d' % ndry )

## Append dry cells to boundary cell list.
    bdcel = np.vstack( (bcel, dcel) )
    print(" Merged boundary and dry cell array shape ", bdcel.shape)

## Size-1 cell increments
    dxlon=0.0135; dylat=0.0135

## Reference or i=j=0 point lon lat.
    zrlon=-10.89625
    zrlat= -7.30095

## Rotated pole lat-lon for AMM15 grid.
    Polon = 177.5
    Polat = 37.5

## Grid related parameters
    zdlnlt = [zrlon, zrlat, dxlon, dylat]

## Initial txtary list as [txt, colr, fontsize].
    fntsz = 10.0
    fntsa = 1.2*fntsz
    fntsb = 1.5*fntsz

    txtary = [ [GridNm+' Grid',   'k', fntsb],
               ['NCel='+str(nc),  'b', fntsa],  
               ['NBdy='+str(nbdy),'r', fntsa],  
               ['NDry='+str(ndry),'r', fntsa] ] 

## Use own color map and defined depth colors 
    colrfile = MyCode+'rgbspectrum.dat'
    colrs = rgbcolor( colrfile )

## Maximum mapping radius.
    radius=10.0

    print (" Draw SMC grid for "+GridNm)

##  European regional plot
    pangle=9.5 
    plon= -1.0
    plat=  2.0
    clrbxy=[ 8.5, -9.5,  0.6,  9.0]
    sztpxy=[ 14.5, 13.0, 5.0, -6.0]
    rngsxy=[-11.0, 11.0,-10.3,10.0]
    papror='portrait'

    print( " Start loop over cells at %s " % datetime.now().strftime('%F %H:%M:%S') )

## Projection pole lon-lat and angle
    rdpols=[radius, pangle, plon, plat]

## Find boundary cell ids for the sub-grid.
    nbdry = smcids(bdcel, cel)
    ncbdy = nbdry[0:nbdy]
    print(GridNm+" len(nbdry), len(ncbdy) and nbdy =", \
                   len(nbdry), len(ncbdy),    nbdy)

## Find boundary cell central SLon SLat in bathy grid. 
    xlon, ylat = smcell(ncbdy, cel, zdlnlt)

##  Convert rotated lat-lon into standard lat-lon values.
    SLat, SLon = eq2lldeg( ylat, xlon, Polat, Polon )

## Create cell vertices
    nvrts, ncels, svrts, scels, nsmrk = smcelvrts( cel, zdlnlt, rdpols, 
                  rngsxy, excids=nbdry)

    print(" nsmrk, nbdy, ndry =", nsmrk, nbdy, ndry )

## Save boundary cell sequential number list for the sub-grid.
    fmt = '%8d '
    Bndylist = Wrkdir+GridNm+'Blst.dat'
## Boundary cell sequential numbers have to increase by 1 for WW3. JGLi06Oct2020
    np.savetxt(Bndylist, np.array(ncbdy)+1, fmt=fmt, header='',  comments='')

## Save the boundary cell SLon, SLat to generate boundary condition in WW3.
    hdr = f'{nbdy:8d} \n'
    fms = '%s '
    Bndylnlt = Wrkdir+GridNm+'Blnlt.dat'
    with open(Bndylnlt, 'w') as flhd:
        flhd.writelines(hdr)
        for j in range(nbdy):
            if( ncbdy[j] >= 0 ):
                flhd.write(f'{SLon[j]:9.3f} {SLat[j]:8.3f}   0.0  0.0  1 \n' )
    print(" Boundary cells saved in "+Bndylnlt )

## Obstruction ratio is all 0 for this grid.
    kobstr = np.zeros((nc), dtype=int)
    hdrline = f"{nc:8d} {1:5d}"
    np.savetxt(Wrkdir+GridNm+'Obst.dat', kobstr, fmt='%4d', header=hdrline, comments='')

## Set plot size and limits and message out anchor point.
    config=np.array([rdpols, sztpxy, rngsxy, clrbxy])
    pzfile=DatSMC+GridNm+'Vrts.npz'

## Store selected north and south verts and cell numbers for swh plots.
## Use the np.savez to save 3/5 variables in one file.  JGLi22Feb2019 
    np.savez( pzfile, nvrt=nvrts, ncel=ncels, cnfg=config) 

## These variables could be loaded back by
#   vrtcls = np.load(DatSMC+GridNm+'Vrts.npz')
#   nvrts = vrtcls['nvrt'] ; ncels = vrtcls['ncel']; config=vrtcls['cnfg']; 
#   svrts = vrtcls['svrt'] ; scels = vrtcls['scel']

## Draw the regional grid plot.
    epsfile=Wrkdir+GridNm+'grd.eps' 
    fig, ax = plt.subplots(figsize=sztpxy[0:2])
    smcgrids(ax, cel, nvrts,ncels,colrs,config, nmark=nsmrk[0]) 

    xydxdy=[sztpxy[2],sztpxy[3], 0.0, -0.5]
    addtexts(ax, xydxdy, txtary)

    plt.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)

    print(" ...... Saving file ", epsfile)
    plt.savefig(epsfile, dpi=None,facecolor='w',edgecolor='w', \
                    orientation=papror)
    plt.close()

    print( " Program finished at %s " % datetime.now().strftime('%F %H:%M:%S') )

## End of main() function.

if __name__ == '__main__':
    main()

## End of program AMM153kmGrids.py.

