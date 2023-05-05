"""
Two functions are defined here: one to calculate the cell centre
longitude and latitude from a given list of cell IDs; and another 
to map a list of longitude and latitude points with SMC cell IDs. 
Use input x0, y0, dlon, dlat and cell array for cell locating.
Return the cell id number if any point is within the cell are.  
Unmatched points are marked by 0 cell numbers. 

The main() function demonstration how the two functions could be 
called as stand alone functions.

              JGLi26Nov2020 

"""

##  idlst is an integer array of cell ids, whose centre lon-lats will be returned
##  as xlon, ylat of float arrays of the same length as idlst. 
##  The cel should be the full SMC grid cell array and zdlonlt should be a 
##  4-elements float array, containing zrlon, zrlat, dxlon, and dylat.

def smcell(idlst, cel, zdlnlt):
    """  
    Calculate cell centre longitude and latitude for a given list of cell IDs. 
    Unmatched cell ids get x0, y0 (ids<0) or dlon, dlat (ids>=ncel).
    """

    import numpy  as np

##  Check idlst and initialise xlon ylat arrays.
    if( len(idlst) > 0 ):
        mt = len(idlst)
        xlon = np.zeros( (mt), dtype=np.float )
        ylat = np.zeros( (mt), dtype=np.float )

    else:
        print (" Empty list, exiting ... ")
        return  0.0, 0.0

##  Extract x0, y0, dlon, dlat
    zlon = zdlnlt[0]; zlat = zdlnlt[1]
    dlon = zdlnlt[2]; dlat = zdlnlt[3]
    print (" zlon, zlat, dlon, dlat =", zlon, zlat, dlon, dlat)
    
##  Work out cell range by given cel array
    ncel = cel.shape[0]

##  Loop over all listed ids and work out cell centre x, y
    im = 0
    for j in range(mt):
        ids = int(idlst[j])
        if( ids < 0 or ids >= ncel ):
            xlon[j] = 0.0; ylat[j] = -90.0
            im += 1
        else:
            xlon[j] = zlon + (float(cel[ids,0]) + 0.5*float(cel[ids,2]))*dlon
            ylat[j] = zlat + (float(cel[ids,1]) + 0.5*float(cel[ids,3]))*dlat

    if( im > 0 ): print(" Missed number of points:", im )

    return  xlon, ylat  


##
def smcmap(xlon, ylat, cel, zdlnlt):
    """  
    Match a given list of lon-lat points to SMC grid cell IDs. 
    Unmatched points are marked with 0 cell IDs.
    """

    import numpy  as np

##  Check input lon-lat list. 
    nx = len(xlon); ny = len(ylat)
    if( nx == ny > 0 ):
        idlst = np.zeros( (nx), dtype=np.int ) -1

    else:
        print (" Empty lon-lats, exiting ... ")
        return 0 

##  Extract zlon, zlat, dlon, dlat
    zlon = zdlnlt[0]; zlat = zdlnlt[1]
    dlon = zdlnlt[2]; dlat = zdlnlt[3]
    print (" zlon, zlat, dlon, dlat =", zlon, zlat, dlon, dlat)
    
##  Convert xlon to ensure xlon >= x0
    xlow = np.where( xlon < x0, xlon+360.0, xlon )

##  Work out cell range by given cel array
    ncel = cel.shape[0]

##  Convert x, y into cell indexes
    ixlw = np.floor( (xlow-zlon)/dlon ).astype(np.int)
    jylt = np.floor( (ylat-zlat)/dlat ).astype(np.int)

##  Show coversion results.
#   print("\n Converted i,j and xlow, ylat:")
#   for i in range(nx): 
#       print(f'{ixlw[i]:8d}{jylt[i]:8d} {xlow[i]:9.3f}{ylat[i]:9.3f}')

##  Loop over all cell array until all points are matched. 
    i = 0
    mfnd = 0
#   print("\n Matched cell array:")
    while (i < ncel and mfnd < nx):
        i0=cel[i,0]; i2=cel[i,2]
        j1=cel[i,1]; j3=cel[i,3]

##  Loop over all listed points to match with the cell  
        for jp in range(nx): 
            if( idlst[jp] < 0 ):
##  Check this point with cell i
                if( i0 <= ixlw[jp] < i0+i2 and 
                    j1 <= jylt[jp] < j1+j3 ):
                    mfnd += 1
                    idlst[jp] = i
#                   print(f'{i:8d}{i0:8d}{j1:8d}{i2:6d}{j3:6d}')
                    break

        i += 1

##  End of while loop.


    if( mfnd < nx ): print(" Missed number of points:", nx-mfnd)

    return idlst 

##
def smcids(idcel, cel):
    """  
    Find cell ids for all cells in idcel from the full cell array cel. 
    Unmatched cells are asigned ids to be 0. 
    """

    import numpy  as np

##  Check idcel to be unempty.
    nids = idcel.shape[0]
    if( nids < 1 ):
        print (" Empty id cell list, exiting ... ")
        return  [-1] 
    else:
        idlst = np.zeros( (nids), dtype=np.int ) -1

##  Work out cell range by given cel array
    ncel = cel.shape[0]

##  Loop over all idcel and work out their cell ids, if found.
    im = 0
    for k in range(nids):
        ix = idcel[k,0]; jy = idcel[k,1]
        for n in range(ncel):
            if( ix == cel[n,0] and jy == cel[n,1] ):
                if( idcel[k,2] == cel[n,2] and idcel[k,3] == cel[n,3] ):
                    idlst[k] = n
                    im += 1
                    break

    print(" Matched number of cells:", im )

    return  idlst 


def main():

    import numpy  as np
    from readcell import readcell   

    DatGMC='../DatGMC/'
    Wrkdir='../tmpfls/'
    SMC_file = DatGMC+'SMC61250Cels.dat'
    Arc_file = DatGMC+'SMC61250BArc.dat'
  
    headrs, cel = readcell( [SMC_file, Arc_file] ) 
    ng = int( headrs[0].split()[0] )
    na = int( headrs[1].split()[0] )
    nc = na + ng

    x0 = 0.0; y0 = 0.0
    dxlon=0.0439453125*2.0
    dylat=0.029296875*2.0

    zdlnlt = [x0, y0, dxlon, dylat]

    idlst = [-1, 0, 10, ng-1, nc-1, nc]

    print("\n Initial ids and cell array if in range:")
    for i in range(len(idlst)):
        j = idlst[i]
        if( 0 <= j < nc ):
            print(f'{j:8d}{cel[j,0]:8d}{cel[j,1]:8d}{cel[j,2]:6d}{cel[j,3]:6d}')
        else:
            print(f'{j:8d}')

#   Convert idlst into xlon, ylat
    xlon, ylat = smcell(idlst, cel, zdlnlt)

#   Convert xlon, ylat back into new idlst
    idsnew = smcmap(xlon, ylat, cel, zdlnlt)

    print("\n Initial ids, xlon, ylat, final ids:")
    for i in range(len(idlst)):
        print(f'{idlst[i]:8d} {xlon[i]:9.3f}{ylat[i]:9.3f}{idsnew[i]:8d}')

#   Find ids of first 10 cells.
    idfst10 = smcids(cel[0:10,:], cel)
    print("\n First 10 cells ids are:", idfst10)

#   End of main function.


if __name__ == '__main__':
    main()

