"""
##  Adapted to generate G6-25kmSMCels from ASCII Bathymetry files. 
##                                   JGLi18Sep2014
##  Adapted for global 6-12-25-50 km SMC grid for CMA model. JGLi10Aug2018 
##  Converted from IDL into Python code.    JGLi14Feb2019 
##  Rectify size-1 cell bug and re-generate 6-50 cels.  JGLi04May2021 
##  Modified for SMCGTools to generate sample 6-50 cells.  JGLi06Oct2021 
##  Updated to read netCDF bathy and generate 2.5-5-10-20-40 km cells.  JGLi10Nov2022 
##  Separate main function from smcellgen function.   JGLi18Nov2022  
##  Updated to generate a global 2.5-5-10-20 km cells.  JGLi14Nov2022 
##  Re-tested with updated smcellgen.py function.     JGLi02Mar2023 
##  Adapted to generate SMC371250 grid with Bathy044_033.nc   JGLi18Mar2023 
##  Adapted to generate SMC61250 grid with Bathy088_059.nc   JGLi28Apr2023 
##
"""

def main():

##  Import relevant modules and functions
    import numpy   as np
    import netCDF4 as nc 

    from smcellgen import smcellgen

    Wrkdir='../tmpfls/'
    bathyf='../Bathys/Bathy088_059deg.nc' 
    Global= True
    Arctic= True
    GridNm='SMC61250'

##  Open and read bathymetry data.
    datas = nc.Dataset(bathyf)

    print(datas)            ## Print information about the dataset.

    print(datas.__dict__)

    print(datas.__dict__['description'])

    for var in datas.dimensions.values(): print(var)

    nlat = datas.dimensions['lat'].size
    nlon = datas.dimensions['lon'].size
    dlat = 180.0 / float(nlat)
    dlon = 360.0 / float(nlon)

    print(" nlat, nlon, dlat, dlon =", nlat, nlon, dlat, dlon )

    xlon = datas.variables['lon'][:]
    ylat = datas.variables['lat'][:]
    Bathy= datas.variables['elevation'][:,:]
    depthmax = np.amin(Bathy)
    hightmax = np.amax(Bathy)
    print(' Bathy range=', depthmax, hightmax)
    print(' Bathy shape=', Bathy.shape )
    print(' x_lon range=', xlon[0], xlon[nlon-1])
    print(' y_lat range=', ylat[0], ylat[nlat-1])

    datas.close()

##  Pack bathy parameters into one list.
    ndzlonlat=[ nlon, nlat, dlon, dlat, xlon[0], ylat[0] ]
##  Check top row bathy values
    print('Bathy[',nlat-1,[f'{i:d}' for i in range(0,nlon,1024)],']')
    print(Bathy[nlat-1,0:nlon:1024])

##  Decide resolution levels and SMC grid i=j=0 lon-lat.
    NLvl = 4
    x0lon = 0.0
    y0lat = 0.0
    depmin= 0.0
    mlvlxy0 = [ NLvl, x0lon, y0lat ]

    smcellgen(Bathy, ndzlonlat, mlvlxy0, FileNm=Wrkdir+GridNm, 
              Global=Global, Arctic=Arctic, depmin=depmin)
    
##  End of main program.

if __name__ == '__main__':
    main()

