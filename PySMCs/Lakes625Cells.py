"""
##  Adapted to generate G6-25kmSMCels from ASCII Bathymetry files. 
##                                   JGLi18Sep2014
##  Adapted for global 6-12-25-50 km SMC grid for CMA model. JGLi10Aug2018 
##  Converted from IDL into Python code.    JGLi14Feb2019 
##  Rectify size-1 cell bug and re-generate 6-50 cels.  JGLi04May2021 
##  Modified for SMCGTools to generate sample 6-50 cells.  JGLi06Oct2021 
##  Updated to read netCDF bathy and generate 2.5-5-10-20-40 km cells.  JGLi10Nov2022 
##  Separate main function from smcellgen function.   JGLi18Nov2022  
##  Modified to generate European 1.25-2.5-5-10 km grid.  JGLi29Nov2022
##  Modified to generate Caspian Sea 2.5-5-10-20 km grid.  JGLi28Feb2023
##  Modified to generate Caspian Lake Superior 3-7-13-25 km grid.  JGLi22Mar2023
##  Modified to generate Caspian Lake Superior 6-12-25 km grid.  JGLi28Apr2023
##
"""

def main():

##  Import relevant modules and functions
    import numpy   as np
    import netCDF4 as nc 

    from smcellgen import smcellgen

    Wrkdir='../tmpfls/'
    bathyf='../Bathys/Bathy088_059deg.nc'
    Global= False
    Arctic= False

##  Open and read bathymetry data.
    datas = nc.Dataset(bathyf)

    print(datas)            ## Print information about the dataset.

    print(datas.__dict__)

    print(datas.__dict__['description'])

    for var in datas.dimensions.values(): print(var)

    nlat = datas.dimensions['lat'].size
    nlon = datas.dimensions['lon'].size
    ylat = datas.variables['lat'][:]
    xlon = datas.variables['lon'][:]
    dlat = (ylat[-1] - ylat[0])/float(nlat-1)
    dlon = (xlon[-1] - xlon[0])/float(nlon-1)
    print(" nlat, nlon, dlat, dlon =", nlat, nlon, dlat, dlon )

    Bathy= datas.variables['elevation'][:,:]
    depthmax = np.amin(Bathy)
    hightmax = np.amax(Bathy)
    print(' Bathy range=', depthmax, hightmax)
    print(' Bathy shape=', Bathy.shape )

    print( " Domain range x0, y0, xn, yn =",xlon[0],ylat[0],xlon[-1],ylat[-1] )

##  Pack bathy parameters into one list.
    ndzlonlat=[ nlon, nlat, dlon, dlat, xlon[0], ylat[0] ]

    datas.close()

##  Check top row bathy values
    print('Bathy[',nlat-1,[f'{i:d}' for i in range(0,nlon,1024)],']')
    print(Bathy[nlat-1,0:nlon:1024])

##  Patch last row values to be exactly the ones in the next inner row.
#   Bathy[nlat-1,:] = Bathy[nlat-2,:]

##  Check last column bathy values
#   print('Bathy[',[f'{i:d}' for i in range(0,nlat,225)],',', nlon-1,']')
#   print(Bathy[0:nlat:225,nlon-1])
#   print('Bathy[',[f'{i:d}' for i in range(0,nlat,225)],',', 0,']')
#   print(Bathy[0:nlat:225,0])

##  Decide resolution levels and SMC grid i=j=0 lon-lat.
    NLvl = 3

##  i/j start point set at zero meridian on the Equator.
    x0lon = 0.0;         y0lat=0.0

##  Lake range area (SW and NE cornner Lon and Lat) and low water level altitude.
    Caspian =[ 45.0, 34.0,  60.0, 50.0]
    Superior=[-93.0, 46.0, -83.5, 49.3]
    Michigan=[-89.0, 41.0, -79.0, 47.0]
    ErieSClr=[-84.0, 40.5, -78.0, 43.5]
    Ontario =[-80.5, 42.5, -75.5, 44.5]
    GridNames=['Casp','Supr','MHur','Erie','Onta']
    LakeRange=[Caspian,Superior,Michigan,ErieSClr,Ontario]
    Lowaterlv=[ -27.5, 183.3,   176.0,   173.5,   74.2   ]

    for grd in range(len(GridNames)):
        GridNm=GridNames[grd]+'6125'
        Ranges=LakeRange[grd]
        depmin=Lowaterlv[grd]
##  Note depmin defines water surface altitude (normally sea level 0.0)
##  and non-zero values are used to define lake surfaces differ from sea level.
##  dshalw defines shallow water altitude and is used to refine shallow ater area
##  so setting dshalw equal to depmin makes all area below depmin to be deep water.
##  Define dshalw = depmin - 60.0 will refine areas less than 60 m below depmin.
        dshalw=depmin
 
        print( GridNm+" range and water level =", Ranges, depmin  ) 

##  Merte lake range with mlvlxy0 array. 
        mlvlxy0 = [ NLvl, x0lon, y0lat ] + Ranges 

        smcellgen(Bathy, ndzlonlat, mlvlxy0, FileNm=Wrkdir+GridNm, 
            Global=Global, Arctic=Arctic, depmin=depmin, dshalw=dshalw) 
    
        print( "GridNm cells saved in ", Wrkdir+GridNm+'cels.dat')

##  End of lake grid loop.
    
##  End of main program.


if __name__ == '__main__':
    main()

