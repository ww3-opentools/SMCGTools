"""
#;;  steromap procedure converts standard latitude longitude to rotated 
#;;  or equatorial latitude longitude with a given new pole position and
#;;  map the rotated lat/lon with a sterographic projection.
#;;  Created on  14 May 2014  by Jian-Guo Li
#;;  Modified on 11 Nov 2014  by Jian-Guo Li
#;;  Converted into a Python function on 5 Dec 2018 by Jian-Guo Li
#
# name:     steromap
#
# purpose:  Converts standard lat/lon to rotated lat/lon and sterographic map
#
# usage:    ELat, ELon, SXxc, SYyc = steromap( SLat, SLon, Polat, Polon, Pangl=Pangl, Onecl=True )
#           Assuming default radius = 10.0, otherwise, add the radius=your_value key word.
#
# input:    SLat, SLon --- Standard lat lon as ndarray in deg
#           Polat, Polon --- New North Pole position in standard lat/lon in deg, scalors.
#           Pangl --- scalor angle in deg from rotated Pole so its projected radius is 10 unit.
#           Onecl --- if True, keep all points in one hemisphere as one cell. 
# output:   ELat, ELon --- Corresponding lat lon in rotated grid in deg, as ndarray.
#           SXxc, SYyc --- Corresponding projected x/y coordinates in range [-10, 10], as ndarray.
#
"""

import numpy as np
from ll2eqdeg import ll2eqdeg

def steromap( SLat, SLon, Polat, Polon, Pangl=90.0, radius=10.0, Onecl=False ):

##  Check Input SLat, Slon elements, should be equal
    nlat=len(SLat)
    nlon=len(SLon)

    if nlat != nlon: 
        print( ' SLat and SLon elements should have equal elements!' )
        print( ' SLat and SLon elements are', nlat, nlat )
        return 

    if Pangl <= 0.0: 
        print( ' Pangl has to be 0 < Pangl <= 90.0!' )
        print( ' Input Pangl is equal to', Pangl )
        return 

##  No need to calculate if North Pole unchanged
    if (Polat == 90.0 and Polon == 0.0): 
        ELat = SLat
        ELon = SLon
    else:
##  Convert slat slon to elat elon with given new pole
        ELat, ELon = ll2eqdeg(SLat, SLon, Polat, Polon)

##  Radian-degree conversion constants.
    d2rad=np.pi/180.0
    r2deg=180.0/np.pi

##  Number of radius from projection point.
    pnrds=4.0

##  Adjusted projected radius so that projected edge radius 
##  to be equal to the original radius. 
    PRadus=radius*(pnrds - 1.0 + np.cos(Pangl*d2rad))/np.sin(Pangl*d2rad)

##  Generate projecting coordiantes. ndarray operations.
    pradmp=PRadus*np.cos(ELat*d2rad)/(pnrds - 1.0 + np.sin(np.absolute(ELat)*d2rad))
    SYyc = pradmp*np.cos(ELon*d2rad)
    SXxc =-pradmp*np.sin(ELon*d2rad)

##  Check if any point is on the southern hemisphere
    indx = ELat < 0.0 
    if( indx.any() ):     ## Any element of this boolean array is true
        if( Onecl ):
##  If it is for one single cell projection, keep all in one hemisphere
            SXxc = -SXxc 
        else:
##  Reverse southern hemisphere x-coordinate for individual points.
            SXxc=np.where(indx, -SXxc, SXxc)

    return ( ELat, ELon, SXxc, SYyc )

##  print('... Finishing steromap.np ...')


