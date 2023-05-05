"""
#; Program to define regular grid parameters for SMC61250 grid.
## Converted from IDL to Python.  JGLi26Feb2019
#; First created:   25 Feb 2015   Jian-Guo Li
#; Last modified:   08 Dec 2022   Jian-Guo Li
##
"""
 
def main():
    """  Generate regular grid parameters for SMC61250 grid. """

    import numpy as np
    from regulrinfo import regulrinfo

##  Define sub-grids and parameters. 
    Wrkdir='./'
    Model='SMC61250'

##  Read SMC61250 grid info from Grindifo.dat 
    zrdlonlat = np.genfromtxt('Gridinfo.dat', dtype=float, skip_header=1)
    print(" Input file zlon zlat dlon dlat = \n", zrdlonlat)

    NRLv= 4
    NLon= 4096
    NLat= 3072
    NLvLonLat=[ NRLv, NLon, NLat ]

    print("\n Model NRLv NLon NLat = ", Model, NRLv, NLon, NLat ) 

    NRLout = regulrinfo( Model, zrdlonlat, NLvLonLat, FLat=-79.8 ) 


##  End of main program.

if __name__ == '__main__':
    main()

