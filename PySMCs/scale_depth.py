"""
#;;  scale_depth function to set up the color scale for depth plot, using the
#;;  first 130 colors from a given color table of 256 colors.
#;;  First Created on  22 Feb 2019  by Jian-Guo Li
#;;  Last Modified on  26 Feb 2019  by Jian-Guo Li
#
# usage:    depth, factr, cstar, marks, ckeys = scale_depth( nclrm=131 )
#
# Bathymetry depth (positive) will need to be converted by
#    ndeps=np.rint( (cstar - np.log10(depth))*factr ).astype(np.int)
# to be consistent with the color key. 
#
"""

def scale_depth( ncstr=0, nclrm=131 ):

    import numpy as np

    depth=np.array([10000,1000,100,10,1], dtype=np.int)
    cstar=np.log10(depth[0]+1000.0)
    factr=(nclrm - 1)/cstar
    marks= ncstr + np.rint( (cstar - np.log10(depth))*factr ).astype(np.int)
    print( nclrm, ' colors with depth marks at ', depth[::-1])

    return ( depth, factr, cstar, marks, ncstr, nclrm )


