"""
##  The smcontour function plots a contour plot of a given smc grid.
##                                JGLi14Sep2022
##  Move figure and ax setup to parent program.  JGLi18May2022
##  Adapted for regular lat-lon box plot.     JGLi14Sep2022
##
"""

def smcontour(ax, lvls, cel, swhs, colrs, config, 
             grid='SMChfd', datx=' ', fontsz=9):

##  Import relevant modules and functions

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from readtext import readtext
    from colrboxy import colrboxy
    from scale_linear import scale_linear

##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

##  Total cell number by cell array
    nc = swhs.shape[0]

##  Local plot configuration parameters
    sztpxy=config[0]
    rngsxy=config[1]
    clrbxy=config[2]
    Arctic=False

##  Polar/Arctic part parameters, if defined.
    if( len(config) >3 ): 
        Arctic=True
        nabgpl=config[3]
        na= int(nabgpl[0])
        nb= int(nabgpl[1])
        nbg=int(nabgpl[2])
        npl=int(nabgpl[3])
        ng= nc - na
        print(' nabgpl defined:', nabgpl )

##  Size-1 cell dlon dlat and reference point zlon zlat.
    if( len(config) >4 ): 
        zdlnlt=config[4]
        zlon= zdlnlt[0]; zlat= zdlnlt[1]
        dlon= zdlnlt[2]; dlat= zdlnlt[3]
        print(' zdlnlt defined:', zdlnlt )

##  Define linear scale to color index conversion parameters and marks.
    factor, marks, ncstr, nclrm = scale_linear(lvls, nclrm=colrs.N)

    nswh0= ncstr
    hgts0= float(lvls[0])
    swhmn= -ncstr/factor + hgts0
    swhmx= float(lvls[-1])

##  Work out max and min values, excluding missing data (-999.0)
    cmax = swhs.max()
    cmin = swhs[ swhs > -999.0 ].min()
    print ( ' swhs range %f, %f' % ( cmin,  cmax) )
    print ( ' Draw range %f, %f' % (swhmn, swhmx) )
    cmxs = ', Vrmax=%6.3f' % cmax
    cmns = ', Vrmin=%7.4f' % cmin

##  Reset missing values (-999.0) or any value < swhmn to be swhmn 
    swhs[ swhs < swhmn ] = swhmn

##  Trim large values into plot range if any
    swhs[ swhs > swhmx ] = swhmx

##  Declare a 2-D array to hold mapped swhs values.
#   nx=int( (rngsxy[1]-rngsxy[0])/dlon ) + 1
    xlond=min([360.0, rngsxy[1]])
    nx=int( (xlond    -rngsxy[0])/dlon ) + 1
    ny=int( (rngsxy[3]-rngsxy[2])/dlat ) + 1
    nx0= int( (rngsxy[0]-zlon)/dlon )
    ny0= int( (rngsxy[2]-zlat)/dlat )
    xlon=(np.arange(nx)+nx0)*dlon
    ylat=(np.arange(ny)+ny0)*dlat

    print( " Plot box nx0 nx ny0 ny :", nx0, nx, ny0, ny )
    print( " Plot box xlon range :", xlon[0], xlon[-1] )
    print( " Plot box ylat range :", ylat[0], ylat[-1] )
    vrxy=np.zeros( (ny, nx) )
    print( " Vrxy 2-D shape = :", vrxy.shape  )

##  Convert swhs into x-y rectangular array.
##  Loop over all cells except Polar Cell, if any.
    for i in range(nc):

        if( (i < ng-nbg) or (i >= ng+nb) ):
            nxc0= cel[i,0] - nx0
            nyc0= cel[i,1] - ny0
            nxc = cel[i,2] + nxc0
            nyc = cel[i,3] + nyc0

##  Check projected cell position within range. 
            if( (0 < nxc <= nx ) and 
                (0 < nyc <= ny ) ):

                vrxy[nyc0:nyc, nxc0:nxc] = swhs[i]

##  Apply period boundary condition if full parallel circle.
    if( int(360.0/dlon) == nx-1 ):
        vrxy[:,nx-1] = vrxy[:,0] 
        print(" Period condition applied on vrxy at column", nx)
 
##  Use selected cells to draw the plot.
#   fig=plt.figure(figsize=sztpxy[0:2])
#   ax=fig.add_subplot(1,1,1)
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    ax.set(**xprop)
    ax.set(**yprop)
##  Work out ticks at 30 degree interval.
#   x0=int(rngsxy[0]/30.0); xn=int( (rngsxy[1]-rngsxy[0])/30.0 )
#   ax.set_xticks( (np.arange(min([13,xn+1]))+x0)*30 )
#   y0=int(rngsxy[2]/30.0); yn=int( (rngsxy[3]-rngsxy[2])/30.0 )
#   ax.set_yticks( (np.arange(min([7,yn+1]))+y0)*30 )
#   ax.set_aspect('equal')
    ax.set_autoscale_on(False)
#   ax.set_axis_off()

##  Draw the contour lines by specified levels. 
    ax.contour(xlon, ylat, vrxy, lvls, colors='k', linewidths=0.5)

##  Put vorticity information on top of the regular plot
#   tpx=(rngsxy[1] - rngsxy[0])*0.5
#   tpy=rngsxy[2]+1.0

#   ax.text(tpx,tpy, "Vorticity * 1.0E5 in Gs T30N160A8K3E5 at "+datx+cmns+cmxs,
#            horizontalalignment='center', fontsize=fontsz+4, color='k' )

##  Save plot as ps file
#   print (" Save the smc grid local plot ... " )
#   plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', \
#               orientation=paprorn,papertype=paprtyp)

##  End of smcontour plot function. ##

