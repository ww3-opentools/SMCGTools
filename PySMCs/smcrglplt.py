"""
##  The smcrgulr function plots a local region of a given smc grid.
##  cel is the full cell array in shape(nc, 5).
##                                JGLi04Mar2019
##  Modified to fill marked cells to be red.  JGLi03Sep2020
##  Move figure and ax setup to parent program.  JGLi18May2022
##  Adapted for regular lat-lon box plot.     JGLi12Sep2022
##
"""

def smcrglplt(ax, hgts, swhs, verts, ncels, colrs, config, 
             grid='SMChfd',datx=' ',panel=' ',Hunit=' ', fontsz=9):

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

##  Define linear scale to color index conversion parameters and marks.
    factor, marks, ncstr, nclrm = scale_linear(hgts, nclrm=colrs.N)

    nswh0= ncstr
    hgts0= float(hgts[0])
    swhmn= -ncstr/factor + hgts0
    swhmx= float(hgts[-1])

##  Work out max and min values, excluding missing data (-999.0)
    cmax = swhs.max()
    cmin = swhs[ swhs > -999.0 ].min()
    print ( ' swhs range %f, %f' % ( cmin,  cmax) )
    print ( ' Draw range %f, %f' % (swhmn, swhmx) )
    cmxs = ', Hmax=%6.3f' % cmax
    cmns = ', Hmin=%7.4f' % cmin

##  Reset missing values (-999.0) or any value < swhmn to be swhmn 
    swhs[ swhs < swhmn ] = swhmn

##  Trim large values into plot range if any
    swhs[ swhs > swhmx ] = swhmx

##  Convert swhs with liniear scale.
    icnf = ncstr+ np.rint( factor*(swhs - hgts0) )
    nswh = np.array( icnf, dtype=int )

##  Use selected cells to draw the plot.
#   fig=plt.figure(figsize=sztpxy[0:2])
#   ax=fig.add_subplot(1,1,1)
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    ax.set(**xprop)
    ax.set(**yprop)
##  Work out ticks at 30 degree interval.
    x0=int(rngsxy[0]/30.0); xn=int( (rngsxy[1]-rngsxy[0])/30.0 )
    ax.set_xticks( (np.arange(min([13,xn+1]))+x0)*30 )
    y0=int(rngsxy[2]/30.0); yn=int( (rngsxy[3]-rngsxy[2])/30.0 )
    ax.set_yticks( (np.arange(min([7,yn+1]))+y0)*30 )
#   ax.set_aspect('equal')
    ax.set_autoscale_on(False)
#   ax.set_axis_off()

##  Create color array for this plot
    pcface = []
    pcedge = []
    for i in ncels:
##  Fill last nmark cells to be red.  JGLi03Sep2020
        if( nswh[i] != nswh0 ): 
            pcface.append(colrs(nswh[i]))
            pcedge.append(colrs(nswh[i]))
        else:
            pcface.append(colrs(255))
            pcedge.append(colrs(0))

##  Create PolyCollection from selected verts and define edge and face color.
    polynorth = PolyCollection(verts)
    polynorth.set_facecolor(pcface)
    polynorth.set_edgecolor(pcedge)
    polynorth.set_linewidth( 0.2 )

##  Draw the selected cells as colored polygons.
    ax.add_collection(polynorth)

##  Draw colorbar inside plot.
    xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks, nclrm=nclrm)
    ax.add_collection(clrply)
    dkx=clrbxy[2]; dky=clrbxy[3]
    for i in range(len(hgts)):
        m = marks[i]
        if( dkx < dky ):    
            ax.text(xkeys[0]+1.15*dkx, ykeys[m], str(hgts[i]),
                verticalalignment='center', fontsize=fontsz, color='b' )
        else:
            ax.text(xkeys[m], ykeys[0]+1.15*dky, str(hgts[i]),
                horizontalalignment='center', fontsize=fontsz, color='b' )

    if( dkx < dky ):    
        ax.text(xkeys[0]+1.2*dkx, (ykeys[marks[2]]+ykeys[marks[1]])*0.5, 'H '+Hunit,
                 rotation=-90,verticalalignment='center', fontsize=fontsz+2, color='k' )
    else:
        ax.text((xkeys[marks[1]]+xkeys[marks[2]])*0.5, ykeys[0]+1.2*dky, 'H '+Hunit,
                 rotation=0,horizontalalignment='center', fontsize=fontsz+2, color='k' )

##  Put vorticity information on top of the regular plot
#   tpx=sztpxy[2] 
#   tpy=sztpxy[3]
#   dpy=3.0 
    tpx=(rngsxy[1] - rngsxy[0])*0.5
    tpy=rngsxy[2]+1.0

    ax.text(tpx,tpy, panel+datx+cmns+cmxs+Hunit,
             horizontalalignment='center', fontsize=fontsz+4, color='k' )

#   ax.text(tpx,-tpy-dpy*1.5, panel,
#            horizontalalignment='left', fontsize=fontsz+4, color='k' )

#   ax.text(tpx, tpy+dpy*0, grid+' Grid',  
#            horizontalalignment='left', fontsize=fontsz+4, color='k' )
#   ax.text(tpx, tpy+dpy*1.5, 'NC='+str(ncgrd), 
#            horizontalalignment='left', fontsize=fontsz+2, color='r' )
#   ax.text(tpx, tpy+dpy*2.5, 'NPol= 2', 
#            horizontalalignment='left', fontsize=fontsz+2, color='r' )

#   ax.text(-tpx,-tpy-dpy*1.5, ritop,
#            horizontalalignment='right', fontsize=fontsz+4, color='k' )

#   ax.text(-tpx, tpy+dpy*0, 'Projection Pole',  
#            horizontalalignment='right', fontsize=fontsz+4, color='k' )
#   ax.text(-tpx, tpy+dpy*1.5, f'PLon={plon:8.2f}', 
#            horizontalalignment='right', fontsize=fontsz+2, color='r' )
#   ax.text(-tpx, tpy+dpy*2.5, f'PLat={plat:7.2f}', 
#            horizontalalignment='right', fontsize=fontsz+2, color='r' )
#   if( Arctic ):
#       ax.text( tpx, tpy+dpy*3.5, 'NA='+str(na), 
#            horizontalalignment='left', fontsize=fontsz+2, color='r' )
#       ax.text( tpx, tpy+dpy*4.5, 'NB='+str(nb), 
#            horizontalalignment='left', fontsize=fontsz+2, color='r' )

##  Save plot as ps file
#   print (" Save the smc grid local plot ... " )
#   plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', \
#               orientation=paprorn,papertype=paprtyp)

##  End of smcfulgrd plot function. ##

