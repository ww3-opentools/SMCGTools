"""
##  The smcrglgrd function plots a local region of a given smc grid.
##  cel is the full cell array in shape(nc, 5).
##                                JGLi04Mar2019
##  Modified to fill marked cells to be red.  JGLi03Sep2020
##  Move figure and ax setup to parent program.  JGLi18May2022
##  Adapted for regular lat-lon box plot.     JGLi12Sep2022
##
"""

def smcrglgrd(ax, cel, verts, ncels, colrs, config, 
             grid='SMC1d', panel=' ', fontsz=9, nmark=0):

##  Import relevant modules and functions

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from readtext import readtext
    from colrboxy import colrboxy
    from scale_depth import scale_depth

##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

##  Total cell number by cell array
    nc = cel.shape[0]

##  Local plot configuration parameters
    sztpxy=config[0]
    rngsxy=config[1]
    clrbxy=config[2]

##  Arctic part cell number if defined.
    Arctic=False
    if( len(config) >3 ): 
        Arctic=True
        nabgpl=config[3]
        na= int(nabgpl[0])
        nb= int(nabgpl[1])
        nbg=int(nabgpl[2])
        npl=int(nabgpl[3])
        ng= nc - na 
        jm=np.max(cel[0:ng, 1])
        j3=cel[ng-1, 3]*3
        j1=cel[ng-1, 3]

##  Define depth to color index conversion parameters and marks.
##  Use only first 131 colors in colrs(0:255). 
    depth, factr, cstar, marks, ncstr, nclrm = scale_depth(nclrm=136)

##  Cell color is decided by its depth value.
    ndeps = np.zeros( (nc), dtype=int )
    for j in range(nc):
        if( cel[j,4] > -11 ): 
            ndeps[j] = ncstr + np.rint( (cstar-np.log10(cel[j,4]+11))*factr ).astype(np.int)
#   ndeps = ncstr + np.rint( (cstar-np.log10(cel[:,4]))*factr ).astype(np.int)
##  Sea level index
    ndep0 = ncstr + np.rint( (cstar-np.log10(11))*factr ).astype(np.int)
    print( " Sea level depth index is ", ndep0)

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
    kend = len(ncels)
    for k in range(kend):
        i = ncels[k] 
##  Mark the arctical map-east reference region by first ring  
##  of its boundary cells at jm-3
        if( Arctic and ( cel[i,1] == jm-j3 or cel[i,1] == -jm+j3-j1 ) ):
            pcface.append(colrs(246))
            pcedge.append(colrs(246))
        elif( Arctic and ( cel[i,1] == jm or cel[i,1] == -jm-j1 ) ):
            pcface.append(colrs(168))
            pcedge.append(colrs(168))
##  Fill last nmark cells to be red.  JGLi03Sep2020
        elif( k >= kend - nmark ):
            pcface.append(colrs(240))
            pcedge.append(colrs(ndeps[i]))
#           pcedge.append(colrs(250))
        else:
#           pcface.append(colrs(ndeps[i]))
            pcface.append(colrs(255))
            pcedge.append(colrs(ndeps[i]))

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
    for i in range(len(depth)):
        m = marks[i]
        if( dkx < dky ):    
            ax.text(xkeys[0]+1.15*dkx, ykeys[m], str(depth[i]),
                verticalalignment='center', fontsize=fontsz, color='b' )
        else:
            ax.text(xkeys[m], ykeys[0]+1.15*dky, str(depth[i]),
                horizontalalignment='center', fontsize=fontsz, color='b' )

    if( dkx < dky ):    
        ax.text(xkeys[0]+1.2*dkx, (ykeys[marks[2]]+ykeys[marks[1]])*0.5, 'Depth m',
                 rotation=-90,verticalalignment='center', fontsize=fontsz+2, color='k' )
    else:
        ax.text((xkeys[marks[1]]+xkeys[marks[2]])*0.5, ykeys[0]+1.2*dky, 'Depth m',
                 rotation=0,horizontalalignment='center', fontsize=fontsz+2, color='k' )

##  Put cell information inside plot
#   tpx=sztpxy[2] 
#   tpy=sztpxy[3]
#   dpy=3.0 
    tpx=(rngsxy[1] - rngsxy[0])*0.5
    tpy=rngsxy[2]+1.0

    ax.text(tpx, tpy, grid+' Grid  NC ='+str(nc),  
             horizontalalignment='center', fontsize=fontsz+3, color='k' )

#   ax.text(tpx,-tpy-dpy*1.5, panel,
#            horizontalalignment='left', fontsize=fontsz+4, color='k' )
#   ax.text(tpx, tpy+dpy*0, grid+' Grid',  
#            horizontalalignment='left', fontsize=fontsz+4, color='k' )
#   ax.text(tpx, tpy+dpy*1.5, 'NC='+str(nc), 
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

##  End of smcrglgrd plot function. ##

