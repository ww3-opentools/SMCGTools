"""
##  Draw local height difference with given polycollections and cell
##  array. Output as A3 ps file.        JGLi16May2022
##  Modified to use ax and heights.     JGLi23May2022
##  Modified to draw Hw difference.     JGLi08Sep2022
#;
"""

def smcfuldif(ax,whtdif,swhs,verts,ncels,colrs,config,
        datx=' ', scheme=' ', rmser=' ', 
        panel=' ', grid=' ', fontsz=9):

##  Import relevant modules and functions

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from readtext import readtext
    from colrboxy import colrboxy
    from colrboxy import colrboxy
    from scale_diff import scale_diff

##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

##  Local plot configuration parameters
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]

##  Maximum mapping radius.
    radius=rdpols[0]

##  Set up difference scale and marks with given colrs.
#   whtdif=[-90, -60, -30, 0, 30, 60, 90]
##  colrs.N returns colrs' total number of colors 256.
    bottom, ceilng, factor, nmarks, ncstr, nclrm = scale_diff(whtdif, nclrm=colrs.N)
    print(' bottom ceilng factor ncstr nclrm = %f %f %f %d %d' % (bottom, ceilng, factor, ncstr, nclrm)) 

##  Some constant variables for plots.
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    
##  Use ijk to count how many times to draw.
    ijk=0

    if True:
##  Work out max and min values, excluding missing data fmdi
        fmdi = -1.0E4
        cmax = swhs.max()
        cmin = swhs[ swhs > fmdi ].min()
        print ( ' swh range %f, %f' % (cmin, cmax) )
        cmxs = 'Difmx = %6.2f m' % cmax
        cmns = 'Difmn =%7.2f m' % cmin
#       cmns = 'Difmn = %10.3E' % cmin

##  Reset minimum values to be bottom.
        swhs[ swhs < bottom ] = bottom

##  Trim large values into plot range if any
        swhs[ swhs > ceilng ] = ceilng

##  Check range of modified difference field.
        cmax = swhs.max()
        cmin = swhs[ swhs > fmdi ].min()
        print ( ' Modified range %f, %f' % (cmin, cmax) )

##  Bottom value use a special color.
        nswh0= ncstr
##  Convert swhs with logarithm scale.
        swhp = np.array( swhs, dtype=float )
        icnf = ncstr + np.rint( factor*(swhp - bottom) )
        nswh = np.array( icnf, dtype=int )

        print (" Drawing panel "+panel )

##  Set up first subplot and axis for northern hemisphere

        ax.set(**xprop)
        ax.set(**yprop)
        ax.set_aspect('equal')
        ax.set_autoscale_on(False)
        ax.set_axis_off()
#       plt.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)

## Prepare PolyCollection for this plot.
        polynorth = PolyCollection(verts)

## Create color array for this plot
        pcface = []
        pcedge = []
        for i in ncels:
            if( nswh[i] != nswh0 ): 
                pcface.append(colrs(nswh[i]))
                pcedge.append(colrs(nswh[i]))
            else:
                pcface.append(colrs(255)) 
                pcedge.append(colrs(0)) 

#   smcpoly.set_color(pcolr)    ## This line defines both edge and face color.
        polynorth.set_facecolor(pcface)
        polynorth.set_edgecolor(pcedge)
        polynorth.set_linewidth( 0.2 )
#       print (" Drawing north hemisphere cells ... ")
        ax.add_collection(polynorth)  

##  Draw colorbar inside plot.
        xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, nmarks)
        ax.add_collection(clrply)
        dkx=clrbxy[2]; dky=clrbxy[3]
        for i in range(len(whtdif)):
            m = nmarks[i]
            if( dkx < dky ): 
                plt.text(xkeys[0]+1.15*dkx, ykeys[m], str(whtdif[i]),
                    verticalalignment='center', fontsize=11, color='b' )
            else:
                plt.text(xkeys[m], ykeys[0]+1.15*dky, str(whtdif[i]),
                    horizontalalignment='center', fontsize=11, color='b' )

#       if( dkx < dky ):
#           plt.text(xkeys[0]+2.0*dkx, ykeys[nmarks[3]], 'WHtDiff m',
#               rotation=90,verticalalignment='center', fontsize=15, color='k' )
#       else:
#           plt.text(xkeys[nmarks[3]], ykeys[0]+2.0*dky, 'WHtDiff m',
#               rotation=0,horizontalalignment='center', fontsize=15, color='k' )

##  Put cell information inside plot
        tpx= sztpxy[2]
        tpy= sztpxy[3]
        dpy= 0.6
        ax.text(tpx,-tpy-dpy*1.5, panel,
             horizontalalignment='left', fontsize=fontsz+4, color='k' )

        ax.text(tpx, tpy+dpy*0., cmns,
             horizontalalignment='left', fontsize=fontsz, color='b' )
        ax.text(tpx, tpy+dpy*1., cmxs,
             horizontalalignment='left', fontsize=fontsz, color='r' )
        ax.text(tpx, tpy+dpy*2.2, datx,
             horizontalalignment='left', fontsize=fontsz+2, color='k' )

        ax.text(-tpx,-tpy-dpy*1.5, rmser,
             horizontalalignment='right', fontsize=fontsz+4, color='k' )

        ax.text(-tpx, tpy+dpy*0., 'WHtDiff m', 
             horizontalalignment='right', fontsize=fontsz, color='b' )
        ax.text(-tpx, tpy+dpy*1., grid, 
             horizontalalignment='right', fontsize=fontsz, color='r' )
        ax.text(-tpx, tpy+dpy*2.2, scheme, 
             horizontalalignment='right', fontsize=fontsz+2, color='k' )


##  End of smcfuldif


