"""
##  Draw local height plot with given polycollections and cell
##  array. Output as A3 ps file.        JGLi16May2022
##  Modified to use ax and heights.     JGLi23May2022
#;
"""

def smcfulplt(ax,hgts,swhs,verts,ncels,colrs,config,
        datx=' ', scheme=' ', rmser=' ', 
        panel=' ', grid=' ', fontsz=9):

##  Import relevant modules and functions

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from readtext import readtext
    from colrboxy import colrboxy
    from scale_linear import scale_linear

##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

##  Local plot configuration parameters
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]

##  Maximum mapping radius.
    radius=rdpols[0]

##  Set up linear height scale and marks with colrs.N.
##  colrs.N returns colrs' total number of colors 256.
#   hgts = [0,1,2,3,4,5,6]
    factor, marks, ncstr, nclrm = scale_linear(hgts, nclrm=colrs.N)

    nswh0= ncstr
    hgts0= float(hgts[0])
    swhmn= -ncstr/factor + hgts0
    swhmx= float(hgts[-1])

##  Some constant variables for plots.
    xprop={'xlim':rngsxy[0:2], 'xlabel':''}
    yprop={'ylim':rngsxy[2:4], 'ylabel':''}
    
##  Use ijk to count how many times to draw.
    ijk=0

    if True:
##  Work out max and min values, excluding missing data (-999.0)
        cmax = swhs.max()
        cmin = swhs[ swhs > -999.0 ].min()
        print ( ' swhs range %f, %f' % ( cmin,  cmax) )
        print ( ' Draw range %f, %f' % (swhmn, swhmx) )
        cmxs = 'Hmx=%6.3f' % cmax
        cmns = 'Hmn=%7.4f' % cmin

##  Reset missing values (-999.0) or any value < swhmn to be swhmn 
        swhs[ swhs < swhmn ] = swhmn

##  Trim large values into plot range if any
        swhs[ swhs > swhmx ] = swhmx

##  Convert swhs with logarithm scale.
        icnf = ncstr+ np.rint( factor*(swhs - hgts0) )
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
        xkeys, ykeys, clrply = colrboxy(clrbxy, colrs, marks)
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

#       if( dkx < dky ):
#           ax.text(xkeys[0]+2.0*dkx, ykeys[marks[3]], 'H m',
#               rotation=90,verticalalignment='center', fontsize=15, color='k' )
#       else:
#           ax.text(xkeys[marks[3]], ykeys[0]+2.0*dky, 'H m',
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

        ax.text(-tpx, tpy+dpy*0., 'Height m', 
             horizontalalignment='right', fontsize=fontsz, color='b' )
        ax.text(-tpx, tpy+dpy*1., grid, 
             horizontalalignment='right', fontsize=fontsz, color='r' )
        ax.text(-tpx, tpy+dpy*2.2, scheme, 
             horizontalalignment='right', fontsize=fontsz+2, color='k' )


##  End of smcfulplt

