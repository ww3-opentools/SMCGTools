"""
##  The smclocal function plots a local region of a given smc grid.
##  cel is the full cell array in shape(nc, 5).
##                                JGLi04Mar2019
##  Modified to fill marked cells to be red.  JGLi03Sep2020
##  Move figure and ax setup to parent program.  JGLi18May2022
##  Change polar region parameter array for npl. JGLi16Sep2022
##
"""

def smcfulgrd(ax, cel, verts, ncels, colrs, config, Arctic=False, 
              grid='SMC1d', buoys='', panel=' ', ritop=' ',
              fontsz=9,  nmark=0):

##  Import relevant modules and functions

    import numpy as np
    import matplotlib.pyplot as plt

    from matplotlib.collections import PolyCollection

    from readtext import readtext
    from steromap import steromap
    from colrboxy import colrboxy
    from scale_depth import scale_depth

##  Degree to radian conversion parameter.
    d2rad=np.pi/180.0

##  Total cell number by cell array
    nc = cel.shape[0]

##  Local plot configuration parameters
    rdpols=config[0]
    sztpxy=config[1]
    rngsxy=config[2]
    clrbxy=config[3]
##  Local or sub-grid cell number if defined.
    if( len(config) >4 ):
        Arctic=True 
        ncabgm=config[4]
        na=int(ncabgm[0])
        nb=int(ncabgm[1])
        nbg=int(ncabgm[2])
        npl=int(ncabgm[3])
        ng= nc - na
        jm=np.max(cel[0:ng, 1])
        j3=cel[ng-1, 3]*3
        j1=cel[ng-1, 3]

##  Maximum mapping radius.
    radius=rdpols[0]
    pangle=rdpols[1]
    plon  =rdpols[2]
    plat  =rdpols[3]

##  Outline circle for pangle>=90.0
    ciran=np.arange(1081)*d2rad/3.0
    xcirc=radius*np.cos(ciran)
    ycirc=radius*np.sin(ciran)

##  Define depth to color index conversion parameters and marks.
##  Use only first 131 colors in colrs(0:255). 
    depth, factr, cstar, marks, ncstr, nclrm = scale_depth(nclrm=136)

##  Cell color is decided by its depth value.
    ndeps = np.zeros( (nc), dtype=np.int )
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
    ax.set_aspect('equal')
    ax.set_autoscale_on(False)
    ax.set_axis_off()
#   plt.subplots_adjust(left=0.0,bottom=0.0,right=1.0,top=1.0)

##  Add a outline circle if pangle >= 90.0
    if( pangle >= 90.0 ):
        ax.plot(xcirc, ycirc, 'b-', linewidth=0.3)

##  Create color array for this plot
    pcface = []
    pcedge = []
    kend = len(ncels)
    for k in range(kend):
        i = ncels[k] 
##  Mark the arctical map-east reference region by first ring  
##  of its boundary cells at jm-3
        if( Arctic and (cel[i,1] == jm-j3 or cel[i,1] == -jm+j3-j1) ):
            pcface.append(colrs(246))
            pcedge.append(colrs(246))
        elif( Arctic and (cel[i,1] == jm or cel[i,1] == -jm-j1) ):
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
    tpx=sztpxy[2] 
    tpy=sztpxy[3]
    dpy=0.6 

    ax.text(tpx,-tpy-dpy*1.5, panel,
             horizontalalignment='left', fontsize=fontsz+4, color='k' )

    ax.text(tpx, tpy+dpy*0, grid+' Grid',  
             horizontalalignment='left', fontsize=fontsz+4, color='k' )
    ax.text(tpx, tpy+dpy*1.5, 'NC='+str(nc), 
             horizontalalignment='left', fontsize=fontsz+2, color='r' )
    ax.text(tpx, tpy+dpy*2.5, 'NPol='+str(npl), 
             horizontalalignment='left', fontsize=fontsz+2, color='r' )

    ax.text(-tpx,-tpy-dpy*1.5, ritop,
             horizontalalignment='right', fontsize=fontsz+4, color='k' )

    ax.text(-tpx, tpy+dpy*0, 'Projection Pole',  
             horizontalalignment='right', fontsize=fontsz+4, color='k' )
    ax.text(-tpx, tpy+dpy*1.5, f'PLon={plon:8.2f}', 
             horizontalalignment='right', fontsize=fontsz+2, color='r' )
    ax.text(-tpx, tpy+dpy*2.5, f'PLat={plat:7.2f}', 
             horizontalalignment='right', fontsize=fontsz+2, color='r' )
    if( Arctic ):
        ax.text( tpx, tpy+dpy*3.5, 'NA='+str(na), 
             horizontalalignment='left', fontsize=fontsz+2, color='r' )
        ax.text(-tpx, tpy+dpy*3.5, 'NB='+str(nb), 
             horizontalalignment='right', fontsize=fontsz+2, color='r' )

#;  Overlay buoy sits on grid map if buoy file is provided.
    if( len(buoys) > 3 ):
        hdr, buoyll = readtext(buoys)
        nmbu=int(hdr[0])
        buoyids=buoyll[:,0].astype(str)
        buoylat=buoyll[:,1].astype(np.float)
        buoylon=buoyll[:,2].astype(np.float)

#; Convert slat slon to elat elon with given new pole
        elat,elon,sxc,syc = steromap(buoylat,buoylon,plat,plon,Pangl=pangle)

#; Mark buoy position on map
        print (' Selected buoys in this plot:')
        for i in range(nmbu):
            if( (elat[i] >= 25.0) and (rngsxy[0] < sxc[i] < rngsxy[1])
                              and (rngsxy[2] < syc[i] < rngsxy[3]) ):
                print (' {:6} {:8.3f} {:8.3f}'.format( buoyids[i], buoylat[i], buoylon[i] ))
                txtsz=int( abs(np.sin(elat[i]*d2rad)*12.0) )
                ax.text(sxc[i], syc[i], 'r',  fontsize=txtsz,
                     horizontalalignment='center', color='r' )
                ax.text(sxc[i], syc[i], '.',  fontsize=txtsz*2,
                     horizontalalignment='center', color='r' )

##  Save plot as ps file
#   print (" Save the smc grid local plot ... " )
#   plt.savefig(psfile, dpi=None,facecolor='w',edgecolor='w', \
#               orientation=paprorn,papertype=paprtyp)

##  End of smcfulgrd plot function. ##

