#!/usr/bin/env python3

import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)

import fieldtools 
import iotools as io
import numpy as np
from domaintools import DomainAnalyzer

import pdb
import scipy.spatial #import Voronoi, voronoi_plot_2d
import matplotlib.pyplot as plt
import matplotlib as mpl


'''Simple example showing use of fieldstools.replicatefields()

'''
def contour_perimeter(contour):
    '''calculate perimeter of contour by suming up the line-segment lengths
        (coppied from domaintools)
    '''
    #assert (np.all(contour[0] == contour[-1])), "Contour must be closed! (1st point == last point)"

    #TODO vectorize this for loop
    p = 0.0
    n=contour.shape[0]
    for i in range(n-1):
       v = contour[i+1] - contour[i] 
       p += np.sqrt(np.square(v).sum())
    return p



def plot_density(fix,ax,infile):
    #infile="density.bin"
    coords, fields = io.ReadBinFile(infile)
    __ndim = len(coords.shape) - 1
    __Nx = coords.shape[:__ndim]
    __hvoxel = np.array([coords[1,0],coords[0,1]])
    __hcell = __hvoxel * __Nx

    nreplicates = (2,2)
    
    if (__ndim != 2):
        raise  NotImplementedError('implement me!')

    repcoords,repfields = fieldtools.replicate_fields(coords,fields,nreplicates)

    # use domain analyzer to get com and contours
    domainanalyzer = DomainAnalyzer(repcoords,repfields)
    domainanalyzer.setDensityThreshold(0.5)
    ndomains, com,area,vol,IQ = domainanalyzer.getDomainStats(plotMesh=True,add_periodic_domains=True)
    contours = domainanalyzer.meshAllDomains()

  
    perimeter = np.zeros(len(contours))
    for i,contour in enumerate(contours):
        perimeter[i] = contour_perimeter(contour)

    # remove the tiny domains from contour
    newcontours = []
    maxperimeter = np.max(perimeter)
    minperimeter = np.min(perimeter)
    for i,contour in enumerate(contours):
        include=True
        if np.all(contour[0] == contour[-1]): # if closed contour
            if (np.abs(perimeter[i]-minperimeter) < 0.1): # and equal to the min
                include=False
        else: 
            include = False # this removes all open contours
        if include:
            newcontours.append(contours[i])
    contours = newcontours 
        



    # sometimes if phiA is big, little domains start to pop up
    # remove these by compating the max and min volumes
    maxvol = np.max(vol)
    minvol = np.min(vol)
    if not np.isclose(minvol,maxvol,rtol=1e-1):
        index = np.isclose(vol,maxvol,rtol=1e-1) 
        com = com[index]        

    # using COM get voronoi cells using scipy
    vor = scipy.spatial.Voronoi(com)
    
    ymax = np.max(coords[:,:,1])
    ymin = np.min(coords[:,:,1])
    
    # now plot
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off() # turn axes off
    ax.set_aspect(1)
    ax.set_ylim(ymin,ymax)

    # plot contours
    for contour in contours:
        ax.plot(contour[:, 0], contour[:, 1], linewidth=2, color='k',ls='--')
   
    # plot voronoi cells
    for i,r in enumerate(vor.ridge_vertices):
        #if r[0] == 0 or r[1] == 0:
        #    continue # first point seems to be bad
        v1 = vor.vertices[r[0]]
        v2 = vor.vertices[r[1]]
        if r[0] == -1:
            continue # for some reason if there's a -1 these ridges are incorrect
        if np.abs(v1[1]) > 3 or np.abs(v2[1]) > 3 :
            continue # manual cutoff to only draw center voronoi cell
        x = [v1[0], v2[0] ]
        y = [v1[1], v2[1] ]
        ax.plot(x,y, color='k',ls=':',marker='')
        #print("i: ",i," ",r[0],r[1],v1,v2,x,y)
    
    # plot center of mass
    #ax.plot(com[:,0],com[:,1], color='red',marker='x',ls='')

    # densities
    xx = repcoords[:,:,0]
    yy = repcoords[:,:,1]
    im=ax.pcolormesh(xx,yy,repfields[:,:,0].T,cmap=mpl.cm.coolwarm,vmin=0,vmax=1)
    #fig.colorbar(im,ax=ax)


if __name__ == "__main__":
    
    #fig = plt.figure(figsize=(5,5),dpi=500)
    #ax = fig.subplots()
    #infile = "density.bin"
    #plot_density(fig,ax,infile)

    fig = plt.figure(figsize=(6,4),dpi=500)

    tau=0.925
    nBAp=9
    phiAs = np.arange(0.34,0.70,0.1)
    phiAs = np.array([[0.3,0.4,0.5],[-1,0.6,0.7]])
    #phiAs = [0.70,0.5,0.6,0.75]

    axes = fig.subplots(2,3)
    for i,axis in np.ndenumerate(axes):
        phiA = phiAs[i]
        if phiA == -1: 
            axes[i].set_axis_off()
        else:
            print(phiA)
            axes[i].set_title('$f_A = $' + f'{phiA:0.1f}')
            infile = f"nA1_nBAp{nBAp}_chiN40.0/tau{tau:0.3f}/phiA{phiA:0.3f}/HEXPhase/density.bin"
            plot_density(fig,axes[i],infile)

    # draw arrow

    #fig.colorbar(im,ax=ax)
    plt.tight_layout()
    plt.savefig('fig_domain_shape_HEX.png')
    plt.savefig('fig_domain_shape_HEX.pdf')
    #plt.show()
    




