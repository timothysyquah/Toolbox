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
from mpl_toolkits.mplot3d import Axes3D
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



def plot_density(infile):
    #infile="density.bin"
    restart = False
    if restart:    
        com = np.load('com.npy')
        coords2d = np.load('coords2d.npy')
        fields2d = np.load('fields2d.npy')
        coords = np.load('coords.npy')
        fields = np.load('fields.npy')
    else:
        coords, fields = io.ReadBinFile(infile)

    __ndim = len(coords.shape) - 1
    __Nx = coords.shape[:__ndim]
    #__hvoxel = np.array([coords[1,0],coords[0,1]])
    __hvoxel = np.array([coords[1,0,0],coords[0,1,0],coords[0,0,1]])
    __hcell = __hvoxel * __Nx


    #nreplicates = (2,2)
    #repcoords,repfields = fieldtools.replicate_fields(coords,fields,nreplicates)

    zslice = 0.5*__hcell[2,2]
    zslice_index = int(0.5*__Nx[2])


    # use domain analyzer to get com and contours
    if not restart:
        domainanalyzer = DomainAnalyzer(coords,fields)
        domainanalyzer.setDensityThreshold(0.5)
        ndomains, com,area,vol,IQ = domainanalyzer.getDomainStats(plotMesh=False,add_periodic_domains=True)
        ##verts,faces,normals,values = domainanalyzer.meshAllDomains(plotfile='tmp.png')

        coords2d = coords[:,:,zslice_index,0:2]
        fields2d = fields[:,:,zslice_index,:]
        np.save('com.npy',com)
        np.save('coords2d.npy',coords2d)
        np.save('fields2d.npy',fields2d)
        np.save('coords.npy',coords)
        np.save('fields.npy',fields)

    domainanalyzer2d = DomainAnalyzer(coords2d,fields2d)
    contours = domainanalyzer2d.meshAllDomains()

    # add extra cells in -/+ x and -/+ y so that voronoi cells are in correct locations 
    com = np.vstack([com,com+[-1,0,0]*np.mat(__hcell), com+[0,-1,0]*np.mat(__hcell), com+[1,0,0]*np.mat(__hcell), com+[0,1,0]*np.mat(__hcell)])

    # using COM get voronoi cells using scipy
    vor = scipy.spatial.Voronoi(com)
    
    
    # now plot
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_axis_off() # turn axes off
    ax.set_aspect(1)
    ymax = np.max(coords2d[:,:,1])
    ymin = np.min(coords2d[:,:,1])
    ax.set_ylim(ymin,ymax)

    # reduce 3d voronoi cells to a 2d slice
    com2d=[]
    points2d=[]
    for i in range(com.shape[0]):
        if np.isclose(com[i,2],zslice):
            com2d.append([com[i,0],com[i,1]])
            points2d.append(i)
            #ax.text(com[i,0],com[i,1],f'{i}')
    com2d = np.array(com2d)

    #plot com2d
    #ax.plot(com2d[:,0],com2d[:,1], color='red',marker='x',ls='')

    for ridge_pair in vor.ridge_points:
        #if ridge_pair[0] in points2d or ridge_pair[1] in points2d:
        if True:
            myridge_verticies = vor.ridge_dict[tuple(ridge_pair)]

            #myverticies = np.array([vor.vertices[i] for i in myridge_verticies if i != -1 and np.all(np.abs(vor.vertices[i]) <100) ])
            myverticies = np.array([vor.vertices[i] for i in myridge_verticies if i != -1 ])
            
            zmin = np.min(myverticies[:,2])
            zmax = np.max(myverticies[:,2])
            #print(ridge_pair,zmin,zmax)
            tol=1e-2
            if np.isclose(zmax,zmin) or ((zmin+tol) < zslice and (zmax-tol) > zslice):

                #ax.plot(myverticies[:,0],myverticies[:,1],color='k',ls='',marker='o')
                #ax.plot(myverticies[:,0],myverticies[:,1],color='k')
                ax.plot(myverticies[:,0],myverticies[:,1],ls='-',color='k')
                ax.plot([myverticies[-1,0],myverticies[0,0]],[myverticies[-1,1],myverticies[0,1]],ls='-',color='k')
    #ax.legend()
  
    # plot contours
    for contour in contours:
        ax.plot(contour[:, 0], contour[:, 1], linewidth=2, color='k',ls='--')

    # plot density fields
    X = coords2d[:,:,0]
    Y = coords2d[:,:,1] 
    z = fields2d[:,:,0]
    ax.imshow(z.T,cmap=plt.cm.coolwarm, vmin=0.0, vmax=1.0,extent=[X.min(), X.max(), Y.min(), Y.max()],interpolation='bicubic', origin='lower')
    #im=ax.pcolormesh(xx,yy,repfields[:,:,0].T,cmap=mpl.cm.coolwarm,vmin=0,vmax=1)

    plt.tight_layout()
    plt.savefig('fig_domain_shape_SIGMA.png')
    plt.savefig('fig_domain_shape_SIGMA.pdf')
    plt.show()
    
    #fig = plt.figure()
    #ax = fig.add_subplot(111, projection='3d')
    ## -----------------------------------------------------------------
    ## plot 3d voronoi cells
    ## (dont delete, might be useful eventually)
    #for i,ridge_vertex in enumerate(vor.ridge_vertices):
    #    first = True
    #    for vert_index in ridge_vertex:
    #        if vert_index == -1: continue
    #        v = vor.vertices[vert_index]
    #        if not first and np.all(np.abs(v) <100):
    #            x = [vprev[0], v[0]]
    #            y = [vprev[1], v[1]]
    #            z = [vprev[2], v[2]]
    #            #if np.isclose(z[0],zslice) and np.isclose(z[1],zslice):
    #            ax.plot(x,y,z)
    #        else:
    #            vfirst = np.copy(v)
    #        first = False
    #        vprev = np.copy(v)
    #    # plot from first to last
    #    x = [vfirst[0], v[0]]
    #    y = [vfirst[1], v[1]]
    #    z = [vfirst[2], v[2]]
    #    ax.plot(x,y,z)
    ##end plot 3d voronoi cells
    ##-----------------------------------------------------------------
    ## plot center of mass
    #ax.plot(com[:,0],com[:,1],com[:,2], color='red',marker='x',ls='')
   


if __name__ == "__main__":
    
    #fig = plt.figure(figsize=(5,5),dpi=500)
    #ax = fig.add_subplot(111, projection='3d')
    infile = "density.bin"
    plot_density(infile)

    #fig = plt.figure(figsize=(6,4),dpi=500)

    #tau=0.925
    #nBAp=9
    #phiAs = np.arange(0.34,0.70,0.1)
    #phiAs = np.array([[0.3,0.4,0.5],[-1,0.6,0.7]])
    ##phiAs = [0.70,0.5,0.6,0.75]

    #axes = fig.subplots(2,3)
    #for i,axis in np.ndenumerate(axes):
    #    phiA = phiAs[i]
    #    if phiA == -1: 
    #        axes[i].set_axis_off()
    #    else:
    #        print(phiA)
    #        axes[i].set_title('$f_A = $' + f'{phiA:0.1f}')
    #        infile = f"nA1_nBAp{nBAp}_chiN40.0/tau{tau:0.3f}/phiA{phiA:0.3f}/HEXPhase/density.bin"
    #        plot_density(fig,axes[i],infile)

    ## draw arrow

    #fig.colorbar(im,ax=ax)
    #plt.tight_layout()
    #plt.savefig('fig_domain_shape_SIGMA.png')
    #plt.savefig('fig_domain_shape_SIGMA.pdf')
    #plt.show()
    




