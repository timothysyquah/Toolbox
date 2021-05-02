#!/usr/bin/env python3

import numpy as np
import pdb
import sys
import os
sys.path.append("/home/tquah/toolbox_github/lib/")
from scipy.fftpack import fftn, fftshift,ifftn

#mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
#libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
#sys.path.append(libpath)
import matplotlib.pyplot as plt
import iotools as io
from domaintools import DomainAnalyzer
from skimage import feature,measure
from scipy import ndimage as ndi
from skimage import filters
from skimage.data import camera
from skimage.util import compare_images
from scipy.stats import gaussian_kde
from scipy.signal import peak_widths
from scipy.interpolate import CubicSpline
from scipy.signal import chirp, find_peaks, peak_widths
from scipy.stats import linregress
from mpl_toolkits.mplot3d import axes3d
from scipy.linalg import lstsq
from skspatial.objects import Points, Plane,Line
from skspatial.plotting import plot_3d

infile =  '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/average100/c_6.2/DensityOperator.dat'
# infile = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/REFINE/average100/c_8.0_SCFT/density.bin'
#infile = 'fields.in'
path = '/home/tquah/Figures/'
#load image
coords, fields = io.ReadDatFile(infile)
# coords, fields = io.ReadBinFile(infile)

dim = coords.shape[-1]
Nx = coords.shape[0:dim]
nfields = fields.shape[-1]
plt.close('all')
image = fields[:,:,:,0]
#get only edeges of 3d
edge_sobel = filters.sobel(image)
plt.figure()
plt.imshow(fields[:,:,:,0][0])
plt.colorbar()
name = '2D_fields.png'
export_path = os.path.join(path,name)
plt.savefig(export_path,dpi=300)



plt.figure()
plt.imshow(edge_sobel[0])
plt.colorbar()

name = '2D_edge.png'
export_path = os.path.join(path,name)
plt.savefig(export_path,dpi=300)


a = 2.0
model = gaussian_kde(np.ravel(edge_sobel))
x = np.linspace(0,1.0,1000)
y = model.pdf(x)
# dx = x[1]-x[0]
# dy = np.diff(y)/dx
# # peak_widths(y,)

# sdiff = np.sign(dy)
# loc = np.where(sdiff[:-1] > sdiff[1:])
#isolate to peaks
peaks, _ = find_peaks(y)
results_half = peak_widths(y, peaks, rel_height=0.1)
width,heights,left,right = peak_widths(y, peaks, rel_height=1.0)

location = [int(left[-1]),int(right[-1])]



plt.figure()
plt.hist(np.ravel(edge_sobel),bins=20,density=True)
plt.plot(x,y)
plt.plot(x[location],y[location],'ok')
plt.ylabel('Count')
plt.xlabel('Intensity')
name = 'hist.png'
plt.tight_layout()
export_path = os.path.join(path,name)
plt.savefig(export_path,dpi=300)



loc =np.where(edge_sobel<x[location][0])
edge_sobel[loc] = 0
loc =np.where(edge_sobel>x[location][0])
edge_sobel[loc] = 1

plt.figure()
plt.imshow(edge_sobel[:,:,20])


labels = measure.label(edge_sobel)
types = np.unique(labels)
# fig, axes = plt.subplots(ncols=3, sharex=True, sharey=True,
#                           figsize=(8, 4))




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

# Grab some test data.

# Plot a basic wireframe.



domainx = np.arange(0,39+1e-6,1)
planes = []
points = []
normal = []
projection = []
for i in range(1,len(types)):
    
    loc = np.where(labels==types[i])
    ax.scatter(loc[0], loc[1], loc[2],label =float(i))
    print(len(loc[0]))
    see = np.vstack(loc).transpose()
    plane = Plane.best_fit(see)
    if len(loc[0]) >1000:
        normalvec = plane.normal
        point = 1e8*np.array(normalvec)

    else:
        plane = Plane(plane.point,normalvec)

    planes.append(plane)
    points.append(np.array(plane.point))
    normal.append(np.array(plane.normal))
    projection.append(np.array(plane.project_point(point))[0])
    xx, yy = np.meshgrid(range(100), range(100))
    pp = np.ones((1,3))*points[i-1]
    d = -pp.dot(normalvec.transpose())
    z = (-normalvec[0] * xx - normalvec[1] * yy - d) * 1. /normalvec[2]
    ax.plot_surface(xx, yy, z, alpha=0.2)



points_array = np.vstack(points)
normal_vector = np.vstack(normal)
index = np.argsort(np.sqrt(np.sum(points_array**2,axis=1)))





ax.legend()
name = 'planes.png'
export_path = os.path.join(path,name)
plt.savefig(export_path,dpi=300)

dappend = []
for j in range(1,len(planes)):
    i = index[j]
    normalline = Line(planes[i-1].point,direction = planes[0].normal) 
    p0 = np.array(planes[i-1].point)
    intersection_point = np.array(planes[i].intersect_line(normalline))
    distance = np.sqrt(np.sum((intersection_point-p0)**2))
    dappend.append(distance)
    
print(np.mean(np.array(dappend)[1:])*4)


