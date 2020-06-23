#!/usr/bin/env python3

import numpy as np
import pdb
import sys
import os
sys.path.append("/home/lequieu/tools/lib/")
#mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
#libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
#sys.path.append(libpath)

import iotools as io

infile = 'density.bin'
#infile = 'fields.in'
coords, fields = io.ReadBinFile(infile)
dim = coords.shape[-1]
Nx = coords.shape[0:dim]
nfields = fields.shape[-1]

downsampling_freq = 2
if downsampling_freq > 1:
    Nx = tuple(np.divide(Nx,downsampling_freq).astype(int))
    coords = coords[::downsampling_freq,::downsampling_freq,::downsampling_freq]
    fields = fields[::downsampling_freq,::downsampling_freq,::downsampling_freq]


maxspecies = np.zeros(Nx)

for (i,j,k) in np.ndindex(Nx):
    argmax = np.argmax(fields[i,j,k])
    maxspecies[i,j,k] = argmax


import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
#ax.set_aspect(1)
xx = coords[:,:,0]
yy = coords[:,:,1]
zz = coords[:,:,2]
#im=ax.pcolormesh(xx,yy,maxspecies.T)

# setup colors of species A and C
colors = np.empty(Nx,dtype=object)
colors[maxspecies == 0] = 'red'
colors[maxspecies == 4] = 'blue'

# only draw voxels for A and C
voxels = np.zeros(Nx,dtype=np.bool)
voxels[maxspecies == 0] = True
voxels[maxspecies == 4] = True

ax.voxels(voxels,facecolors=colors,edgecolor='k')
#ax.scatter(xx,yy,zz,c=maxspecies.T==1)

#im=ax.pcolormesh(xx,yy,fields[:,:,0].T)
#fig.colorbar(im,ax=ax)
plt.show()


