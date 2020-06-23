#!/usr/bin/env python3

import sys
import os
sys.path.append("/home/mkamata/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)


import fieldtools 
import iotools as io
import numpy as np
import pdb


import argparse as ap
parser = ap.ArgumentParser(description='Script to quickly plot a 2d density of a triblock')
parser.add_argument('--nreplicates',default=[1,1],nargs=2,help='number of times to replicate box in each dimension')
parser.add_argument('--input','-i',default='density.bin',type=str,help='input density fields')
parser.add_argument('--output','-o',default='density.png',type=str,help='Output h5 file to be read by MPB code')
parser.add_argument('--quiet','-q',default=False,action='store_true',help='should I show interactive plot window?')
parser.add_argument('--smooth',default=False,action='store_true',help='should I smooth the densities')
# Parse the command-line arguments
args=parser.parse_args()



infile = args.input
#infile = 'fields.in'
nreplicates = tuple([int(i) for i in args.nreplicates])
interpolate_flag=args.smooth
pngfilename=args.output
quiet_flag=args.quiet

if quiet_flag and (not pngfilename):
  raise RuntimeError("Quiet flag set and no pngfilename specified. No output will be generated")

coords, fields = io.ReadBinFile(infile)
repcoords,repfields = fieldtools.replicate_fields(coords,fields,nreplicates)
dim = repcoords.shape[-1]
Nx = repcoords.shape[0:dim]
#pdb.set_trace()

assert(dim==2), "Script only works for 2d plots"

maxspecies = np.zeros(Nx)

for (i,j) in np.ndindex(Nx):
    argmax = np.argmax(repfields[i,j])
    maxspecies[i,j] = argmax/2.0

#multiply maxspecies by 2 and add 1 to get different blocks to show up as Blue, green red
maxspecies = 2.0*maxspecies + 1.0

import matplotlib.pyplot as plt
import matplotlib as mpl

fig, ax = plt.subplots()
ax.set_aspect(1)
xx = repcoords[:,:,0]
yy = repcoords[:,:,1]
if not interpolate_flag:
  # this uses pcolormesh
  im=ax.pcolormesh(xx,yy,maxspecies.T,cmap=mpl.cm.Paired,vmin=0,vmax=11) # set vmax to manually pick out color from "Paired" colormap
else:
  # this uses imshow, which permits interpolation
  # can break on braid
  colors = [(50/256.,160/256.,44/256.0,i) for i in np.linspace(0,1,3)] # green 
  cmap = mpl.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=1)
  im=ax.imshow(repfields[:,:,2],cmap=cmap,vmin=0,vmax=1)

  colors = [(227/256.,26/256.,28/256.0,i) for i in np.linspace(0,1,3)] # red
  cmap = mpl.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=2)
  im=ax.imshow(repfields[:,:,0],interpolation='gaussian',cmap=cmap,vmin=0,vmax=1)

  colors = [(31/256.,120/256.,180/256.0,i) for i in np.linspace(0,1,3)] # blue
  cmap = mpl.colors.LinearSegmentedColormap.from_list('mycmap', colors, N=2)
  im=ax.imshow(repfields[:,:,4],interpolation='gaussian',cmap=cmap,vmin=0,vmax=1)


#im=ax.pcolormesh(xx,yy,repfields[:,:,0].T)
#fig.colorbar(im,ax=ax)
plt.axis('off')
if pngfilename:
  plt.savefig(pngfilename,dpi=500)
if not quiet_flag:
  plt.show()

