#!/usr/bin/env python3

import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)
import iotools as io

import h5py
import numpy as np
import pdb


# TODO: setup argparse
import argparse as ap
parser = ap.ArgumentParser(description='This scrips converts a PolyFTS density file into a hdf5 file \
            that can be read by the MPB code to input the dielectric as a function of space')
parser.add_argument('--epsilons',default=[],nargs="+",help='epsilon of each species in the system')
parser.add_argument('--input',default='',type=str,help='Input fields file from PolyFTS. Can be binary or dat')
parser.add_argument('--output',default='epsilon.h5',type=str,help='Output h5 file to be read by MPB code')
# Parse the command-line arguments
args=parser.parse_args()

# setup and parse input file name
# use iotools to load density.bin to get coords and fields
if args.input != '':
  infile = args.input
  suffix = infile.split('.')[-1]
  if suffix == 'bin':
    coords, fields = io.ReadBinFile(infile)
    fields_have_imag = True
  elif suffix == 'dat':
    coords, fields = io.ReadDatFile(infile)
    fields_have_imag = False
  else:
    raise RuntimeError(f"invalid suffix of args.input: {suffix}")
else: # args.input not specified, try density.bin and density.dat
   
  if os.path.isfile('density.bin'):
    coords, fields = io.ReadBinFile('density.bin')
    fields_have_imag = True
  elif os.path.isfile('density.dat'):
    coords, fields = io.ReadDatFile('density.dat')
    fields_have_imag = False
  else:
    raise RuntimeError("Default inputs density.bin and density.dat not found")

#h5_filename='epsilon.h5'
h5_filename=args.output

# this is so we get the right field index
if fields_have_imag:
  factor=2
else :
  factor = 1

ndim = len(coords.shape) - 1
Nx = coords.shape[:ndim]

nfields = fields.shape[-1]
nspecies = int(nfields / factor) #TODO is this right?

# set epsilon of each species
if args.epsilons == []: # if empty use defaults
  if nspecies == 2:
    eps_of_species = np.array([13.0, 1.0])
  elif nspecies == 3:
    eps_of_species = np.array([13.0, 1.0, 6.0])
else: # else use values from comand args
  assert(len(args.epsilons) == nspecies)
  eps_of_species = [float(i) for i in args.epsilons]

# create epsilon(r) field
epsilon = np.zeros(Nx)
for i in range(nspecies):
  if ndim == 2:
    epsilon += eps_of_species[i] * fields[:,:,factor*i] #TODO confirm if this should be x2 (if imag vs real)
  elif ndim == 3:
    epsilon += eps_of_species[i] * fields[:,:,:,factor*i] #TODO confirm if this should be x2 (if imag vs real)

# write to h5 file
f = h5py.File(h5_filename, "w")
# set epsilon as data of "dataset"
dset = f.create_dataset("dataset", Nx, dtype='f', data=epsilon)
f.close()


#f2 = h5py.File(h5_filename, 'r')
#pdb.set_trace()
