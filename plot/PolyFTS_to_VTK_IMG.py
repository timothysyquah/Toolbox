#!/usr/bin/env python
import sys
import struct
import numpy as np
from os import path
import pdb
# Import PolyFTS IO library
sys.path.insert(0, path.split(sys.path[0])[0]) # Add search path for PolyFTSIO module, which is located in the parent directory
from PolyFTSIO import *

def get_PolyFTS_to_VTK_IdxMap(M,Nx):
    idx=np.zeros(M,dtype=np.uint64)
    Dim = len(Nx)
    if Dim == 1:
        for ix in range(Nx[0]):
            idx[ix] = ix
    elif Dim == 2:
        #looks good
        m=0
        for iy in range(Nx[1]):
          for ix in range(Nx[0]):
            idx[m] = ix*Nx[1] + iy
            m+=1
    elif Dim == 3:
        m=0
        for iz in range(Nx[2]):
          for iy in range(Nx[1]):
            for ix in range(Nx[0]):
                idx[m] = ix*Nx[1]*Nx[2] + iy*Nx[2] + iz
                #idx[m] = iz*Nx[0]*Nx[1] + iy*Nx[0] + ix
                m+=1
    return idx

# For command-line runs, build the relevant parser
import argparse as ap
parser = ap.ArgumentParser(description='Convert PolyFTS Field file to Legacy VTK format')
parser.add_argument('infile',metavar='inputfile',nargs='?',default='./density.bin',help='Input filename containing unformatted Field data')
# Parse the command-line arguments
args=parser.parse_args(sys.argv[1:])
# Generate the output file name automatically
outfile,ext = path.splitext(args.infile)
outfile = outfile + ".vtk"

# Read input file
fielddata = PolyFTSFieldReader()
fielddata.readFields(args.infile,False)

# Generate the mapping from PolyFTS (z fastest in 3D) to VTK (x fastest)
IdxMap = get_PolyFTS_to_VTK_IdxMap(fielddata.NPW,fielddata.griddim)
# Remap field samples from PolyFTS order to VTK order
#AllFields = fielddata.AllFields[:,IdxMap]
# For some reason the above one-liner doesn't work - data is permuted, but is only obvious in non-orthorhombic cells
AllFields = fielddata.AllFields[:,:]
AllFields[:,:] = AllFields[:,IdxMap]

# for debugging
#AllFields[0] = np.array([float(i)/M for i in range(M)])
#np.savetxt('idx.dat',IdxMap)
#np.savetxt('fields.dat',AllFields.transpose())
#pdb.set_trace()

# Legacy VTK file.
print "Outputting to Legacy VTK formatted file {}".format(outfile)
o = open(outfile,"w")

if fielddata.orthorhombic:
  o.write("# vtk DataFile Version 3.0\n")
  o.write("PolyFTS field data\n")
  o.write("ASCII\n")
  o.write("DATASET STRUCTURED_POINTS\n")
  spacing=np.zeros(fielddata.Dim)
  if fielddata.Dim == 1:
      spacing[0] = fielddata.AllMeshCoords[0][1]
      o.write("DIMENSIONS {} 1 1\n".format(*fielddata.griddim))
      o.write("ORIGIN 0\n")
      o.write("SPACING {} 0 0\n".format(*spacing))
  elif fielddata.Dim == 2:
      spacing[0] = fielddata.AllMeshCoords[0][fielddata.griddim[1]]
      spacing[1] = fielddata.AllMeshCoords[1][1]
      o.write("DIMENSIONS {} {} 1\n".format(*fielddata.griddim))
      o.write("ORIGIN 0 0 0\n")
      o.write("SPACING {} {} 0\n".format(*spacing))
  elif fielddata.Dim == 3:
      spacing[0] = fielddata.AllMeshCoords[0][fielddata.griddim[1]*fielddata.griddim[2]]
      spacing[1] = fielddata.AllMeshCoords[1][fielddata.griddim[2]]
      spacing[2] = fielddata.AllMeshCoords[2][1]
      o.write("DIMENSIONS {} {} {}\n".format(*fielddata.griddim))
      o.write("ORIGIN 0 0 0\n")
      o.write("SPACING {} {} {}\n".format(*spacing))
  o.write("POINT_DATA {0}\n".format(fielddata.NPW))
  for i in range(fielddata.nfields):
      o.write("SCALARS field{0} float 1\n".format(i))
      o.write("LOOKUP_TABLE default\n")
      np.savetxt(o, AllFields[i], fmt="%.11f")
else:
  o.write("# vtk DataFile Version 3.0\n")
  o.write("PolyFTS field data\n")
  o.write("ASCII\n")
  o.write("DATASET STRUCTURED_GRID\n")
  if fielddata.Dim == 1:
      o.write("DIMENSIONS {} 1 1\n".format(*fielddata.griddim))
  elif fielddata.Dim == 2:
      o.write("DIMENSIONS {} {} 1\n".format(*fielddata.griddim))
  elif fielddata.Dim == 3:
      o.write("DIMENSIONS {} {} {}\n".format(*fielddata.griddim))
  o.write("POINTS {} float\n".format(fielddata.NPW))
  # Remap the mesh coordinates to VTK order
  AllCoords = fielddata.AllMeshCoords[:,IdxMap]
  # Append the missing dimensions to the coordinates
  if fielddata.Dim < 3:
    AllCoords = np.append(AllCoords, np.zeros([1,AllCoords.shape[1]]), 0)
  if fielddata.Dim < 2:
    AllCoords = np.append(AllCoords, np.zeros([1,AllCoords.shape[1]]), 0)
  #np.savetxt('coords.dat',AllCoords.transpose()) # Debug
  np.savetxt(o, AllCoords.transpose(), fmt="%.11f")
  o.write("\nPOINT_DATA {0}\n".format(fielddata.NPW))
  for i in range(fielddata.nfields):
  #for i in range(1):
      # write as scalar
      #o.write("SCALARS field{0} float 1\n".format(i)) #STARTED TO SEG FAULT WHEN TRYING TO READ SCALARS
      #o.write("LOOKUP_TABLE default\n")
      #np.savetxt(o, AllFields[i], fmt="%14.11f")
      #np.savetxt(o, np.vstack([AllCoords,AllFields[i]]).transpose(), fmt="%14.11f")

      # write as vector
      o.write("VECTORS field{0} float\n".format(i)) #writing as vector fixed the seg fault
      tmp=np.zeros((3,fielddata.NPW))
      tmp[0,:] = fielddata.AllFields[i]
      np.savetxt(o, tmp.transpose(), fmt="%14.11f")
o.close()

