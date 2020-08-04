#!/usr/bin/env python

from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np

Nx=8
NPW=Nx*Nx*Nx
dummyfielddata = np.zeros((2,NPW)) # Dummy data for two fields
dummyfielddataimpart = np.zeros((2,NPW)) # Dummy data for two fields
# Distinguish components by magnitude and sign, for testing
dummyfielddata[0] = np.random.rand(NPW)
dummyfielddata[1] = 10.*np.random.rand(NPW)
dummyfielddataimpart[0] = -np.random.rand(NPW)
dummyfielddataimpart[1] = -10.*np.random.rand(NPW)
writePolyFTSBinFile("test.bin", [Nx,Nx,Nx], [[4.0,0.,0.],[0.,3.,0.],[1.,0.,5.]], True, False, dummyfielddata, dummyfielddataimpart)
fields = PolyFTSFieldReader()
try:
    fields.readFields("test.bin",True)
except:
    print "*** Error during field reading ***"

print "Version = ",fields.version
print " Dim = ",fields.Dim
print " Cell = ",fields.hcell
print " # fields = ",fields.nfields
print " grid = ",fields.griddim
print " Total NPW = ",fields.NPW
print "\n\n MESH:"
print np.transpose(fields.AllMeshCoords)
#for i in range(fields.Dim):
#    print "direction ",i," = ",fields.AllMeshCoords[i]
print " Orthorhombic? ",fields.orthorhombic
print " Mesh spacing = ",fields.spacing

print "\n Fields by index:"
for i in range(fields.nfields):
    print fields.AllFields[i]
    print fields.AllFieldsImPart[i]
    print "Num points = ",len(fields.AllFields[i])
# TODO: turn this into a proper unit test - store out data in variables, write it out, read it back in and compare.
# Then clean up files.

fpre = fields.AllFields
fafter = dummyfielddata



fnorm1 = fpre[0,:]/np.linalg.norm(fpre[0,:])
fnorm2 = fafter[0,:]/np.linalg.norm(fafter[0,:])
check = np.dot(fnorm1,fnorm2.transpose())