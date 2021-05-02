#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  9 01:17:11 2021

@author: tquah
"""

import os

import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools/')
import iotools as io
import numpy as np
from copy import deepcopy
from PolyFTSIO import PolyFTSFieldReader,PolyFTSFieldWriter
def getfield_stats(field):
    return np.mean(field),np.std(field)



multiplier = 0.01
os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/seeded/0.01_Amplitude_REFINE/fields')

listofdir = os.listdir()


for i in range(len(listofdir)):
    outfile = 'reduced_'+listofdir[i]
    field = PolyFTSFieldReader()
    field.readFields(listofdir[i],True)
    field_Re = field.AllFields*multiplier
    field_Im = field.AllFieldsImPart*multiplier
    hcell = field.hcell
    AllFields = field.AllFields[:,:]
    NxNyNz = np.array(field.griddim)    
    PolyFTSFieldWriter.writePolyFTSBinFile(outfile,NxNyNz, hcell, True, False, field_Re, field_Im)


path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/seeded/0.01_Amplitude_REFINE/reduced_fields/L_10.0/fields.bin'
field = PolyFTSFieldReader()
field.readFields(path,True)
field_Re = field.AllFields*multiplier
field_Im = field.AllFieldsImPart*multiplier
hcell = field.hcell
AllFields = field.AllFields[:,:]
NxNyNz = np.array(field.griddim)    
