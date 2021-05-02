#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 28 00:49:50 2021

@author: tquah
"""
import sys
import os
#sys.path.append("/home/lequieu/Work/tools/lib/")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools/')
import numpy as np
from copy import deepcopy
from PolyFTSIO import *
from scipy.fftpack import ifftn,fftn


def getfield_stats(field):
    return np.mean(field),np.std(field)


#Doesnt seem to be doable!

os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/seeded/CL_TEST/investigate_amplitude')


# infilelist = 'amp_0.001/fields.bin'
# fielddata = PolyFTSFieldReader()
# fielddata.readFields(infilelist,False)

# AllFields = fielddata.AllFields[:,:]


# coord_local, fields_local = PolyFTSFieldReader
# fieldsk = fftn(fields_local)


# infilelist = 'amp_0.001/fields_k.bin'

# coord_local, fieldsk_local = io.ReadBinFile(infilelist)




