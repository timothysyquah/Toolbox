#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Dec 22 23:24:21 2020

@author: tquah
"""
import os
import sys
#sys.path.append("/home/lequieu/Work/tools/lib")
mypath = os.path.realpath(sys.argv[0])#the absolute path to this script
libpath= '/'.join(mypath.split('/')[0:-2])+'/lib' # remove script name and domain_analysis directory
sys.path.append(libpath)
import numpy as np
import logging
import pdb
import viztools
import iotools 
import argparse
import re, glob
import matplotlib.pyplot as plt
from copy import deepcopy
#from scipy.interpolate import 
#writePolyFTSDatFile(outfilename, griddim, hcell, complexdata, kspacedata, fielddata, fielddataimpart=[])
import_path = '/home/tquah/TestDIR/readwrite_different_resolutions/GYR_fields.in'
export_path = '/home/tquah/TestDIR/readwrite_different_resolutions/fields.out'
field = iotools.ReadBinFile(import_path)