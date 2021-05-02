#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 20 20:48:16 2021

@author: tquah
"""
import sys
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools/PolyFTSIO')
from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np
import argparse
from copy import deepcopy
if __name__ == "__main__":
    #input parameters
    parser = argparse.ArgumentParser(description='Weakens PolyFTS fields_k for use in CL simulations')
    parser.add_argument('-i', '--infile', action='store', default='fields_k.bin',help='field input in kspace')
    parser.add_argument('-o', '--outfile', action='store', default='fields.out',help='field output in kspace')
    parser.add_argument('-m', '--multiplier', action='store', default=0.01,type=float,help='multiplier to weaken fields')
    args = parser.parse_args()
    
    #field paramters read in
    field = PolyFTSFieldReader()
    field.readFields(args.infile,True)
    field_Re = field.AllFields
    field_Im = field.AllFieldsImPart
    hcell = field.hcell
    print('Assumes Cubic Cell')
    npw = field.NPW
    NxNyNz = [np.int(np.round(field.NPW**(1/3))),np.int(np.round(field.NPW**(1/3))),np.int(np.round(field.NPW**(1/3)))]  
    #multiplier fields
    field_Re_local = deepcopy(field_Re)
    field_Im_local = deepcopy(field_Im)
    field_Re_local = field_Re_local*args.multiplier
    field_Im_local = field_Im_local*args.multiplier
    #write output
    writePolyFTSBinFile(args.outfile,NxNyNz, hcell, True, True, field_Re_local, field_Im_local)
