#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 19:17:05 2020

@author: tquah
"""
####Note this needs a python2 env
#import relevent packages
from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np
import argparse
import os
from Fields_Check import fields_compare
#Function 
def return_fields(path):
    field = PolyFTSFieldReader()
    field.readFields(path,True)
    field_Re = field.AllFields
    field_Im = field.AllFieldsImPart
    hcell = field.hcell
    return field_Re,field_Im,hcell
#import path 
path = '/home/tquah/SEEDS/Polymer/SEEDS_2spec1red/HEXPhase/fields_k.bin'
export_path = '/home/tquah/SEEDS/Polymer/SEEDS_2spec1red/HEXPhase/fields.in'
NxNyNz = [32,32] #npw
#get real and imaginary fields you'll notice real fields don't do much
field_Re,field_Im,hcell = return_fields(path)
new_field = np.zeros(np.shape(field_Re))
new_field[:,0] = -field_Re[:,1]
new_field[:,1] = -field_Re[:,0]
writePolyFTSBinFile(export_path,NxNyNz, hcell, True, True, new_field, field_Im)