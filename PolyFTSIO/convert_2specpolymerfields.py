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
import os
import glob
#Function 
def return_fields(path):
    field = PolyFTSFieldReader()
    field.readFields(path,True)
    field_Re = field.AllFields
    field_Im = field.AllFieldsImPart
    hcell = field.hcell
    npw = field.griddim
    return field_Re,field_Im,hcell,npw
#import path 
#path = '/home/tquah/EXPORT_KNOT/SEEDS_2spec1red'
#os.chdir(path)
Phases = glob.glob('*Phase')
#Phases = ['A15Phase']


#export_path = '/home/tquah/SEEDS/Polymer/SEEDS_2spec1red/HEXPhase/fields.in'
#get real and imaginary fields you'll notice real fields don't do much

WDIR = os.getcwd()
for phase in Phases:
    os.chdir(phase)
    field_Re,field_Im,hcell,npw = return_fields('fields_save.in')
    new_field = np.zeros(np.shape(field_Re))
    new_field[0,:] = -field_Re[1,:]
    new_field[1,:] = -field_Re[0,:]
    new_field_im = np.zeros(np.shape(field_Im))
    new_field_im[0,:] = -field_Im[1,:]
    new_field_im[1,:] = -field_Im[0,:]

    writePolyFTSBinFile('fields.in',npw, hcell, True, True, new_field, new_field_im)
    os.chdir(WDIR)
    