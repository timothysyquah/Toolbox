#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 20 19:17:05 2020

@author: tquah
"""

from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np
import argparse
import os
from Fields_Check import fields_compare

def return_fields(path):
    field = PolyFTSFieldReader()
    field.readFields(path,True)
    field_Re = field.AllFields
    field_Im = field.AllFieldsImPart
    return field_Re,field_Im

#Nref = 100
#Nref_desire = 1 
#path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt2Spec/spheresSigma_P4_2overmnm_Aminor_fields.in'
#path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt2Spec/networkO70_Fdddconv_Aminor_fields.in'
#path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt2Spec/networkO70_Fdddconv_Aminor_fields.in'
path1 = '/home/tquah/SEEDS/Polymer/DifferenceBetweenSeedTypes/spheresBCCprim_Im-3m_Aminor_fields_2spec.in'
path2 = '/home/tquah/SEEDS/Polymer/DifferenceBetweenSeedTypes/spheresBCCprim_Im-3m_Aminor_fields.in'

#export_path = '/home/tquah/SEEDS/Polymer/BCCPhase/fields.in'



field1_Re,field1_Im = return_fields(path1)
field2_Re,field2_Im = return_fields(path2)


field_diff = field1_Re[1,:]+field2_Re[0,:]
field_diff = field1_Re[0,:]+field2_Re[0,:]

#field_adjust = Nref/Nref_desire
#hcell_adjust = np.sqrt(Nref/Nref_desire)

    
#field_Re = field.AllFields/field_adjust
#field_Im = field.AllFieldsImPart/field_adjust
#hcell = field.hcell/hcell_adjust
#NxNyNz = [16,16,16]
#writePolyFTSBinFile(export_path,NxNyNz, hcell, True, True, field_Re, field_Im)
#
#fieldtest = PolyFTSFieldReader()
#fieldtest.readFields(export_path,True)


#
