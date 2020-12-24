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
Nref = 100
Nref_desire = 1 
#path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt2Spec/spheresSigma_P4_2overmnm_Aminor_fields.in'
#path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt2Spec/networkO70_Fdddconv_Aminor_fields.in'
#path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt2Spec/networkO70_Fdddconv_Aminor_fields.in'
path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt/TwoSpecies/spheresBCCprim_Im-3m_Aminor_fields.in'

export_path = '/home/tquah/SEEDS/Polymer/BCCPhase/fields.in'
field_adjust = Nref/Nref_desire
hcell_adjust = np.sqrt(Nref/Nref_desire)


field = PolyFTSFieldReader()
field.readFields(path,True)
field_Re = field.AllFields/field_adjust
field_Im = field.AllFieldsImPart/field_adjust
hcell = field.hcell/hcell_adjust
NxNyNz = [16,16,16]
writePolyFTSBinFile(export_path,NxNyNz, hcell, True, True, field_Re, field_Im)

fieldtest = PolyFTSFieldReader()
fieldtest.readFields(export_path,True)


#
