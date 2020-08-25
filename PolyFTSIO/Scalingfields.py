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
path = '/home/tquah/PolyFTS_ALL/PolyFTS/seeds/BlockPolymerMelt2Spec/spheresSigma_P4_2overmnm_Aminor_fields.in'
export_path = '/home/tquah/SEEDS/SigmaPhase/nref1.0/NscA12.0_NscB28.0/fA0.33000/SigmaPhase/fields.in'
field_adjust = Nref/Nref_desire
hcell_adjust = np.sqrt(Nref/Nref_desire)


field = PolyFTSFieldReader()
field.readFields(path,True)
field_Re = field.AllFields/field_adjust
field_Im = field.AllFieldsImPart/field_adjust
hcell = field.hcell/hcell_adjust
NxNyNz = [64,64,32]
writePolyFTSBinFile(export_path,NxNyNz, hcell, True, True, field_Re, field_Im)

fieldtest = PolyFTSFieldReader()
fieldtest.readFields(export_path,True)


#
