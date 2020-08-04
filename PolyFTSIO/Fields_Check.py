#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 13:46:23 2020

@author: tquah
"""


from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np
import parser


#paths to fields 1 and 2
Fields_1 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields.in'
Fields_2 = '/home/tquah/IMPORT_BRAID/fieldstest/GYR/fields_k.bin'





fields = PolyFTSFieldReader()
try:
    fields.readFields(Fields_1,True)
except:
    print "*** Error during field reading ***"

f1 = fields.AllFields

fields2 = PolyFTSFieldReader()
try:
    fields2.readFields(Fields_2,True)
except:
    print "*** Error during field reading ***"

f2 = fields2.AllFields

f1norm1 = f1[0,:]/np.linalg.norm(f1[0,:])
f2norm1 = f2[0,:]/np.linalg.norm(f2[0,:]) 
check = np.dot(f1norm1,f2norm1)
print(check)