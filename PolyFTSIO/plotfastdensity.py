#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  8 10:35:38 2020

@author: tquah
"""
from PolyFTSFieldWriter import *
from PolyFTSFieldReader import *
import numpy as np
import argparse
import os
import matplotlib.pyplot as plt
from Fields_Check import fields_compare
from matplotlib import cm



path = '/home/tquah/Projects/sweep-asym-armlength_corrected/chiAB_0.0289/NscA_2.0_NscB_38.0/fA0.53305/HEXPhase/density.bin'


field = PolyFTSFieldReader()
field.readFields(path,True)
coords = field.AllMeshCoords
density = field.AllFields[0,:]
colors = cm.ScalarMappable('jet').to_rgba(density,norm=False)



plt.scatter(coords[0,:],coords[1,:],density)