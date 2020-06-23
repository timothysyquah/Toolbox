#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 11:15:31 2020

@author: tquah
"""

import os
import numpy as np

path = '/home/tquah/IMPORT_BRAID/N_4_f_0.3_AGYR32_density.dat'
path_export = '/home/tquah/ParaviewInport/GyroidProblems/N4f0.3_AGYR32.csv'
data = np.loadtxt(path,skiprows=6)
np.savetxt(path_export,data,delimiter = ',')