#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar 24 12:15:25 2020

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt

IMPORT_PATH = '/home/tquah/IMPORT_BRAID/L_159/LAM.out'

op = open(IMPORT_PATH,'r')
splitfile = op.read().split('\n')
op.close()
Q = []
for i in range(0,len(splitfile),1):
    r = list(filter(lambda x: 'Partition Function, Q' in x, splitfile[i]))[0]
    Q.append(r)
# Box_initial = float(r[r.find("(")+1:r.find(")")])
