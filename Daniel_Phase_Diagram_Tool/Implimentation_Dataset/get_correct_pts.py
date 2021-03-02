#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 22 10:38:50 2021

@author: tquah
"""

import numpy as np
from scipy.stats import linregress
import matplotlib.pyplot as plt
import os
from copy import deepcopy
from itertools import permutations, repeat
phases = ['DIS.dat','SIGMA.dat','BCC.dat','HEX.dat','A15.dat']
phases = ['SIGMA.dat','A15.dat','BCC.dat']
fA = [0.1865,0.27,0.14]
Abratio = [1.0,1.68,1.9,2.7]


tuples = [list(zip(fA, p)) for p in permutations(Abratio)]

pts = []
count = 0
for lst in tuples:
    for tupe in lst:
       pts.append(tupe) 
pts = list(set(pts))
pts = [list(a) for a in pts]


def return_likely_point(pt,data):
    r = np.sum((data[:,0:2]-pt)**2,axis=1)
    print()
    loc = np.where(r==np.min(r))[0]
    return data[loc]



def get_intersection(m1,m2,b1,b2):
    if (m1-m2)==0:
        return np.nan
    else:
        return (b2-b1)/(m1-m2)







# data = []
IDIR = os.getcwd()

realpts_energy = []
i = 0
for phase in phases:
    os.chdir('PHASE_FREE_ENERGY_all')

    data = np.loadtxt(phase)
    realpts_temp = []
    for pt in pts:
        realpts_temp.append(return_likely_point(pt,data))
    realpts_energy.append(np.vstack(realpts_temp))
            
        
    os.chdir(IDIR)
    i+=1




