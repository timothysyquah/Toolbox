#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Dec 25 22:45:17 2020

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
import glob
import os
import pandas as pd
import regex as re
from resubmit_ODTlike_calculations import *

os.chdir('/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_asymBCC_fix/')
wdir = os.getcwd()
dirs = glob.glob("chiAB_*/Nsc*/fA*")



def unique_category_values(directory,file='diagnose.dat'):
    data_path = os.path.join(directory,file)
    df = pd.read_csv(data_path,delimiter=' ',header=None)
    phases = list(df[0])
    phaselist = []
    for i in range(0,len(phases)):
        temp_tuple = (phases[i],float(df[1][i]),int(df[2][i]),float(np.round(df[3][i],3)))
        phaselist.append(temp_tuple)
        
    category = directory.split('/')
    category_values = []
    for i in range(0,len(category)):
        extract_values = re.findall('[\d]*[.][\d]+' ,category[i])
        if len(extract_values)==1:
            category_values.append(float(extract_values[0]))
        elif len(extract_values)>1:
            t = ()
            for tup in extract_values:
                t = t + (float(tup),)
            category_values.append(t)
    category_values.append(phaselist)
    return category_values


