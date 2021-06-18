#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 17 15:44:43 2021

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
import pandas as pd
os.chdir('/media/tquah/TOSHIBA EXT/Projects/eps_development/BCC-DIS-bB_SC_1-1.5')

filelist = glob.glob('./**/F0_phases.dat',recursive=True)

datalist = []
for file in filelist:
    print(file)
    b = float(file.split('/')[1][2:])
    df = pd.read_csv(file,delimiter=' ',header=None)
    array = np.array([b,df[1][0],df[1][1]])
    datalist.append(array)
    
mainarray = np.vstack(datalist)
index = np.argsort(mainarray[:,0])
mainarray = mainarray[index,:]
mainarray[:,1] = mainarray[:,1]-mainarray[:,2]
mainarray[:,2] = mainarray[:,2]-mainarray[:,2]



plt.figure()
plt.plot(mainarray[:,0],mainarray[:,1],'-ok')
plt.plot(mainarray[:,0],mainarray[:,2],'-r^')
plt.ylabel(r'$\beta F-\beta F_{DIS}$')
plt.xlabel(r'$b_{B,SC}$')

plt.tight_layout()