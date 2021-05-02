#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 11 13:15:56 2021

@author: tquah
"""

#maybe see if we can display in 3D
import pandas as pd
import numpy as np
from mpl_toolkits import mplot3d
import matplotlib.pyplot as plt

path = '/media/tquah/TOSHIBA EXT/Projects/PhaseDiagram/Experimental_Dataset_v2.csv'
df = pd.read_csv(path,skiprows=np.arange(1,37.5,1,dtype=int))
#path to symmetrical bottlebrush

#path to asymmtrical bottlebrush


#path to asymmetrical bottlebrush overall







chiN = np.array([97.3,78,121.9,109,105.9, 102.6, 99.3, 97.9, 94.9])
alpha2 = (np.array(df['NscB'])+1)/np.array(df['NscA'])+1
alpha = np.sqrt(alpha2)
fA = np.array(df['fA'])
phases = list(df['Phase'])
uniquephases = list(set(phases))
colors = ['k','r','g','b']
fig = plt.figure()
ax = plt.axes(projection='3d')

for i in range(len(uniquephases)):
    indices = [j for j, x in enumerate(phases) if x == uniquephases[i]]
    print(indices)
    phases_itter = [phases[j] for j in indices]
    print(phases_itter)
    
    ax.scatter3D(fA[indices], chiN[indices],alpha[indices],c= colors[i])
