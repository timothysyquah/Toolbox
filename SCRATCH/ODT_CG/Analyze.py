#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 25 16:49:33 2021

@author: tquah
"""

import numpy as np
import glob
import itertools
import pandas as pd
import matplotlib.pyplot as plt
plt.close('all')

#note some datasets are flucutating over a nonphysical 

filelist = glob.glob('*.dat')
marker = itertools.cycle(('s', '^', 'P', 'o', '*')) 
fig = plt.figure(figsize = (12,6))
ax = plt.gca()
for i in range(len(filelist)):
    
    newmarker = next(marker)
    array = np.loadtxt(filelist[i])
    header = list(pd.read_csv(filelist[i],sep = " "))[1:]
    
    color = itertools.cycle(('lightcoral', 'peru', 'darkorange', 'gold','lawngreen','aquamarine','teal','skyblue','dodgerblue','rebeccapurple','mediumorchid')) 

    Cunique = np.unique(array[:,0])
    for j in range(len(Cunique)):
        newcolor = next(color)

        loc = np.where(array[:,0]==Cunique[j])[0]
        reorder1 = np.argsort(array[loc,1])
        reorder2 = np.argsort(array[loc,4])
        reorder3 = np.argsort(array[loc,6])

        ax.errorbar(array[loc,1][reorder1]*100,array[loc,2][reorder1],label = f'C = {Cunique[j]*10}',marker = newmarker,color = newcolor,linestyle = '--',yerr=array[loc,3])
        ax.errorbar(array[loc,4][reorder2]*100,array[loc,2][reorder2],label = f'C = {Cunique[j]*10}',marker = newmarker,color = newcolor,linestyle = '--',yerr=array[loc,3])
        ax.errorbar(array[loc,6][reorder3]*100,array[loc,2][reorder3],label = f'C = {Cunique[j]*10}',marker = newmarker,color = newcolor,linestyle = '--',yerr=array[loc,3])
        


ax.legend(loc=(1.01, 0.01), ncol=4)
plt.tight_layout()
