#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  6 17:42:02 2021

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from scipy.interpolate import CubicSpline

def Reorder(listofarray,crit_array):
    sorting = np.argsort(listofarray[crit_array])
    newarryalist = []
    for array in listofarray:
        newarray = array[sorting]
        newarray = newarray[np.logical_not(np.isnan(newarray))]

        newarryalist.append(newarray)
    return newarryalist
def RemoveNan(listofarray):
    np.where
    


path_to_data = '/home/tquah/Projects/JournalData/zhulina.csv'



df = pd.read_csv(path_to_data)

#check if plot is correct

plt.close('all')
plt.figure()
header = list(df)
xfull = np.linspace(0.05,0.95)
for i in range(0,int(len(header)/2),1):
    
    x = np.array(df[header[2*i]])
    y = df[header[2*i+1]]
    newxy = Reorder([1-x,y],0)
    cs = CubicSpline(newxy[0],newxy[1],extrapolate=True)
    
    plt.plot(newxy[0],newxy[1],'--k')
    plt.plot(xfull,cs(xfull),'-r',alpha = 0.5)
plt.ylim(0,6)
plt.xlabel('$f_A$')
plt.ylabel(r'$ \beta^{1/2} \frac{\eta_B}{\eta_A}$')
plt.tight_layout()


plt.figure()
header = list(df)
xfull = np.linspace(0.00,0.98,1000)
for i in range(0,int(len(header)/2),1):
    
    x = np.array(df[header[2*i]])
    y = df[header[2*i+1]]
    newxy = Reorder([1-x,y],0)
    cs = CubicSpline(newxy[0],newxy[1],extrapolate=True)
    if i==0:

        plt.plot(xfull,cs(xfull),'--r',label = 'CL-Scaling')
        plt.plot(xfull,cs(xfull)**2,'-k',label = 'BB-Scaling')
    else:
        plt.plot(xfull,cs(xfull),'--r')
        plt.plot(xfull,cs(xfull)**2,'-k')

plt.legend()
plt.ylim(1,3.5)
plt.xlabel('$f_A$')
plt.ylabel(r'$\epsilon \sim \left( \frac{N_{sc,B}}{N_{sc,A}} \right)^{1/2} $')
plt.tight_layout()

