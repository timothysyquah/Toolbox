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
from itertools import combinations
phases = ['DIS.dat','SIGMA.dat','BCC.dat','HEX.dat','A15.dat']
phases = ['SIGMA.dat','A15.dat','HEX.dat']

pts = [[0.26,2.7],[0.25,2.2],[0.26,2.2],[0.24,1.9],[0.25,1.9],[0.26,1.9]]
pairlist = [[0,1],[0,2],[1,3],[1,4],[1,5],[2,5]]
tol = 50
plt.close('all')
def return_likely_point(pt,data):
    r = np.sum((data[:,0:2]-pt)**2,axis=1)
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






#norm to A15
normalize = deepcopy(realpts_energy[1][:,-1])
realpts_energy[0][5,2] = 3.9413196247e-03
realpts_energy[0][4,2] = 3.9088089519e-03

# for i in range(len(realpts_energy)):
#     realpts_energy[i][:,-1] = (realpts_energy[i][:,-1]-normalize)

for i in range(len(pairlist)):
    plt.figure()
    title = f'Pair{pairlist[i][0]}:{pairlist[i][1]}'
    plt.title(title)
    linregress_var = []
    for j in range(len(realpts_energy)):
        
        # x1 = np.array([realpts_energy[j][pairlist[i][0]][0],realpts_energy[j][pairlist[i][1]][0]])
        # x2 = np.array([realpts_energy[j][pairlist[i][0]][1],realpts_energy[j][pairlist[i][1]][1]])
        # x3 = np.sqrt((x1[1]-x1[0])**2+(x2[1]-x2[0])**2)
        # x = np.array([0,x3])
        x = np.array([realpts_energy[j][pairlist[i][0]][1],realpts_energy[j][pairlist[i][1]][1]])
        y = np.array([realpts_energy[j][pairlist[i][0]][-1],realpts_energy[j][pairlist[i][1]][-1]])
        error = np.sum(y**2)
        m,b,r,p,std = linregress(x,y)
        linregress_var.append(np.array([m,b,np.min(x),np.max(x)]))
        if error<tol:
            plt.plot(x,y,label = phases[j])
    
    
    combos = list(combinations(np.arange(0,len(realpts_energy)-1+1e-6,1,dtype = int),2))
    
    
    for j in range(len(combos)):
        
        intersection = get_intersection(linregress_var[combos[j][0]][0],\
                                        linregress_var[combos[j][1]][0],\
                                        linregress_var[combos[j][0]][1],\
                                        linregress_var[combos[j][1]][1])
        yint = linregress_var[combos[j][0]][0]*intersection+\
            linregress_var[combos[j][0]][1]
        print(intersection)
        if intersection>linregress_var[combos[j][0]][2] and intersection<linregress_var[combos[j][0]][3]:
            plt.plot(intersection,yint,'ok')
    
    plt.savefig(title+'.svg',dpi = 300)
    plt.legend()