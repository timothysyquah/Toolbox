#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 12:00:09 2020

@author: tquah
"""

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import linregress
import matplotlib.pyplot as plt
import glob
from natsort import natsorted, ns
import os 
def powerlawscale(x,a,b):
    return a*x**b
import_path = '/home/tquah/IMPORT_BRAID/2020-03-05/LD_data.dat/Bottlebrushes/DGC**/*.dat'
export_path = '/home/tquah/Presentations/2020-03-06/'
name = 'linear_DvsL.pdf'
data = []

armlength = []

for file in glob.iglob(import_path,recursive=True):
    data.append(np.loadtxt(file))
    groups = file.split('/')
    specs = groups[-1].split('_')
    armlength.append(specs[-4]+' '+specs[-3])
sorted_armlength = natsorted(armlength, alg=ns.IGNORECASE)

colors = ['r','g','b','k','m','c']

plt.figure(0)
for j in range(0,len(data),1):
    i = armlength.index(sorted_armlength[j])
    popt,pcov = curve_fit(powerlawscale,data[i][-3:,0],data[i][-3:,1])
    print(popt)
    
    m,b,r,p,se = linregress(np.log10(data[i][-2:,0]),np.log10(data[i][-2:,1]))
    xx = np.log10(np.linspace(10,100,100))
    # if i!=len(armlength)-1:
    plt.plot(xx,m*xx+b,c=colors[j])
    plt.scatter(np.log10(data[i][:,0]),np.log10(data[i][:,1]),c=colors[j],label = armlength[i]+'=%0.3f'%m)
    
    print(armlength[i])
    
    print(m)
    print(b)
    print(r**2)
    
#    print(popt)
plt.legend()
plt.xlabel('$\log_{10}(N)$')
plt.ylabel('$\log_{10}(D_0\sqrt{6}/b_{ref})$')
epath =os.path.join(export_path,name)
plt.savefig(epath,dpi = 300)