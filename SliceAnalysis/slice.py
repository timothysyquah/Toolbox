#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 12 15:04:10 2021

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

import os
import glob, os

plt.figure()
count = 0

phases= ['BCC','$\sigma$','A15','HEX']
label = ['$b_{bb,A}=1.24$','$b_{SC,B}=1.24$','$\lambda=2$','Base']
for file in glob.glob("*.dat"):
    print(file)
    data = np.loadtxt(file)
    for i in range(0,len(data[:,0])-1):
        plt.plot(data[i,0]*np.ones(10),count+np.linspace(0,1,10),'k') 
        xmean = (data[i,0]+data[i+1,0])/2
        plt.text(xmean,count+0.5,phases[i],ha = 'center')
    
    plt.text(0.00,count+0.5,s = label[count],c = 'r')
    count+=1
plt.xlim(0.00,0.35)
plt.yticks([])
plt.title(r'$\epsilon \approx 2.7$')
plt.xlabel(r'$f_A$')
plt.tight_layout()
plt.savefig('SliceEps.pdf',dpi = 300)