#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 13:19:06 2021

@author: tquah
"""
import itertools
import numpy as np
import matplotlib.pyplot as plt
# plt.close('all')
array = np.loadtxt('resultsQ_linear100.dat')

array = np.loadtxt('resultsQ_point.dat')
aunique = np.unique(array[:,1])
zetaunique = np.unique(array[:,2])

for i in range(0,len(aunique)):
    marker = itertools.cycle(('s', 'p', '^', 'o', '*')) 
    color = itertools.cycle(('r', 'g', 'b', 'k')) 

    plt.figure()
    
    for j in range(0,len(zetaunique)):
        
        tempmarker = next(marker)
        tempcolor = next(color)
        loca = np.where(aunique[i]==array[:,1])[0]
        locz = np.where(zetaunique[j]==array[:,2])[0]
        rows = np.intersect1d(loca,locz)
        plt.errorbar(array[rows,0],array[rows,3],yerr = array[rows,4],marker = tempmarker,color =tempcolor,linestyle = '',label = f'$\zeta = {zetaunique[j]}$, $a = {aunique[i]}$')
    plt.xscale('log')
    # plt.yscale('log')
    plt.ylim(-0.1,1.1)
    plt.xlabel('$C_{ref}$')
    plt.ylabel(r'$Q[w]$')
    plt.legend()
    plt.tight_layout()
    plt.savefig(f'/home/tquah/Figures/a{aunique[i]}_Q.pdf',dpi = 300)
    