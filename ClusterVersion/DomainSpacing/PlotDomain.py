#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 11:21:17 2021

@author: tquah
"""


import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
plt.close('all')
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
# rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'
rcParams['font.style'] = 'normal'
rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 15
rcParams['axes.titlesize'] = 20

os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/rho_1_z_2_point_invest')

cl = np.loadtxt('CL_domain_data.dat')
scft = np.loadtxt('SCFT_domain_data.dat')

def powlaw(x,a,b):
    return np.power(x,b)*a

def powlaw2(x,a):
    b = 2/3
    return np.power(x,b)*a



bounds = np.array([-3,None])
plt.loglog(cl[:,2],cl[:,3],marker = '^',color = 'r',linewidth=0,label = 'FTS-CL')
var1,pcov1 = curve_fit(powlaw,cl[bounds[0]:bounds[1],2],cl[bounds[0]:bounds[1],3])


plt.loglog(scft[:,2],scft[:,3],marker = 's',color = 'k',linewidth=0,label = 'SCFT')
var2,pcov2 = curve_fit(powlaw,scft[bounds[0]:bounds[1],2],scft[bounds[0]:bounds[1],3])




x = np.linspace(np.min(scft[:,2]),np.max(scft[:,2]),100)
y = powlaw(x,var2[0],var2[1])
plt.loglog(x,y,'--k')

x = np.linspace(np.min(cl[:,2]),np.max(scft[:,2]),100)
y = powlaw(x,var1[0],var1[1])
plt.loglog(x,y,'-.r')
bounds = np.array([-2,None])
var1,pcov1 = curve_fit(powlaw,cl[bounds[0]:bounds[1],2],cl[bounds[0]:bounds[1],3])
y = powlaw(x,var1[0],var1[1])
plt.loglog(x,y,':r')

var1,pcov1 = curve_fit(powlaw2,cl[bounds[0]:bounds[1],2],cl[bounds[0]:bounds[1],3])
y = powlaw2(x,5)
plt.loglog(x,y,'-k')

plt.xticks([100,140,200])
plt.yticks([80,100, 140,170])
plt.legend()
plt.xlabel('$N_{bb}$')
plt.ylabel('$D(l)$')
plt.tight_layout()
plt.savefig('/home/tquah/Figures/Domain.svg',dpi = 300)