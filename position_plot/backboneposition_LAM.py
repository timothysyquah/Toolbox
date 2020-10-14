#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 12:51:38 2020

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'

rcParams['axes.labelsize'] =20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 11


color = ['tab:blue','tab:orange','tab:red','tab:olive','tab:brown']

path = '/home/tquah/Projects/positions/Test_1/density_chain0.dat'

data = np.loadtxt(path)
# desired_list = [1,3,5,6,8,10]
desired_list = [6,8,10]

# alpha = [0.5,0.75,1.0,1.0,0.75,0.5]
# color = ['r','r','r','b','b','b']
color = ['r','r','r']
alpha = [1.0,0.75,0.5]


shape = np.shape(data)
plt.close('all')
plt.figure()
count = 0
for i in range(1,shape[1]):
    if i<shape[1]-2:
        
        if i in desired_list:
            
            cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
            integral = cs.integrate(0,1)
            x = np.linspace(0,1,1000)
            y = cs(x)/integral
        
            plt.plot(x,y,label = str(i),alpha = alpha[count],color = color[count])
            count +=1
# plt.legend()
color = ['b']*3

path = '/home/tquah/Projects/positions/Test_linear_1/density_chain0.dat'

data = np.loadtxt(path)

shape = np.shape(data)
# plt.close('all')
# plt.figure()
count=0
for i in range(1,shape[1]):
    if i<shape[1]:
        if i in desired_list:

            cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
            integral = cs.integrate(0,1)
            x = np.linspace(0,1,1000)
            y = cs(x)/integral
            plt.plot(x,y,'--',alpha = alpha[count],color = color[count])
            count +=1





















# plt.legend()
plt.xlabel('$x^* = x/L$')
plt.ylabel(r'$\rho(x^*) / \int dx^* \rho(x^*)$')
























plt.tight_layout()
path = '/home/tquah/Projects/positions/Test_1/density.dat'
data = np.loadtxt(path)

shape = np.shape(data)
plt.figure()
for i in range(1,shape[1]):
        
    cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
    
    integral = cs.integrate(0,1)
    x = np.linspace(0,1,1000)
    y = cs(x)/integral
    plt.plot(x,y,label = 'bottlebrush')

path = '/home/tquah/Projects/positions/Test_linear_1/density.dat'
data = np.loadtxt(path)

shape = np.shape(data)
for i in range(1,shape[1]):
        
    cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
    
    integral = cs.integrate(0,1)
    x = np.linspace(0,1,1000)
    y = cs(x)/integral
    plt.plot(x,y,label ='linear')
plt.legend()