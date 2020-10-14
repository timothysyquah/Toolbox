#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 15 12:51:38 2020

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline


color = ['tab:blue','tab:orange','tab:red','tab:olive','tab:brown']

path = '/home/tquah/Projects/positions/Test/density_chain0.dat'

data = np.loadtxt(path)

shape = np.shape(data)
plt.close('all')
plt.figure()
for i in range(1,shape[1]):
    if i<shape[1]-2:
        
        cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
        integral = cs.integrate(0,1)
        x = np.linspace(0,1,1000)
        y = cs(x)/integral
        
        plt.plot(x,y,label = str(i))
# plt.legend()

path = '/home/tquah/Projects/positions/Test_linear/density_chain0.dat'

data = np.loadtxt(path)

shape = np.shape(data)
# plt.close('all')
# plt.figure()
for i in range(1,shape[1]):
    if i<shape[1]:
        
        cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
        integral = cs.integrate(0,1)
        x = np.linspace(0,1,1000)
        y = cs(x)/integral
        plt.plot(x,y,'--')
plt.legend()
plt.xlabel('$x^* = x/L$')
plt.ylabel(r'$\rho(x^*) / \int dx \rho(x^*)$')


plt.tight_layout()
path = '/home/tquah/Projects/positions/Test/density.dat'
data = np.loadtxt(path)

shape = np.shape(data)
plt.figure()
for i in range(1,shape[1]):
        
    cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
    
    integral = cs.integrate(0,1)
    x = np.linspace(0,1,1000)
    y = cs(x)/integral
    plt.plot(x,y,label = 'bottlebrush')

path = '/home/tquah/Projects/positions/Test_linear/density.dat'
data = np.loadtxt(path)

shape = np.shape(data)
for i in range(1,shape[1]):
        
    cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
    
    integral = cs.integrate(0,1)
    x = np.linspace(0,1,1000)
    y = cs(x)/integral
    plt.plot(x,y,label ='linear')
plt.legend()