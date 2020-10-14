#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 14:47:27 2020

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

path = '/home/tquah/Projects/positions/lowchicomp/Test_1/density_chain0.dat'

data = np.loadtxt(path)
# desired_list = [1,3,5,6,8,10]
desired_list = [1,3,5]

# alpha = [0.5,0.75,1.0,1.0,0.75,0.5]
# color = ['r','r','r','b','b','b']
color = ['r','r','r']
alpha = [1.0,0.66,0.33]


shape = np.shape(data)
plt.close('all')
plt.figure()
count = 0
for i in range(1,shape[1]):
    if i in desired_list:
            
        cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
        integral = cs.integrate(0,1)
        x = np.linspace(0,1,1000)
        y = cs(x)/integral
    
        plt.plot(x,y,label = str(i),alpha = alpha[count],color = color[count])
        count +=1


path = '/home/tquah/Projects/positions/lowchicomp/Side_1/density_chain0.dat'
data = np.loadtxt(path)
shape = np.shape(data)
alpha = [0.33,0.66,1.0]

count = 0
desired_list = [3,7,11]
color = ['b']*len(desired_list)
# alpha = [1.]*len(desired_list)
for i in range(1,shape[1],1):
    if i in desired_list:
        
        cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
        integral = cs.integrate(0,1)
        x = np.linspace(0,1,1000)
        y = cs(x)/integral
    
        plt.plot(x,y,'--',alpha = alpha[count],color = color[count])
        count +=1
        
        
        
color = ['g']*3
alpha = [1.0,0.66,0.33]

path = '/home/tquah/Projects/positions/lowchicomp/Test_linear_1/density_chain0.dat'

data = np.loadtxt(path)

shape = np.shape(data)
# plt.close('all')
# plt.figure()
desired_list = [1,3,5]

count=0
for i in range(1,shape[1]):
    if i in desired_list:
        print(i)
        cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
        integral = cs.integrate(0,1)
        x = np.linspace(0,1,1000)
        y = cs(x)/integral
        plt.plot(x,y,'-.',alpha = alpha[count],color = color[count])
        count +=1

plt.xlabel('$x^* = x/L$')
plt.ylabel(r'$\rho(x^*) / \int dx^* \rho(x^*)$')
plt.tight_layout()
plt.savefig('/home/tquah/Presentations/makefigures/weakcompare_chaindensity.png')


path = '/home/tquah/Projects/positions/lowchicomp/Test_linear_1/density.dat'
data = np.loadtxt(path)

shape = np.shape(data)
plt.figure()

for i in range(1,shape[1]):
        
    cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
    
    integral = cs.integrate(0,1)
    x = np.linspace(0,1,1000)
    y = cs(x)/integral
    if i==1:
        plt.plot(x,y,'r',label = 'bottlebrush')
    else:
        plt.plot(x,y,'r')

path = '/home/tquah/Projects/positions/lowchicomp/Test_1/density.dat'
data = np.loadtxt(path)

shape = np.shape(data)
for i in range(1,shape[1]):
        
    cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
    
    integral = cs.integrate(0,1)
    x = np.linspace(0,1,1000)
    y = cs(x)/integral
    if i==1:
        plt.plot(x,y,'b',label = 'linear')
    else:
        plt.plot(x,y,'b')

plt.legend()
plt.xlabel('$x^* = x/L$')
plt.ylabel(r'$\rho(x^*) / \int dx^* \rho(x^*)$')
plt.tight_layout()
plt.savefig('/home/tquah/Presentations/makefigures/weakcompare_fulldensity.png')
