#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 14:47:27 2020

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
plt.close('all')
from matplotlib.ticker import AutoMinorLocator

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
from matplotlib.ticker import MultipleLocator

import os 

def get_density_profiles(mainpath,desired_list):
    #determine the location of the period
    datalist = []
    
    path_maindensity = os.path.join(mainpath,'density.dat')
    data = np.loadtxt(path_maindensity)
    fAmin = np.min(data[:,1])
    cs = CubicSpline(data[:,0],data[:,1]-fAmin)
    x = np.linspace(np.min(data[:,0]),np.max(data[:,0]),1000)
    roots = cs.roots()
    loc = np.where((roots>np.min(data[:,0])) & (roots<np.max(data[:,0])))[0]
    domain_min = np.min(roots[loc[1]])
    domain_max = np.max(roots[loc[3]])
    loc = np.where((data[:,0]>domain_min)&(data[:,0]<domain_max))[0]   
    data_domain = data[loc,:]
    data_domain[:,0] = data_domain[:,0]-np.min(data_domain[:,0])
    cs = CubicSpline(data_domain[:,0]/np.max(data_domain[:,0]),data_domain[:,1])
    integral = 1#cs.integrate(0,1)
    x = np.linspace(0,1,10000)
    y = cs(x)/integral
    datalist.append(np.vstack((x,y)).transpose())
    path_data = os.path.join(mainpath,'density_chain0.dat')

    all_data = np.loadtxt(path_data)

    loc = np.where((all_data[:,0]>domain_min)&(all_data[:,0]<domain_max))[0]
    data = all_data[loc,:]
    data[:,0] = data[:,0]-np.min(data[:,0])
    shape = np.shape(data)
    for i in range(1,shape[1]):
        
        if i in desired_list:
                
            cs = CubicSpline(data[:,0]/np.max(data[:,0]),data[:,i])
            integral = cs.integrate(0,1)
            x = np.linspace(0,1,10000)
            y = cs(x)/integral
            datalist.append(np.vstack((x,y)).transpose())
    return datalist            
            
os.chdir('/media/tquah/TOSHIBA EXT/Projects/positions/sidearmposition/DGC/lowchicomp')

color = ['b','r','g']
linestyle = ['dotted','dashed','dashdot']
alpha = [1.0,0.66,0.33]
path = 'Bottlebrush_Sidechains/density.dat'
desired_list = [[3,7,11],[1,3,5],[1,3,5]]
listofdirectories = ['Bottlebrush_Sidechains','Bottlebrush_Backbones','Linear_Backbones']
fig, ax = plt.subplots()
general_density = []
for i in range(len(listofdirectories)):    
    
    datalist = get_density_profiles(listofdirectories[i],desired_list[i])
    for j in range(len(datalist)):
        if j==0:
            general_density.append(datalist[j])
        else:
            ax.plot(datalist[j][:,0],datalist[j][:,1],color = color[i],alpha = alpha[j-1],linestyle=linestyle[i],linewidth=3.0)

ml = MultipleLocator(10)
minor_locator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator)

plt.xticks(np.arange(0.1,1.0+1e-6,0.2))
plt.xlabel(r'$\bar{x}$')
plt.ylabel(r'$\rho(\bar{x}) / \int d\bar{x} \rho(\bar{x})$')
plt.tight_layout()
################################################################################################################################################


fig, ax = plt.subplots()
ax.plot(general_density[1][:,0],general_density[1][:,1],color = 'b',linestyle='dashed',linewidth=3.0)
# ax.plot(1-general_density[1][:,0],general_density[1][:,1],color = 'r',linestyle='dashed',linewidth=3.0)

ax.plot(general_density[2][:,0],general_density[2][:,1],color = 'g',linestyle='dashdot',linewidth=3.0)
# ax.plot(1-general_density[2][:,0],general_density[2][:,1],color = 'r',linestyle='dashed',linewidth=3.0)

ml = MultipleLocator(10)
minor_locator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator)

plt.xticks(np.arange(0.1,1.0+1e-6,0.2))

plt.xlabel(r'$\bar{x}$')
plt.ylabel(r'$\rho(\bar{x})$')
plt.tight_layout()







###############################################################################################################################################
os.chdir('/media/tquah/TOSHIBA EXT/Projects/positions/sidearmposition/DGC/highchicomp')

color = ['b','r','g']
linestyle = ['dotted','dashed','dashdot']
alpha = [1.0,0.66,0.33]
path = 'Bottlebrush_Sidechains/density.dat'
desired_list = [[3,7,11],[1,3,5],[1,3,5]]
listofdirectories = ['Bottlebrush_Sidechains','Bottlebrush_Backbones','Linear_Backbones']
fig, ax = plt.subplots()
general_density = []

for i in range(len(listofdirectories)):    
    
    datalist = get_density_profiles(listofdirectories[i],desired_list[i])
    for j in range(len(datalist)):
        if j==0:
            general_density.append(datalist[j])

        else:
            ax.plot(datalist[j][:,0],datalist[j][:,1],color = color[i],alpha = alpha[j-1],linestyle=linestyle[i],linewidth=3.0)

ml = MultipleLocator(10)
minor_locator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator)

plt.xticks(np.arange(0.1,1.0+1e-6,0.2))
plt.xlabel(r'$\bar{x}$')
plt.ylabel(r'$\rho(\bar{x}) / \int d\bar{x} \rho(\bar{x})$')
plt.tight_layout()
#################################################################################################
fig, ax = plt.subplots()
ax.plot(general_density[0][:,0],general_density[0][:,1],color = 'b',linestyle='dashed',linewidth=3.0)
ax.plot(general_density[2][:,0],general_density[2][:,1],color = 'g',linestyle='dashdot',linewidth=3.0)

ml = MultipleLocator(10)
minor_locator = AutoMinorLocator(2)
ax.xaxis.set_minor_locator(minor_locator)

plt.xticks(np.arange(0.1,1.0+1e-6,0.2))

plt.xlabel(r'$\bar{x}$')
plt.ylabel(r'$\rho(\bar{x})$')
plt.tight_layout()


