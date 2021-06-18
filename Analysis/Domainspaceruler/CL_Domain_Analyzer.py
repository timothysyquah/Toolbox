#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 11 11:26:33 2021

@author: tquah
"""

import numpy as np 
import matplotlib.pyplot as plt
import glob
import matplotlib
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline
import re
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'

rcParams['axes.labelsize'] = 23
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 20
from matplotlib.lines import Line2D
import os
def powerlawscale(x,a,b):
    return a*x**b


pathimport = '/home/tquah/Projects/CL_analysis/C_2'
exportpath='/home/tquah/Presentations/Candidacy/Figures'
exportpath='/home/tquah/Figures'

exportname = 'CLSCFT_new.pdf'
exportnamegamma = 'gamma_new.pdf'


fullpathname = os.path.join(exportpath,exportname)
fullpathnamegamma = os.path.join(exportpath,exportnamegamma)

filelist = glob.glob(pathimport+'/nsc*.dat')





plt.close('all')
fig, ax = plt.subplots(figsize=(5 , 5))
ax = plt.subplot(111)
plt.plot([], [])


marker = ['o','^']

def filelistorder(filelist):
    place = []
    for i in range(0,len(filelist),1):
        splitvalue = filelist[i].split('/')
        number = re.findall(r'\d+', splitvalue[-1])[0]
        print(number)
        place.append([float(number),i])
    array = np.vstack(place)
    print(array)
    return np.argsort(array[:,0])
    
print(filelist)
order = filelistorder(filelist)
print(order)
arraylist = []
for i in range(0,len(filelist)):
    dataarray = np.loadtxt(filelist[order[i]])
    popt,pcov = curve_fit(powerlawscale,dataarray[-2:,2],dataarray[-2:,3])
    xx = np.linspace(np.min(dataarray[:,2]),np.max(dataarray[:,2]),1000)
    ax.loglog(xx,powerlawscale(xx,popt[0],popt[1]),c = 'r')

    ax.scatter(dataarray[:,2],dataarray[:,3],marker = marker[i],color = 'r')
    arraylist.append(dataarray[:,2:])
    
filelist = glob.glob(pathimport+'/SCFT*.dat')
print(filelist)
order = filelistorder(filelist)
print(order)

for i in range(0,len(filelist)):
    dataarray = np.loadtxt(filelist[order[i]])
    sort = np.argsort(dataarray[:,0])
    dataarray = dataarray[sort,:]
    
    
    popt,pcov = curve_fit(powerlawscale,dataarray[-2:,0],dataarray[-2:,1])
    xx = np.linspace(np.min(dataarray[:,0]),np.max(dataarray[:,0]),1000)
    ax.loglog(xx,powerlawscale(xx,popt[0],popt[1]),c = 'k')

    ax.scatter(dataarray[:,0],dataarray[:,1],marker = marker[i],color = 'k')
    arraylist.append(dataarray)

    
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xticks([40, 60, 100, 200])
ax.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.set_yticks([60,100, 150, 200])
ax.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())

plt.xlabel('$N_{bb}$')
plt.ylabel('$D/l$')

legendhandles = [Line2D([0], [0], marker = 'o',markerfacecolor = 'k',label = r'$N_{sc} = 20$',color = 'w',linestyle='None',markersize=10,alpha = 0.5),\
                 Line2D([0], [0], marker = '^',markerfacecolor = 'k',label = r'$N_{sc} = 40$',color = 'w',linestyle='None',markersize=10,alpha = 0.5),\
                 Line2D([0], [0], color = 'r',label = r'FTS-CL',linewidth = 3.0),\
                 Line2D([0], [0], color = 'k',label = r'SCFT',linewidth = 3.0)]



plt.tight_layout()
ax.legend(handles =legendhandles,loc='upper left',frameon=False)
plt.savefig(fullpathname,dpi = 300)
#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################
fig, ax = plt.subplots(figsize=(5 , 5))
ax = plt.subplot(111)
plt.plot([], [])
color = ['-r','-.r','-k','-.k']
for i in range(len(arraylist)):
    newarray = arraylist[i]
    arrayfit = (newarray[0:-1,0]+newarray[1:,1])/2

    xx = np.linspace(np.min(arrayfit),np.max(arrayfit),1000)
    cs = CubicSpline(np.log(arraylist[i][:,0]),np.log(arraylist[i][:,1]))
    Dspace_dX = cs.derivative(1)
    ax.plot(xx,Dspace_dX(np.log(xx)),color[i])
    
legendhandles = [Line2D([0], [0], color = 'k',linestyle = '-',label = r'$N_{sc}=20$',linewidth = 3.0,alpha = 0.5),\
             Line2D([0], [0], color = 'k',linestyle = '-.',label = r'$N_{sc}=40$',linewidth = 3.0,alpha = 0.5),\
             Line2D([0], [0], color = 'r',label = r'FTS-CL',linewidth = 3.0),\
             Line2D([0], [0], color = 'k',label = r'SCFT',linewidth = 3.0)]

plt.plot(np.linspace(-100,1000,1000),2/3*np.ones(1000),'--k')
ax.set_xlim(50,200)
    
ax.legend(handles =legendhandles,loc='lower right',frameon=False)

plt.xlabel('$N_{bb}$')
plt.ylabel('$\gamma$')
plt.tight_layout()
plt.savefig(fullpathnamegamma,dpi = 300)