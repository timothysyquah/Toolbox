#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul  4 14:41:42 2021

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'

rcParams['axes.labelsize'] = 23
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 15


def Cref_to_v0(Cref):
    return 1/(Cref*6**(3/2))

def alpha_ret(b,Nsc,z,v0):
    b3v0=(b**3)/v0
    b3v0Nsc_12=b3v0/np.sqrt(Nsc)
    if b3v0>1:
        if z> b3v0Nsc_12 and z<b3v0:
            alpha = (z*np.sqrt(Nsc))**3/b3v0/(1+z*Nsc)**2
        elif z> b3v0:
            alpha = np.sqrt(b3v0)*np.sqrt(Nsc*z)**3/(1+z*Nsc)**2
        else:
            alpha = 0
            print('Warning: Comblike Scaling')
    else:
        if z> b3v0Nsc_12 and z<(b3v0)**2:
            alpha = (z*np.sqrt(Nsc))**3/b3v0/(1+z*Nsc)**2
        
        elif z>(b3v0)**2 and z<np.sqrt(b3v0):
            alpha = z**(5/2)*Nsc**(3/2)/(1+z*Nsc)**2
        
        elif z> np.sqrt(b3v0):
            alpha = np.sqrt(b3v0)*np.sqrt(Nsc*z)**3/(1+z*Nsc)**2
    return alpha

def Nbar(Nbb,Nsc,z,cref,b):
    v0 = Cref_to_v0(cref)
    alpha = alpha_ret(b,Nsc,z,v0)
    return alpha*Nbb

def chiN_ODT(Nbar_array):
    return 10.495+41*np.power(Nbar_array,-1/3)+123*np.power(Nbar_array,-0.56)


Nbb = np.linspace(30,100,1000)#np.array([40,60,80,100])
b = 1
Nsc = 40
z = 0.5
chi = 0.5
Ntot = Nbb*(z*Nsc+1)
chiNtot = Ntot*chi
Nbar_array = Nbar(Nbb,Nsc,z,1/(6**(3/2)),b)
chiNODT = chiN_ODT(Nbar_array)
plt.close('all')
plt.figure()
plt.plot(Nbb,chiNtot,'r',label = rf'$\chi = {chi}$')
plt.plot(Nbb,chiNtot/2,'b',label = rf'$\chi = {chi/2}$')
plt.plot(Nbb,chiNODT,'k',label = rf'ODT')
plt.ylabel(r'$\chi N$')
plt.xlabel(r'$N_{bb}$')
plt.legend()
plt.tight_layout()
# plt.savefig('/home/tquah/Figures/ODT.pdf',dpi = 300)

Nbb = np.array([40,60,80,100])
b = 1
Nsc = 10
z = 2
chi = 0.2
Ntot = Nbb*(z*Nsc+1)
chiNtot = Ntot*chi
Nbar_array = Nbar(Nbb,Nsc,z,1/(6**(3/2)),b)
chiNODT = chiN_ODT(Nbar_array)
C = np.sqrt(Nbar_array)/6**3
Clin = np.sqrt(Ntot)*1/(6*np.sqrt(6))
