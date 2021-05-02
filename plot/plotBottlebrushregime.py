#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 25 14:01:17 2021

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt

Nsc = np.linspace(5,50,100)
zstar = 1/np.sqrt(Nsc)
z1 = 1/Nsc
plt.close('all')

plt.figure()
plt.plot(zstar,Nsc,'--k',label = '$z_1$ - Loosely grafted combs')
plt.plot(z1,Nsc,'-.k',label = '$z^{\star}$ - Densely grafted combs')

plt.xlabel('$z$')
plt.ylabel('$N_{sc}$')
plt.xlim(0.05,0.5)
plt.ylim(5,50)
plt.legend()

plt.figure()
cn = np.linspace(0,2,100)
zstar = cn*20**(-1/2)

gzero = np.where(cn<=1)[0]



zstarstar = np.ones_like(cn)*cn
zstarstar[gzero] = cn[gzero]*cn[gzero]
#zstarstar2 = cn
zstarstarstar=cn




plt.plot(cn,zstar,'.r',label = '$z^{\star}$ - Loosely grafted combs')
plt.plot(cn,zstarstar,'-.g',label = '$z^{\star \star}$ - Densely grafted combs')
plt.plot(cn,zstarstarstar,'-.b',label = '$z^{\star \star \star}$ - Brush')



plt.xlabel('$c/n$')
plt.ylabel('$z$')
plt.legend()
plt.savefig('/home/tquah/Figures/zvsc.pdf',dpi = 300)


plt.figure()


Nbb = 40
Nsc = 20
b = 1
l= 1
v =1 
 
z = np.linspace(0,2,100)

c = (l*z*Nsc)**(3/2)/v**(1/2)/(1+z*Nsc)/(1+z*Nsc)*Nbb
plt.plot(z,c)

c = z**(5/2)*Nsc**(3/2)/(1+z*Nsc)/(1+z*Nsc)*Nbb
plt.plot(z,c)

c = v*z**(3)*Nsc**(3/2)/(1+z*Nsc)/(1+z*Nsc)/(b*l)**3/2*Nbb
plt.plot(z,c)


