#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 23 15:48:29 2021

@author: tquah
"""

import os
import numpy as np
import matplotlib.pyplot as plt

os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/SCFT_into_CL/CL/CL_0.01/REFINE')

fig = plt.figure()
ax = plt.gca()


data = np.loadtxt('L_98/SKAA.dat')


ax.scatter(data[:,0],data[:,1])
ax.plot(data[:,0],data[:,1],label ='L=98')

ax.set_yscale('symlog')
ax.set_xlim(0,1)
ax.set_xlabel('$q(l)$')
ax.set_ylabel('$S_{AA}(q)/CN$')




data = np.loadtxt('L_100/SKAA.dat')

ax.scatter(data[:,0],data[:,1])
ax.plot(data[:,0],data[:,1],label ='L=100')


data = np.loadtxt('L_102/SKAA.dat')


ax.scatter(data[:,0],data[:,1])
ax.plot(data[:,0],data[:,1],label ='L=102')

ax.set_yscale('symlog')
ax.set_xlim(0,1)
ax.set_xlabel('$q(l)$')
ax.set_ylabel('$S_{AA}(q)/CN$')


data = np.loadtxt('L_104/SKAA.dat')


ax.scatter(data[:,0],data[:,1])
ax.plot(data[:,0],data[:,1],label ='L=104')

ax.set_yscale('symlog')
ax.set_xlim(0,0.8)
ax.set_xlabel('$q(l)$')
ax.set_ylabel('$S_{AA}(q)/CN$')

# data = np.loadtxt('L_106/SKAA.dat')


# ax.scatter(data[:,0],data[:,1])
# ax.plot(data[:,0],data[:,1],label ='L=106')

# ax.set_yscale('symlog')
# ax.set_xlim(0,1)
# ax.set_xlabel('$q(l)$')
# ax.set_ylabel('$S_{AA}(q)/CN$')

plt.legend()
plt.savefig('/home/tquah/Figures/SKAA.pdf',dpi=300)