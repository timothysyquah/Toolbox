#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun  2 10:20:44 2021

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import glob
os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/C_2/CL/f0.5/chi0.1/nsc40.0/nbb79.0/LAM3DPhase')

os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/C_2/CL/f0.5/chi0.1/nsc40.0/nbb129.0/LAM3DPhase')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/finite_size_effects/nx_double/CL/f0.5/chi0.1/nsc40.0/nbb79.0/LAM3DPhase')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/finite_size_effects/nynz_646464/CL/f0.5/chi0.1/nsc40.0/nbb79.0/LAM3DPhase')

plt.close('all')
fig = plt.figure()
ax = plt.gca()
listdir = glob.glob('L*')
listdir = [listdir[0]]
color = ['r']

for i in range(0,len(listdir),1):
    fp = float(listdir[i][2:])
    data = np.loadtxt(f'{listdir[i]}/SKAA.dat')

    ax.scatter(data[:,0],data[:,1],color = color[i])
    ax.plot(data[:,0],data[:,1],label =f'L={fp}',color = color[i],alpha = 0.5)
    
# ax.set_xlim(0,2)

ax.set_yscale('symlog')
ax.set_xlabel('$q(l)$')
ax.set_ylabel('$S_{AA}(q)/CN$')



plt.legend()
plt.tight_layout()
# plt.savefig('/home/tquah/Figures/SKAA_cutoff_improve.pdf',dpi=300)


# plt.close('all')
# fig = plt.figure()
# ax = plt.gca()
os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/finite_size_effects/nx_notdouble/CL/f0.5/chi0.1/nsc40.0/nbb79.0/LAM3DPhase')
listdir = glob.glob('L*')
listdir = [listdir[0]]

color = ['k']

for i in range(0,len(listdir),1):
    fp = float(listdir[i][2:])
    data = np.loadtxt(f'{listdir[i]}/SKAA.dat')

    ax.scatter(data[:,0],data[:,1],color = color[i],marker = '^')
    ax.plot(data[:,0],data[:,1],label =f'L={fp}, $\Delta x=1$',color = color[i],alpha = 0.5)

    
# ax.set_xlim(0,2)

ax.set_yscale('symlog')
ax.set_xlabel('$q(l)$')
ax.set_ylabel('$S_{AA}(q)/CN$')



plt.legend()
plt.tight_layout()
plt.savefig('/home/tquah/Figures/SKAA_cutoff_compare_0.pdf',dpi=300)


# plt.close('all')
# listdir = glob.glob('L*')



# for i in range(0,len(listdir),1):
#     fig = plt.figure()
#     ax = plt.gca()

#     fp = float(listdir[i][2:])
#     data = np.loadtxt(f'{listdir[i]}/operators.dat')

#     # ax.scatter(data[:,0],data[:,4])
#     ax.plot(data[:,0],data[:,4],label ='$\sigma_{x,x}$')
#     ax.plot(data[:,0],data[:,6],label ='$\sigma_{y,y}$')
#     ax.plot(data[:,0],data[:,8],label ='$\sigma_{z,z}$')
#     ax.set_xlabel('Step')
#     ax.set_ylabel('$\Im(\sigma_{i,i})$')
#     plt.legend()
#     plt.tight_layout()
#     plt.savefig(f'/home/tquah/Figures/operators{fp}.pdf',dpi=300)

# ax.set_xlim(0,1)
# ax.set_yscale('symlog')



