#!/usr/bin/env python3

import numpy as np
import os
import glob
import matplotlib.pyplot as plt


IDIR = os.getcwd()
os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/prototypedir/NSC20/nbb9/CL/TimeConvergeData')
file_list = glob.glob('sym*.dat')

plt.close('all')
dimension = [r'$L=32 \frac{b}{\sqrt{6}}$',r'$L=64 \frac{b}{\sqrt{6}}$',r'$L=16 \frac{b}{\sqrt{6}}$']
color = ['r','k','b']
symbol = ['o','^','s']
fig = plt.figure()
ax = plt.gca()
datalist = []
for i in range(len(file_list)):
    data = np.loadtxt(file_list[i])
    if i==2:
        
        datalist.append(data[:-1])
    else:
        datalist.append(data)
    ax.scatter(1/data[:,0],data[:,1],label = dimension[i],marker=symbol[i],color = color[i])
    ax.errorbar(1/data[:,0],data[:,1],yerr=data[:,2],ls='none')
    
    
ax.set_xscale('log')
ax.set_xlabel('$1/\Delta t$')
ax.set_ylabel(r'$Re(\left< H \right>)$')
ax.legend()
plt.tight_layout()
plt.savefig('sym1.pdf',dpi = 300)

fig = plt.figure()

ax = plt.gca()
for i in range(len(file_list)):
    # ax.scatter(1/datalist[i][:,0],datalist[i][:,1]-datalist[1][-1,1]*np.ones_like(datalist[i][:,1]))
    ax.scatter(1/datalist[i][:,0],datalist[i][:,1]-4.9999999948e-02*np.ones_like(datalist[i][:,1]),label = dimension[i],marker=symbol[i],color = color[i] )

    ax.errorbar(1/datalist[i][:,0],datalist[i][:,1]-4.9999999948e-02*np.ones_like(datalist[i][:,1]),yerr=datalist[i][:,2],ls='none')
    
ax.legend()

ax.set_xscale('log')
ax.set_xlabel('$1/\Delta t$')
ax.set_ylabel(r'$Re(\left< H_{CL} \right>)-Re(\left< H_{SCFT} \right>)$')
plt.tight_layout()

plt.savefig('sym2.pdf',dpi = 300)
os.chdir(IDIR)