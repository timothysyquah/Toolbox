#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr  4 11:50:14 2020

@author: tquah
"""

import matplotlib.pyplot as plt
from matplotlib import cm
import numpy as np
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401 unused import
import matplotlib as mpl
path_export = '/home/tquah/ParaviewInport/GyroidProblems/N4f0.3_AGYR32.csv'

data = np.loadtxt(path_export,skiprows = 1,delimiter = ',')

data[:,3] = data[:,3]/np.max(data[:,3])
data[:,4] = data[:,4]/np.max(data[:,3])
data[:,5] = data[:,5]/np.max(data[:,3])

a = np.array([10,10,10,0,0,0])
data = np.vstack((data,a))

a = np.array([-10,-10,-10,1,1,1])
data = np.vstack((data,a))

alphaset = 0.1
colorplot1 = []
colorplot2= []
ind1 = []
ind2 = []
for i in range(0,len(data[:,4])):
    if data[i,3] or data[i,5]>1e-5:
        if data[i,3]>data[i,5]:
            colorplot1.append(data[i,3])
            ind1.append(i)
        if data[i,3]>data[i,5]:
            colorplot2.append(data[i,5])
            ind2.append(i)

plt.close('all')
m= 'o'
fig = plt.figure(figsize=plt.figaspect(0.25))
ax = fig.add_subplot(1, 3, 1, projection='3d')
norm = mpl.colors.Normalize(vmin=0.,vmax=1.)
colors_1 = cm.Reds(colorplot1)
colors_2 = cm.Blues(colorplot2)
ax.scatter(data[ind1,0], data[ind1,1], data[ind1,2], marker=m,alpha = alphaset,c = colors_1)
ax.scatter(data[ind2,0], data[ind2,1], data[ind2,2], marker=m,alpha = alphaset,c = colors_2)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

ax = fig.add_subplot(1, 3, 2, projection='3d')
norm = mpl.colors.Normalize(vmin=0.,vmax=1.)
colors_1 = cm.Reds(colorplot1)
ax.scatter(data[ind1,0], data[ind1,1], data[ind1,2], marker=m,alpha = alphaset,c = colors_1)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])

ax = fig.add_subplot(1, 3, 3, projection='3d')
norm = mpl.colors.Normalize(vmin=0.,vmax=1.)
colors_2 = cm.Blues(colorplot2)
ax.scatter(data[ind2,0], data[ind2,1], data[ind2,2], marker=m,alpha = alphaset,c = colors_2)
ax.grid(False)
ax.set_xticks([])
ax.set_yticks([])
ax.set_zticks([])
plt.savefig('test.png',dpi = 300)
# fig = plt.figure()
# ax = fig.add_subplot(111, projection='3d')
# colorplot = np.zeros_like(data[:,4])
# for i in range(0,len(data[:,4])):
#     if data[i,4]>data[i,5]:
#         colorplot[i] = -data[i,4]
#     if data[i,5]>data[i,4]:
#         colorplot[i] = data[i,5]


# colors_1 = cm.bwr(colorplot)

# ax.scatter(data[:,0], data[:,1], data[:,2], marker=m,alpha = alphaset,c = colors_1)
# ax.grid(False)
# ax.set_xticks([])
# ax.set_yticks([])
# ax.set_zticks([])
# plt.axis('off')
