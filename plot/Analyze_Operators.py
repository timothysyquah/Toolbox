#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 16 14:46:14 2021

@author: tquah
"""


import matplotlib.pyplot as plt
import os
import glob
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from scipy.optimize import bisect

plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
# rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'
rcParams['font.style'] = 'normal'

rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 11
\

def Create_CStress_Tensor(df):
    L = df['L']
    H = df['Hamiltonian.Real']
    list_of_stress = ['StressXX.Real','StressYY.Real','StressZZ.Real',\
                  'StressXY.Real','StressYZ.Real','StressXZ.Real']
    list_of_stress_error = ['StressXX.Real_Error','StressYY.Real_Error','StressZZ.Real_Error',\
                   'StressXY.Real_Error','StressYZ.Real_Error','StressXZ.Real_Error']

    stress_tensor_loc = [[[0,0]],[[1,1]],[[2,2]],\
                     [[0,1],[1,0]],[[1,2],[2,1]],[[0,2],[2,0]]]    
    cauchy_tensor = np.zeros((3,3))
    cauchy_tensor_error = np.zeros((3,3))
    print(list(df))
    for i in range(len(list_of_stress)):
        for j in range(len(stress_tensor_loc[i])):
            cauchy_tensor[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress[i]]
            cauchy_tensor_error[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress_error[i]]

    eigenval,eigenvect = np.linalg.eig(cauchy_tensor)

    return L,H,cauchy_tensor,eigenval,eigenvect,cauchy_tensor_error


def random_Cauchy_tensor(cauchytensor,cauchyerror):
    newcauchy_tensor = np.zeros_like(cauchytensor)
    
    stress_tensor_loc = [[[0,0]],[[1,1]],[[2,2]],\
                     [[0,1],[1,0]],[[1,2],[2,1]],[[0,2],[2,0]]]    
    for i in range(len(stress_tensor_loc)):
        for j in range(len(stress_tensor_loc[i])):
            print(i,j)
    


def Analyze_Principle_ST(eigenvalues,relative_tolerance):
    loc = [[0,1],[0,2],[1,2]]
    check = np.zeros(3)
    for i in range(len(loc)):
        check[i] = np.isclose(eigenvalues[loc[i]][0],eigenvalues[loc[i]][1],rtol = relative_tolerance)
    return check


def Primary_Stress_Tensor(infile):
    df = pd.read_csv(infile,delimiter=" ")
    shape = df.shape
    # primary_
    data_list = []
    outputlist = []
    for i in range(shape[0]):
        L, H, cauchy_tensor, eigenval, eigenvect,cauchy_tensor_error = Create_CStress_Tensor(df.iloc[i])
        check = Analyze_Principle_ST(eigenval,0.01)
        outputlist.append(np.hstack((np.array([L,H]),eigenval,check)))
        data_list.append([L,cauchy_tensor,cauchy_tensor_error])
    outarray = np.vstack(outputlist)
    return outarray,data_list

    
def IntersectionL1L2(x,yarray1,yarray2):
    diff = yarray1-yarray2
    interpfxn = interpolate.interp1d(x,diff)
    root = bisect(interpfxn,np.min(x),np.max(x))
    return root

def Predict_D0(array):
    roota = IntersectionL1L2(array[:,0],array[:,2],array[:,3])
    rootb = IntersectionL1L2(array[:,0],array[:,2],array[:,4])
    return (roota+rootb)/2

    
    
    

infile = 'operators_result.dat'

# Check_Stress_Operator(infile)

# Isolate_Other(array)




plt.close('all')

pathlist = ['/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_64_64_64/Period_1',\
          '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_32_32_32/Period_1',\
              '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2/Period_1']

larray = [64,32,16]
fig, ax = plt.subplots(1, 3,sharex=True, sharey=True,figsize=(10,5))

for j in range(len(pathlist)): 
    os.chdir(pathlist[j])
    array,datalist = Primary_Stress_Tensor(infile)
    ax[j].title.set_text(f'$L_y=L_z = {larray[j]}$')
    
    colors = ['g','b','r']
    shape = ['^','d','s']
    text = ['$\sigma_1$','$\sigma_2$','$\sigma_3$']
    yaverage = np.linspace(np.min(array[:,2])*1.5,np.max(array[:,2])*1.5,10)
    for i in range(3):
        ax[j].plot(array[:,0],array[:,2+i],colors[i])
        ax[j].scatter(array[:,0],array[:,2+i],marker=shape[i],color=colors[i],label =text[i] )
    
    # plt.ylim((np.min(array[:,2])*1.1,np.max(array[:,2])*1.1))
    meanD0 = Predict_D0(array)
    
    xaverage = np.ones_like(yaverage)*meanD0
    ax[j].plot(xaverage,yaverage,'--k')

plt.setp(ax, xticks=array[::2,0])
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.ylabel('$\sigma_i$')
plt.xlabel('$L(l)$')
ax[0].legend()
plt.tight_layout()
plt.savefig('/home/tquah/Figures/domainstress.png',dpi=300)



plt.close('all')

pathlist = ['/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2/Period_1',\
              '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_long/Period_1',\
                   '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_dt/Period_1']

textlabel = ['Time Steps = $5x10^5$','Time Steps = $1x10^6$','Time Steps = $1x10^6$*']
fig, ax = plt.subplots(1, 3,sharex=True, sharey=True,figsize=(10,5))

for j in range(len(pathlist)): 
    os.chdir(pathlist[j])
    array,datalist = Primary_Stress_Tensor(infile)
    ax[j].title.set_text(textlabel[j])
    
    colors = ['g','b','r']
    shape = ['^','d','s']
    text = ['$\sigma_1$','$\sigma_2$','$\sigma_3$']
    yaverage = np.linspace(np.min(array[:,2])*1.5,np.max(array[:,2])*1.5,10)
    for i in range(3):
        ax[j].plot(array[:,0],array[:,2+i],colors[i])
        ax[j].scatter(array[:,0],array[:,2+i],marker=shape[i],color=colors[i],label =text[i] )
    
    # plt.ylim((np.min(array[:,2])*1.1,np.max(array[:,2])*1.1))
    meanD0 = Predict_D0(array)
    
    xaverage = np.ones_like(yaverage)*meanD0
    ax[j].plot(xaverage,yaverage,'--k')

plt.setp(ax, xticks=array[::2,0])
fig.add_subplot(111, frameon=False)
plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
plt.ylabel('$\sigma_i$')
plt.xlabel('$L(l)$')
ax[0].legend()
plt.tight_layout()
plt.savefig('/home/tquah/Figures/longaveragestress.png',dpi=300)


# plt.close('all')

# pathlist = ['/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_long/Period_2',\
#               '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period/Period_2',\
#                    '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2/Period_2']

# textlabel = ['Time Steps = $1x10^6$','Time Steps = $5x10^5$','Time Steps = $5x10^5$*']
# fig, ax = plt.subplots(1, 3,sharex=True, sharey=False,figsize=(10,5))

# for j in range(len(pathlist)): 
#     os.chdir(pathlist[j])
#     array,datalist = Primary_Stress_Tensor(infile)
#     ax[j].title.set_text(textlabel[j])
    
#     colors = ['g','b','r']
#     shape = ['^','d','s']
#     text = ['$\sigma_1$','$\sigma_2$','$\sigma_3$']
#     yaverage = np.linspace(np.min(array[:,2])*1.5,np.max(array[:,2])*1.5,10)
#     for i in range(3):
#         ax[j].plot(array[:,0],array[:,2+i],colors[i])
#         ax[j].scatter(array[:,0],array[:,2+i],marker=shape[i],color=colors[i],label =text[i] )
    
#     # plt.ylim((np.min(array[:,2])*1.1,np.max(array[:,2])*1.1))
#     meanD0 = Predict_D0(array)
    
#     xaverage = np.ones_like(yaverage)*meanD0
#     ax[j].plot(xaverage,yaverage,'--k')

# plt.setp(ax, xticks=array[::2,0])
# fig.add_subplot(111, frameon=False)
# plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# plt.ylabel('$\sigma_i$')
# plt.xlabel('$L(l)$')
# ax[0].legend()
# plt.tight_layout()
# plt.savefig('/home/tquah/Figures/periodstress.png',dpi=300)
