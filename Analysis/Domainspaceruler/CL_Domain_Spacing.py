#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May  5 12:50:08 2021

@author: tquah
"""

import matplotlib.pyplot as plt
import os
import glob
import numpy as np
import pandas as pd
import scipy.interpolate as interpolate
from scipy.optimize import bisect,newton,curve_fit
from scipy.stats import linregress,norm
import re 

plt.close('all')
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
# rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'
rcParams['font.style'] = 'normal'

rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 15

def ListPattern(string,list_):
    r = re.compile(string)
    newlist = list(filter(r.match,list_))
    return newlist

def extract_value(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string) #im lazy so i made a function 


def Create_CStress_Tensor(df):
    L = df['L']
    H = df['Hamiltonian.Real']
    list_of_stress = ['StressXX.Real','StressYY.Real','StressZZ.Real',\
                  'StressXY.Real','StressYZ.Real','StressXZ.Real']
    list_of_stress_error = ['StressXX.Real_Error','StressYY.Real_Error','StressZZ.Real_Error',\
                   'StressXY.Real_Error','StressYZ.Real_Error','StressXZ.Real_Error']
    list_of_stress_Var = ['StressXX.Real_Var','StressYY.Real_Var','StressZZ.Real_Var',\
                   'StressXY.Real_Var','StressYZ.Real_Var','StressXZ.Real_Var']

    stress_tensor_loc = [[[0,0]],[[1,1]],[[2    ,2]],\
                     [[0,1],[1,0]],[[1,2],[2,1]],[[0,2],[2,0]]]    
    cauchy_tensor = np.zeros((3,3))
    cauchy_tensor_error = np.zeros((3,3))
    cauchy_tensor_var = np.zeros((3,3))

    # print(list(df))
    for i in range(len(list_of_stress)):
        for j in range(len(stress_tensor_loc[i])):
            cauchy_tensor[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress[i]]
            cauchy_tensor_error[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress_error[i]]
            cauchy_tensor_var[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress_Var[i]]

    # eigenval,eigenvect = np.linalg.eig(cauchy_tensor)
    eigenval = np.array([cauchy_tensor[0,0],cauchy_tensor[1,1],cauchy_tensor[2,2]])
    eigenvect = 0 
    var_diagonal = np.array([cauchy_tensor_var[0,0],cauchy_tensor_var[1,1],cauchy_tensor_var[2,2]])
    sem_diagonal = np.array([cauchy_tensor_error[0,0],cauchy_tensor_error[1,1],cauchy_tensor_error[2,2]])

    return L,H,cauchy_tensor,eigenval,eigenvect,var_diagonal,sem_diagonal


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
        L, H, cauchy_tensor, eigenval, eigenvect,var_diagonal,sem_diagonal = Create_CStress_Tensor(df.iloc[i])
        check = Analyze_Principle_ST(eigenval,0.05)
        outputlist.append(np.hstack((np.array([L,H]),eigenval,check,var_diagonal,sem_diagonal)))
        data_list.append([L,cauchy_tensor,var_diagonal,sem_diagonal])
    outarray = np.vstack(outputlist)
    return outarray,data_list
def powerlawscale(x,a,b):
    return a*x**b

    
def Bisection(x,yarray1,yarray2):
    diff = yarray1-yarray2
    interpfxn = interpolate.interp1d(x,diff)
    root = bisect(interpfxn,np.min(x),np.max(x))
    return root

def Newton(x,yarray1,yarray2):
    diff = yarray1-yarray2
    interpfxn = interpolate.interp1d(x,diff)
    root = newton(interpfxn,np.mean(x))
    return root

def Linregress(x,yarray1,yarray2):
    m1, b1, r1, p1, se1 = linregress(x,yarray1)
    m2, b2, r2, p2, se2 = linregress(x,yarray2)
    return (b2-b1)/(m1-m2)
    
    

def Predict_D0(array,Intersection):
    try:
        roota = Intersection(array[:,0],array[:,2],array[:,3])
        rootb = Intersection(array[:,0],array[:,2],array[:,4])
        return (roota+rootb)/2
    except:
        return 0

# def Linregress(x,yarray1,yarray2):

def ParsePath(path,parameters):
    pathsplit = path.split('/')
    returnlist = []
    for i in range(len(parameters)):
        strg = ListPattern(parameters[i],pathsplit)
        value = extract_value(strg[0])
        returnlist.append(value)
    return np.hstack(returnlist)
    
    
    

infile = 'operators_result.dat'

# Check_Stress_Operator(infile)

# Isolate_Other(array)
desired_parameters = ['chi','nsc','nbb']

colors = ['g','b','r']
shape = ['^','d','s']
confidence = 0.05
confidence_interval = norm.ppf(1-confidence/2)
plotdebug = True
# plt.close('all')

#change path
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/main_test/gausswidth_2.0_invzeta_1.0_3D_CL')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/C_2/junkCL')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/C_2/')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/Unsmeared')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/finite_size_effects/nynz_646464')

# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/finite_size_effects/nx_notdouble')
# # os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/Unsmeared')
# # os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/main_test/c_1')

os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/oldC/C_2')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/gausswidth_2.0_invzeta_1.0_3D_CL')

IDIR = os.getcwd()
# print(os.listdir())
file_list = glob.glob('CL/**/nsc40.0/nbb*/**/operators_result.dat',recursive=True)


#dataarray 

arraylist = []

#



for i in range(len(file_list)):
    absinfile = os.path.join(IDIR,file_list[i])
    df = pd.read_csv(absinfile,delimiter=" ")
    if len(df.index)>1:
        print(file_list[i])
        print(absinfile)
        array_unsort,datalist = Primary_Stress_Tensor(absinfile)
        # if len(np.shape(array))==1:
        reorder = np.argsort(array_unsort[:,0])
        array = array_unsort[reorder,:]
        #     continuefig,ax = plt.figure()

        parray = ParsePath(file_list[i],desired_parameters)
        temparray = np.zeros((1,len(desired_parameters)+1))
        temparray[0,0:-1] = parray

        if plotdebug:
            plt.figure()
            plt.title('Nbb '+str(temparray[0][2]+1))

            for j in range(3):
                plt.errorbar(array[:,0],array[:,2+j],yerr = array[:,11+j]*confidence_interval, marker=shape[j],color=colors[j],linestyle='none',capsize=5.0)
                # plt.xlabel("$L_x(l)$")
                # plt.ylabel("$\sigma_{i,i}$")
                # plt.tight_layout()
                # if j ==0:
                #     ax[1].errorbar(array[:,0],array[:,2+j],yerr = array[:,8+j]*confidence_interval, marker=shape[j],color=colors[j],linestyle='none',capsize=5.0)
                #     ax[0].errorbar(array[:,0],array[:,2+j],yerr = array[:,11+j]*confidence_interval, marker=shape[j],color=colors[j],zorder = 10,linestyle='none',capsize=5.0)
                # else:
                #     ax[1].errorbar(array[:,0],array[:,2+j],yerr = array[:,8+j]*confidence_interval, marker=shape[j],color=colors[j],linestyle='none',capsize=5.0)
                #     ax[0].errorbar(array[:,0],array[:,2+j],yerr = array[:,11+j]*confidence_interval, marker=shape[j],color=colors[j],zorder = 10,linestyle='none',capsize=5.0)
                plt.ylabel("$\sigma_{i,i}$")
                plt.tight_layout()
                plt.savefig(f'/home/tquah/Figures/nbb{str(temparray[0][2]+1)}nsc40_00.pdf',dpi = 300)

            # plt.errorbar(array[:,0],array[:,8+j])

        meanD0 = Predict_D0(array,Newton)

        if meanD0==0:
            print('skiping...')
            continue
        # if i==2:
        # plt.plot(meanD0,array[0,4],'ok')

        #     continue

        
        temparray[0,-1] = meanD0
        arraylist.append(temparray)
                
# # plt.xlabel('$L_x(l)$')
# # plt.ylabel('$\sigma_{i,i}$')
# import matplotlib.patches as mpatches

# l1 = mpatches.Patch(color='g', label='$\sigma_{x,x}$')
# l2 = mpatches.Patch(color='r', label='$\sigma_{y,y}$')
# l3 = mpatches.Patch(color='b', label='$\sigma_{z,z}$')

# lines = [l1,l2,l3]
# labels = [line.get_label() for line in lines]

# # ax[0].legend(lines,labels)

# plt.xlabel("$L_x(l)$")
# plt.ylabel("$\sigma_{i,i}$",labelpad=15)

# plt.tight_layout()
# plt.savefig('/home/tquah/Figures/nbb80nsc40.pdf',dpi = 300)

# fullarray = np.vstack(arraylist)
# reorder = np.argsort(fullarray[:,2])
# fullarray[:,2] = fullarray[:,2]+1

# def reorder_array(arrayin,column):
#     reorder = np.argsort(arrayin[:,column])
#     return arrayin[reorder,:] 

# # SCFTPath = '2021-05-10/nsc40.0_chi0.1_f0.5_LD.dat'
# SCFTPath = '2021-05-10/nsc40.0_chi0.1_f0.5_LD.dat'

# scftarray = np.loadtxt(SCFTPath)
# scftarray[:,0] = scftarray[:,0]+1
# scftarray = reorder_array(scftarray,0)


# plt.figure()
# newarray = reorder_array(fullarray,2)

# maxvalue = np.max((newarray[:,2]))
# loc = np.where(scftarray[:,0]<=maxvalue)[0]
# # loc = np.array([1,3,6,9,13,17,21])
# scftarray_plot = scftarray[loc,:]



# xspace = np.linspace(np.min(newarray[:,2]),np.max((newarray[:,2])),100)
# cs = interpolate.CubicSpline(np.log(newarray[:,2]),np.log(newarray[:,3]))
# Dspace = cs(np.log(xspace))

# cs_scft = interpolate.CubicSpline(np.log(scftarray_plot[:,0]),np.log(scftarray_plot[:,1]))
# Dspace_scft = cs_scft(np.log(xspace))




# plt.plot(xspace,np.exp(Dspace),'r')
# plt.plot(xspace,np.exp(Dspace_scft),'--k')

# plt.loglog(newarray[:,2],newarray[:,3],'or')
# plt.loglog(scftarray_plot[:,0],scftarray_plot[:,1],'^k')

# plt.xticks([50,75,100,150,200])
# plt.yticks([60,70,80,100,120,150,200])
# plt.ylabel('$D/l$')
# plt.xlabel('$N_{bb}$')
    
# plt.tight_layout()
# plt.savefig('/home/tquah/Figures/40_NbbvsD_c2.pdf',dpi = 300)
# Dspace_dX = cs.derivative(1)

# arrayfit = (newarray[0:-1,2]+newarray[1:,2])/2

# DspaceDX = Dspace_dX(np.log(arrayfit))

# Dspacescft_dX = cs_scft.derivative(1)
# DspacescftDX = Dspacescft_dX(np.log(arrayfit))

# plt.figure()
# plt.plot(arrayfit,DspaceDX,'r')
# plt.plot(arrayfit,DspacescftDX,'k')
# plt.plot(xspace,np.ones_like(xspace)*2/3,'--k')
# # for point in newarray[:,2]:
# #     plt.axvline(point,color= 'k',linestyle='--',alpha = 0.5)



# plt.ylabel('$\gamma$')
# plt.xlabel('$N_{bb}$')
# plt.tight_layout()
# plt.savefig('/home/tquah/Figures/gamma_40c2.pdf',dpi = 300)

# popt,pcov = curve_fit(powerlawscale,newarray[-2:,2],newarray[-2:,3])


#         # meanD0 = Predict_D0(array)
#         # print(meanD0)


# file_list = glob.glob('CL/**/nsc40.0/nbb*/**/SKAA.dat',recursive=True)
# plt.figure()
# for i in range(0,len(file_list),1):
#     dataarray = np.loadtxt(file_list[i])
#     plt.plot(dataarray[:,0],dataarray[:,1])    
    
    
    # array,datalist = Primary_Stress_Tensor(absinfile)


# # pathlist = ['/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_64_64_64/Period_1',\
# #           '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_32_32_32/Period_1',\
# #               '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2/Period_1']

# # larray = [64,32,16]
# # fig, ax = plt.subplots(1, 3,sharex=True, sharey=True,figsize=(10,5))

# # for j in range(len(pathlist)): 
# #     os.chdir(pathlist[j])
# #     array,datalist = Primary_Stress_Tensor(infile)
# #     ax[j].title.set_text(f'$L_y=L_z = {larray[j]}$')
    
# #     colors = ['g','b','r']
# #     shape = ['^','d','s']
# #     text = ['$\sigma_1$','$\sigma_2$','$\sigma_3$']
# #     yaverage = np.linspace(np.min(array[:,2])*1.5,np.max(array[:,2])*1.5,10)
# #     for i in range(3):
# #         ax[j].plot(array[:,0],array[:,2+i],colors[i])
# #         ax[j].scatter(array[:,0],array[:,2+i],marker=shape[i],color=colors[i],label =text[i] )
    
# #     # plt.ylim((np.min(array[:,2])*1.1,np.max(array[:,2])*1.1))
# #     meanD0 = Predict_D0(array)
    
# #     xaverage = np.ones_like(yaverage)*meanD0
# #     ax[j].plot(xaverage,yaverage,'--k')

# # plt.setp(ax, xticks=array[::2,0])
# # fig.add_subplot(111, frameon=False)
# # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# # plt.ylabel('$\sigma_i$')
# # plt.xlabel('$L(l)$')
# # ax[0].legend()
# # plt.tight_layout()
# # plt.savefig('/home/tquah/Figures/domainstress.png',dpi=300)



# # plt.close('all')

# # pathlist = ['/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2/Period_1',\
# #               '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_long/Period_1',\
# #                    '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_dt/Period_1']

# # textlabel = ['Time Steps = $5x10^5$','Time Steps = $1x10^6$','Time Steps = $1x10^6$*']
# # fig, ax = plt.subplots(1, 3,sharex=True, sharey=True,figsize=(10,5))

# # for j in range(len(pathlist)): 
# #     os.chdir(pathlist[j])
# #     array,datalist = Primary_Stress_Tensor(infile)
# #     ax[j].title.set_text(textlabel[j])
    
# #     colors = ['g','b','r']
# #     shape = ['^','d','s']
# #     text = ['$\sigma_1$','$\sigma_2$','$\sigma_3$']
# #     yaverage = np.linspace(np.min(array[:,2])*1.5,np.max(array[:,2])*1.5,10)
# #     for i in range(3):
# #         ax[j].plot(array[:,0],array[:,2+i],colors[i])
# #         ax[j].scatter(array[:,0],array[:,2+i],marker=shape[i],color=colors[i],label =text[i] )
    
# #     # plt.ylim((np.min(array[:,2])*1.1,np.max(array[:,2])*1.1))
# #     meanD0 = Predict_D0(array)
    
# #     xaverage = np.ones_like(yaverage)*meanD0
# #     ax[j].plot(xaverage,yaverage,'--k')

# # plt.setp(ax, xticks=array[::2,0])
# # fig.add_subplot(111, frameon=False)
# # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# # plt.ylabel('$\sigma_i$')
# # plt.xlabel('$L(l)$')
# # ax[0].legend()
# # plt.tight_layout()
# # plt.savefig('/home/tquah/Figures/longaveragestress.png',dpi=300)

# plt.close('all')

# pathlist = ['/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_64_64_64/Period_1',
#               '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_64_64_64/Period_2']

# textlabel = ['$L_y = L_z = 64$ $l$ : Period = 1','$L_y = L_z = 64$ $l$: Period = 2']
# fig, ax = plt.subplots(1, 2,sharex=False, sharey=True,figsize=(10,5))

# for j in range(len(pathlist)): 
#     os.chdir(pathlist[j])
#     array,datalist = Primary_Stress_Tensor(infile)
#     ax[j].title.set_text(textlabel[j])
    
#     colors = ['k','b','m']
#     shape = ['^','d','s']
#     text = ['$\sigma_1$','$\sigma_2$','$\sigma_3$']
#     yaverage = np.linspace(-10,10,10)
#     for i in range(3):
#         ax[j].plot(array[:,0],array[:,2+i],colors[i])
#         ax[j].scatter(array[:,0],array[:,2+i],marker=shape[i],color=colors[i],label =text[i] )
    
#     # plt.ylim((np.min(array[:,2])*1.1,np.max(array[:,2])*1.1))
#     meanD0 = Predict_D0(array)
#     print(meanD0/2)
#     xaverage = np.ones_like(yaverage)*meanD0
#     ax[j].plot(xaverage,yaverage,'--k')
#     # ax[j].set_xticks(array[::2,0])
#     ax[j].set_xlim(np.min(array[:,0])-1,np.max(array[:,0])+1)
#     ax[j].set_ylim(np.min(array[:,2])-0.001,np.max(array[:,2])+0.001)
# plt.tight_layout()

# fig.add_subplot(111, frameon=False)
# plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# plt.ylabel('$\sigma_i$')
# plt.xlabel('$L_x$ $(l)$')
# ax[0].legend()
# plt.tight_layout()
# plt.savefig('/home/tquah/Figures/646464_Period_1_2.png',dpi=300)

# # plt.close('all')

# # pathlist = ['/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period_long/Period_2',\
# #               '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2_period/Period_2',\
# #                    '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/DougMethod/Doug_Stress_Method/Study_Increase_Resolution_Along_z/L_npw_2/Period_2']

# # textlabel = ['Time Steps = $1x10^6$','Time Steps = $5x10^5$','Time Steps = $5x10^5$*']
# # fig, ax = plt.subplots(1, 3,sharex=True, sharey=False,figsize=(10,5))

# # for j in range(len(pathlist)): 
# #     os.chdir(pathlist[j])
# #     array,datalist = Primary_Stress_Tensor(infile)
# #     ax[j].title.set_text(textlabel[j])
    
# #     colors = ['g','b','r']
# #     shape = ['^','d','s']
# #     text = ['$\sigma_1$','$\sigma_2$','$\sigma_3$']
# #     yaverage = np.linspace(np.min(array[:,2])*1.5,np.max(array[:,2])*1.5,10)
# #     for i in range(3):
# #         ax[j].plot(array[:,0],array[:,2+i],colors[i])
# #         ax[j].scatter(array[:,0],array[:,2+i],marker=shape[i],color=colors[i],label =text[i] )
    
# #     # plt.ylim((np.min(array[:,2])*1.1,np.max(array[:,2])*1.1))
# #     meanD0 = Predict_D0(array)
    
# #     xaverage = np.ones_like(yaverage)*meanD0
# #     ax[j].plot(xaverage,yaverage,'--k')

# # plt.setp(ax, xticks=array[::2,0])
# # fig.add_subplot(111, frameon=False)
# # plt.tick_params(labelcolor='none', which='both', top=False, bottom=False, left=False, right=False)
# # plt.ylabel('$\sigma_i$')
# # plt.xlabel('$L(l)$')
# # ax[0].legend()
# # plt.tight_layout()
# # plt.savefig('/home/tquah/Figures/periodstress.png',dpi=300)
