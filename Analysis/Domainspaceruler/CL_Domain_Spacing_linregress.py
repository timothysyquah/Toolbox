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

# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/C_2')
# os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/gausswidth_2.0_invzeta_1.0_3D_CL')
os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/C_1')
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
                plt.savefig(f'/home/tquah/Figures/nbb{str(temparray[0][2]+1)}nsc40_01.pdf',dpi = 300)

            # plt.errorbar(array[:,0],array[:,8+j])

        meanD0 = Predict_D0(array,Linregress)

        if meanD0==0:
            print('skiping...')
            continue
        # if i==2:
        # plt.plot(meanD0,array[0,4],'ok')

        #     continue

        
        temparray[0,-1] = meanD0
        arraylist.append(temparray)
arraycombine = np.vstack(arraylist)
