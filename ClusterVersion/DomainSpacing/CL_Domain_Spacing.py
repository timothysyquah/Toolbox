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
from scipy.stats import linregress
import re 
import argparse
from scipy.optimize import newton,bisect
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
rcParams['axes.titlesize'] = 20
def ListPattern(string,list_):
    r = re.compile(string)
    newlist = list(filter(r.match,list_))
    return newlist
def extract_value(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string) #im lazy so i made a function 
def Create_CStress_Tensor(df):
    L = df['L']
    H = df['Hamiltonian.Real']
    Herr = df['Hamiltonian.Real_Error']
    list_of_stress = ['StressXX.Real','StressYY.Real','StressZZ.Real',\
                  'StressXY.Real','StressYZ.Real','StressXZ.Real']
    list_of_stress_error = ['StressXX.Real_Error','StressYY.Real_Error','StressZZ.Real_Error',\
                   'StressXY.Real_Error','StressYZ.Real_Error','StressXZ.Real_Error']

    stress_tensor_loc = [[[0,0]],[[1,1]],[[2    ,2]],\
                     [[0,1],[1,0]],[[1,2],[2,1]],[[0,2],[2,0]]]    
    cauchy_tensor = np.zeros((3,3))
    cauchy_tensor_error = np.zeros((3,3))

    # print(list(df))
    for i in range(len(list_of_stress)):
        for j in range(len(stress_tensor_loc[i])):
            cauchy_tensor[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress[i]]
            cauchy_tensor_error[stress_tensor_loc[i][j][0],stress_tensor_loc[i][j][1]] = df[list_of_stress_error[i]]

    # eigenval,eigenvect = np.linalg.eig(cauchy_tensor)
    eigenval = np.array([cauchy_tensor[0,0],cauchy_tensor[1,1],cauchy_tensor[2,2]])
    eigenvect = 0 
    sem_diagonal = np.array([cauchy_tensor_error[0,0],cauchy_tensor_error[1,1],cauchy_tensor_error[2,2]])
    return L,H,Herr,cauchy_tensor,eigenval,eigenvect,sem_diagonal
def Analyze_Principle_ST(eigenvalues,relative_tolerance):
    loc = [[0,1],[0,2],[1,2]]
    check = np.zeros(3)
    for i in range(len(loc)):
        check[i] = np.isclose(eigenvalues[loc[i]][0],eigenvalues[loc[i]][1],rtol = relative_tolerance)
    return check
def Primary_Stress_Tensor(infile):
    df = pd.read_csv(infile,delimiter=" ",index_col=False)
    shape = df.shape
    # primary_
    outputlist = []
    for i in range(shape[0]):
        L, H,Herr, cauchy_tensor, eigenval, eigenvect,sem_diagonal = Create_CStress_Tensor(df.iloc[i])
        check = Analyze_Principle_ST(eigenval,0.05)
        outputlist.append(np.hstack((np.array([L,H,Herr]),eigenval,check,sem_diagonal)))
    outarray = np.vstack(outputlist)
    return outarray
def bracket(array):
    c1 = np.sign(array[:,2]-array[:,0])
    cross = np.where(c1[:-1]-c1[1:]>0)[0]
    if len(cross)>1:
        print('Warning...Multiple Crossing Points...Picking First One')
    
    if len(cross)==0:
        print('No Crossing Point')
        return np.array([0,None])
    else:
        pointcross = cross[0]    

        if pointcross==0:
            bracket = np.array([0,3])
        elif pointcross==len(c1):
            bracket = np.array([-3,None])
        else:
            bracket = np.array([pointcross-1,pointcross+2])
        return bracket
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
def Predict_D0(x0,array,Intersection):
    try:
        roota = Intersection(x0,array[:,0],array[:,1])
        rootb = Intersection(x0,array[:,0],array[:,2])
        return (roota+rootb)/2
    except:
        return 0

def ParsePath(path,parameters):
    pathsplit = path.split('/')
    returnlist = []
    for i in range(len(parameters)):
        strg = ListPattern(parameters[i],pathsplit)
        value = extract_value(strg[0])
        returnlist.append(value)
    return np.hstack(returnlist)
FUNCTION_MAP = {'Bisection' : Bisection,\
                    'Newton' : Newton,\
                    'Linregress' : Linregress}
if __name__ == '__main__':
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to sweep bottlebrushes')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default="CL/**/chi0.2/nsc10.0/nbb*/**/",help='list of directories that contain each phase point')
    parser.add_argument('-f','--file',action = 'store',default = 'operators_result.dat', help = 'File to read with averages and error', type = str)
    parser.add_argument('-s','--plotdebug',action = 'store',default = False, help = 'Stress Plot', type = bool)
    parser.add_argument('-hh','--plothamiltonian',action = 'store',default = False, help = 'Hamiltonian Plot', type = bool)
    parser.add_argument('-ep','--exportpath_plot',action = 'store',default = './', help = 'Path to Save Figure', type = str)
    parser.add_argument('-e','--exportname',action = 'store',default = 'domain_data.dat', help = 'Path to Save Data', type = str)
    parser.add_argument('-dp', '--desired_parameters', action='store',nargs='+', default=['chi','nsc','nbb'],help='Desired parameters',type = str)
    parser.add_argument('-c','--command',default = 'Linregress', choices=FUNCTION_MAP.keys())
    args = parser.parse_args()
    path = os.path.join(args.dirs,args.file)
    colors = ['g','b','r']
    shape = ['^','d','s']
    IDIR = os.getcwd()
    file_list = glob.glob(path,recursive=True)
    arraylist = []
    func = FUNCTION_MAP[args.command]

    if len(file_list)==0:
        print('Warning...No Files')
    else:
        for i in range(len(file_list)):
            absinfile = os.path.join(IDIR,file_list[i])
            df = pd.read_csv(absinfile,delimiter=" ",index_col=False)
            if len(df.index)>1:
                array_unsort = Primary_Stress_Tensor(absinfile)
                # if len(np.shape(array))==1:
                reorder = np.argsort(array_unsort[:,0])
                array = array_unsort[reorder,:]
                #     continuefig,ax = plt.figure()
                parray = ParsePath(file_list[i],args.desired_parameters)
                temparray = np.zeros((1,len(args.desired_parameters)+1))
                temparray[0,0:-1] = parray
                bounds = bracket(array[:,3:6])
                meanD0 = Predict_D0(array[bounds[0]:bounds[1],0],array[bounds[0]:bounds[1],3:6],func)
                if args.plotdebug:
                    plt.figure()
                    plt.title('Nbb '+str(temparray[0][2]+1))
                    for j in range(3):
                        plt.errorbar(array[:,0],array[:,3+j],yerr = array[:,9+j], marker=shape[j],color=colors[j],linestyle='none',capsize=5.0)
                        plt.xlabel(r'Cell Size $(L_x)$')
                        plt.ylabel(r"Stress Components ($\left < \sigma_{i,i} \right>$)")
                    plt.tight_layout()
                    plt.scatter(meanD0,np.mean(array[bounds[0]:bounds[1],4]),color = 'k')
                    name = f'nbb{str(temparray[0][2]+1)}nsc{str(temparray[0][1])}.pdf'
                    exportpath = os.path.join(args.exportpath_plot,name)
                    plt.savefig(exportpath,dpi = 300)
                if args.plothamiltonian:
                    plt.figure()
                    plt.title('Nbb '+str(temparray[0][2]+1))
                    plt.errorbar(array[:,0],array[:,1],yerr = array[:,2], marker=shape[j],color=colors[j],linestyle='none',capsize=5.0)
                    plt.xlabel(r'Cell Size $(L_x)$')
                    plt.ylabel(r"$\left< H[w_+,w_-] \right>$")
                    plt.tight_layout()
                    name = f'H_nbb{str(temparray[0][2]+1)}nsc{str(temparray[0][1])}.pdf'
                    exportpath = os.path.join(args.exportpath_plot,name)
                    plt.savefig(exportpath,dpi = 300)
                if meanD0==0:
                    print('skiping...')
                    continue
                temparray[0,-1] = meanD0
                arraylist.append(temparray)
        data_array = np.vstack(arraylist)
        data_array[:,2] = data_array[:,2]+1
        reorder = np.argsort(data_array[:,2])
        data_array = data_array[reorder,:]
        np.savetxt(args.exportname,data_array)