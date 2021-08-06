#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 11:53:52 2021

@author: tquah
"""

import numpy as np
import re
import pandas as pd 
from scipy import interpolate
from scipy.optimize import bisect,newton
from scipy.stats import linregress

def powlaw(x,a,b):
    return (a*x**b)


def convert_string_float(lst):
    return [float(i) for i in lst]

def determine_directory_structure(listoffile):
    length = []
    directory_num = []
    for file in listoffile:
        param_list = re.findall(r"[-+]?\d*\.\d+|\d+", file)
        directory_num.append(convert_string_float(param_list))
        length.append(len(param_list))
        
        
    directory_level = max(set(length), key = length.count)
    
    purged_list =[i for i in directory_num if len(i)==directory_level]
    purged_file_list =[i for i in listoffile if len(re.findall(r"[-+]?\d*\.\d+|\d+", i))==directory_level]
    array = np.vstack(purged_list)
    dyanamic = []
    for i in range(0,directory_level,1):
        if len(np.unique(array[:,i]))>1:
            dyanamic.append(i)
    return dyanamic,purged_list,array,purged_file_list


def get_domain_spacing_SCFT_simple(file):
    path_list = file.split('/')
    pathstatus = ''
    for i in range(0,len(path_list)-1,1):
        pathstatus+=path_list[i]+'/'
    
    statuspath = pathstatus+'STATUS'

    co = open(statuspath,'r')
    statusread=int(co.read())
    co.close()

    # operator_path = pathstatus+'operators.dat'
    # load_data = np.loadtxt(operator_path)[:,1]
    # print(load_data)
    # store = np.min(load_data)

    if int(statusread)!=2:
        # print("Convergance Fail")
        return False
    else:
        fo = open(file,'r')
        content = fo.read().split('\n')
        fo.close()
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        result = float(r[r.find("(")+1:r.find(")")])
        return [result]



def fill_dictionaries(listoffile,desiredparameters,paramarray ,extractmethod,cutoffpt = 300):
    count=0
    data_dict = dict()
    for file in listoffile:
        t = ()
        for i in range(0,len(desiredparameters)-1,1):
            var = paramarray[count,desiredparameters[i]]
            t = t + (var,)
        if get_domain_spacing_SCFT_simple(file)!=False:
            if t not in list(data_dict):
                    data_dict[t] = []
                    data_dict[t].append([paramarray[count,desiredparameters[-1]]]+get_domain_spacing_SCFT_simple(file))
            else:
                    data_dict[t].append([paramarray[count,desiredparameters[-1]]]+get_domain_spacing_SCFT_simple(file))            
        count+=1    
        
    new_data_dict = dict()
    for header in list(data_dict):
        array = np.vstack(data_dict[header])
        newsort = np.argsort(array[:,0])
        newarray  = array[newsort,:]
        cutoff = np.where(newarray[:,0]<cutoffpt)[0]
        new_data_dict[header] = newarray[cutoff,:]
    return new_data_dict

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

    return L,H,cauchy_tensor,eigenval,eigenvect,sem_diagonal    

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
        L, H, cauchy_tensor, eigenval, eigenvect,sem_diagonal = Create_CStress_Tensor(df.iloc[i])
        check = Analyze_Principle_ST(eigenval,0.05)
        outputlist.append(np.hstack((np.array([L,H]),eigenval,check,sem_diagonal)))
        data_list.append([L,cauchy_tensor,sem_diagonal])
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
    


