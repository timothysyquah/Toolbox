#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov  5 11:35:01 2020

@author: tquah
"""



import numpy as np
import glob,re
import os,sys
import matplotlib.pyplot as plt
from collections import defaultdict
import pickle
from scipy.interpolate import CubicSpline,interp1d,UnivariateSpline,barycentric_interpolate
from copy import deepcopy
import scipy as sp
import scipy.stats
import itertools

def data_extraction(files,filename,keyword,rule):
    dir_dict = dict()
    for file in files:
        dir_split = file.split('/')
        op = open(file,'r')
        dataread = op.read().splitlines()
        op.close()
        dictionary_list = []
        for i in range(0,len(keyword)-1):
            dir1loc = dir_split.index([s for s in dir_split if keyword[i] in s][0])
            dir1=dir_split[dir1loc]
            listofvalues = re.findall("\d+\.\d+", dir1)
            converted = [float(i) for i in listofvalues]
            process_rule = rule[i](np.array(converted))
            dictionary_list.append(process_rule)
    
        dir2loc = dir_split.index([s for s in dir_split if keyword[-1] in s][0])
        dir2=dir_split[dir2loc]
        innerstore = float(re.findall("\d+\.\d+", dir2)[0])
        
        for i in range(0,len(dataread),1):
            newlist = []
            splitdata = dataread[i].split(' ')
            phase = splitdata[0][:-5]
            newlist+=dictionary_list
            newlist.append(phase)      
            if int(splitdata[2])==2 or int(splitdata[2])==0 or int(splitdata[2])==3:
                if tuple(newlist) in dir_dict:
                    temparray = np.array([innerstore,float(splitdata[1])]) 
                    dir_dict[tuple(newlist)] = np.vstack((dir_dict[tuple(newlist)],temparray))
                else:
                    dir_dict[tuple(newlist)] = np.array([innerstore,float(splitdata[1])])
            else:
                print('Simulations have STATUS 0 1 3...revist extractF0.dat...')
                break
    return dir_dict



















def obtain_unique_dims(dictionary):
    
    keys = list(dictionary)
    var = len(keys[0])
    returnlist = []
    for i in range(0,var,1):
        temp = []

        for key in keys:
            temp.append(key[i])
        returnlist.append(sorted(list(set(temp))))
    return returnlist

def Clean_Array(Array):
    unique_array = np.unique(Array[:,0:2],axis=0)
    # print(unique_array)
    shape = np.shape(unique_array)
    return_array = np.zeros((shape[0],3))
    for i in range(0,shape[0]):
        loc = np.where((unique_array[i,0]==Array[:,0])&(unique_array[i,1]==Array[:,1]))[0]
        if len(loc)==1:
            return_array[i,:] = Array[loc,:]
        if len(loc)>1:
            temp_compare = np.zeros(len(loc))
            count = 0
            for j in loc:
                temp_compare[count] = Array[j,2]
                count+=1
            choice = np.where(temp_compare==np.min(temp_compare))[0]
            if len(choice)>0:
                return_array[i,:] = Array[loc[choice[0]],:]
            else:
                return_array[i,:] = Array[loc[choice],:]
    return return_array

def Combine_Phases(dictionary,phase):  
    count = 0          
    header = list(dictionary)
    for key in header:
        if key.find(phase)!=-1:
#            print(key)
            if count == 0:
                temp_array = deepcopy(dictionary[key])
                count+=1
            else:
                temp_array = np.vstack((temp_array,dictionary[key]))
    return temp_array


# path = "/media/tquah/TOSHIBA EXT/Projects/chiN_60_asymdir"
# path = "/media/tquah/TOSHIBA EXT/Projects/sweep-asym-armlength_BCC_fix"

path = "/media/tquah/TOSHIBA EXT/Projects/sweep-asym-armlength_BCC_fix"

path = '/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_corrected'

export_path = "/home/tquah/Figures/11-11-2020-PhasePlotterTool/"

pwd = os.getcwd()
os.chdir(path)

directories = 'chiAB_0.0289*/Nsc*/fA0.*'
filename = 'F0_phases.dat'
keyword = ['Nsc','f']

list_of_directory = glob.glob(os.path.join(directories,filename), recursive = True)


tuplerule = lambda x: np.sqrt((x[1]+1)/(x[0]+1))
ABratio = lambda x : x*1
chiN = lambda x : 2100*x
fA = lambda x : x*1

rules = [tuplerule,fA]

data_dict = data_extraction(list_of_directory,filename,keyword,rules)

# phaselist = ['LAM','HEX','BCC','DIS','GYR']
phaselist = ['BCC','HEX','A15','SIGMA','DIS']

        
listofvar = obtain_unique_dims(data_dict)

header = list(data_dict)
#print(listofvar)

phase_dict = dict()
for key in header:
    shape = np.shape(data_dict[key])
    if len(shape)==2:
        temp_array = np.zeros((shape[0],3))
        temp_array[:,0] = data_dict[key][:,0]
        temp_array[:,1] = key[0]
        temp_array[:,2] = data_dict[key][:,1]

    else:
        temp_array = np.zeros((1,3))
        temp_array[:,0] = data_dict[key][0]
        temp_array[:,1] = key[0]
        temp_array[:,2] = data_dict[key][1]


    if key[1] in list(phase_dict):
        phase_dict[key[1]] = np.vstack((phase_dict[key[1]],temp_array))
    else:        
        phase_dict[key[1]] = temp_array


#only to clean up SIGMA shinanigans 
temp_array = Combine_Phases(phase_dict,'SIGMA')
clean_array = Clean_Array(temp_array)
header = list(phase_dict)
for key in header:
    if key.find('SIGMA')!=-1:
        phase_dict.pop(key, None)
phase_dict['SIGMA'] = clean_array



alldir = os.listdir()
if 'PHASE_FREE_ENERGY' not in alldir:
    os.mkdir('PHASE_FREE_ENERGY')
os.chdir('PHASE_FREE_ENERGY')

for phase in list(phase_dict):
    print(phase)
    np.savetxt(phase+'.dat',phase_dict[phase],header = '# fA epsilon free_energy',delimiter = ' ')




os.chdir(pwd)
