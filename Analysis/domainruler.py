#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 19 10:29:17 2021

@author: tquah
"""

import os
import numpy as np
import glob
import re
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['axes.labelsize'] = 15
rcParams['xtick.labelsize'] = 12
rcParams['ytick.labelsize'] = 12
rcParams['legend.fontsize'] = 12

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
    operator_path = pathstatus+'operators.dat'
    store = np.min(np.loadtxt(operator_path)[:,1])


    if int(statusread)!=2:
        # print("Convergance Fail")
        return False
    else:
        fo = open(file,'r')
        content = fo.read().split('\n')
        fo.close()
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        result = float(r[r.find("(")+1:r.find(")")])
        
        return [result,store]



def fill_dictionaries(listoffile,desiredparameters,paramarray ,extractmethod,cutoffpt = 500):
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



IDIR = os.getcwd()

os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/')

desired_filename = 'LAM.out'

list_of_files = glob.glob('**/'+desired_filename,recursive=True)


important_directory,purged_list,paramarray,purged_file_list = determine_directory_structure(list_of_files)

data_dictionary = fill_dictionaries(purged_file_list,important_directory,paramarray,get_domain_spacing_SCFT_simple)

plt.close('all')
#
fig = plt.figure(figsize=(10,10))

for header in list(data_dictionary):
    if header[-1]==20:
        if float(header[1])>0.5:
            print(header[1])

            ax = plt.gca()
        
            ax.scatter(data_dictionary[header][:,0],data_dictionary[header][:,1],marker = 'o',label = f'$a = {header[0]}$ & $\zeta = {1/header[1]}$')
            ax.set_yscale('log')
            ax.set_xscale('log')
            params,pcov = curve_fit(powlaw,data_dictionary[header][-2:,0],data_dictionary[header][-2:,1],p0 = [1,2/3])
            print(params[1])
plt.xlabel('$N_{bb}$')
plt.ylabel('$D_0 \sqrt{6}/b_{ref}$')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')
plt.tight_layout()

plt.savefig('domainspacing.pdf',dpi = 300)

fig = plt.figure()
for header in list(data_dictionary):
    if header[-1]==10:
        if float(header[1])>0.5:
            print(header[1])

            ax = plt.gca()
            intersection,i1,i2 = np.intersect1d(data_dictionary[header][:,0], data_dictionary[(0,0,20.0)][:,0],return_indices=True)
            ax.scatter(data_dictionary[header][i1,0],data_dictionary[header][i1,2]-data_dictionary[(0,0,20.0)][i2,2],marker = 'o',label = f'$a = {header[0]}$ & $\zeta N = {header[1]}$')
plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

plt.xlabel('$N_{bb}$')
plt.ylabel('$F-F_{a=0,\zeta N = 0}$')
plt.tight_layout()
plt.savefig('free_energy.pdf',dpi = 300)


# fig = plt.figure()
# for header in list(data_dictionary):
#     if header[-1]==20.0:
#         print(header)
#         ax = plt.gca()
#         intersection,i1,i2 = np.intersect1d(data_dictionary[header][:,0], data_dictionary[(header[0],0,20.0)][:,0],return_indices=True)
#         ax.scatter(data_dictionary[header][i1,0],data_dictionary[header][i1,2]-data_dictionary[(header[0],0,20.0)][i2,2],marker = 'o',label = f'$a = {header[0]}$ & $\zeta N = {header[1]}$')
# plt.legend(bbox_to_anchor=(1.05, 1), loc='upper left')

# plt.xlabel('$N_{bb}$')
# plt.ylabel('$F-F_{\zeta N = 0}$')
# plt.tight_layout()
# plt.savefig('free_energy_1.pdf',dpi = 300)
