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
rcParams['axes.labelsize'] = 18
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 15

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



IDIR = os.getcwd()

os.chdir('/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/main_test/')

desired_filename = 'LAM.out'

list_of_files = glob.glob('**/'+desired_filename,recursive=True)


important_directory,purged_list,paramarray,purged_file_list = determine_directory_structure(list_of_files)

data_dictionary = fill_dictionaries(purged_file_list,important_directory,paramarray,get_domain_spacing_SCFT_simple)

plt.close('all')
#
color = ['k','r','b','g']
marker = ['^','o','+','h']

fig = plt.figure(figsize=(6,6))
ax = plt.gca()
arraylist = []
for i in range(len(data_dictionary)):
    array = data_dictionary[list(data_dictionary)[i]]
    plt.scatter(array[:,0],array[:,1],c= color[i],marker= marker[i],label = f'$N_{{sc}} = {list(data_dictionary)[i][0]}$')
    params,pcov = curve_fit(powlaw,array[-2:,0],array[-2:,1])
    print(params)
    arraylist.append(np.hstack((array,np.ones((len(array),1))*list(data_dictionary)[i][0].transpose())))
plt.legend()
ax.set_xscale('log')
ax.set_yscale('log')

ax.set_xticks([40,50,60,70,80,90, 100])
ax.set_yticks([40,50,60,70,80,90, 100,110,120,130])


plt.xlabel('$N_{bb}$')
plt.ylabel('$D_0/l$')
plt.tight_layout()
plt.savefig('/home/tquah/Figures/largedomain.pdf',dpi=300)
vstack_array = np.vstack(arraylist)
f = open('domain.dat','w')
f.write('Nbb Do0 Nsc')
f.write('\n')
np.savetxt(f,vstack_array)
f.close()


