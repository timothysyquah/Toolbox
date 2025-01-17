#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Aug  7 14:50:24 2020

@author: tquah
"""


import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
import pandas as pd
from scipy.signal import savgol_filter
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'

rcParams['axes.labelsize'] =20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 11

def file_sorter(file_list, filetype):
    list_files = []
    for file in file_list:
        if file[-4:] == filetype:
            list_files.append(file)
    return sorted(list(set(list_files)))

    
def phase_list_parser(data):
    phase_list = []
    phase_loc = []
    count = 0
    for lines in data:
        if len(lines) > 0:
            if lines[0] == "#":
                if len(lines[1:]) <= 12:
                    phase_list.append(lines[1:])
                    phase_loc.append(count)
        count += 1

    return list(phase_list), phase_loc


def whitespace_remover(list_file):
    while ("" in list_file):
        list_file.remove("")
    return list_file


def extract_data(phase_list, data_from_file, phase_loc, delimiter=' '):
    data_dict = dict()
    for i in range(0, len(phase_list), 1):
        if i == len(phase_list) - 1:
            data_line_list = data_from_file[phase_loc[i] + 1:]
        else:
            data_line_list = data_from_file[phase_loc[i] + 1:phase_loc[i + 1]]
        data_array = np.zeros((len(data_line_list), 2))

        for j in range(0, len(data_line_list), 1):
            data_split = data_line_list[j].split(delimiter)
            data_array[j, 0] = np.float(data_split[0])
            data_array[j, 1] = np.float(data_split[1])

        data_dict[phase_list[i]] = data_array
    return data_dict

def dictionary_depth(d):
    if isinstance(d, dict):
        return 1+(max(map(dictionary_depth, d.values())) if d else 0)
    return 0

def return_full_list(d):
    full_phase_list = []
    full_fA_list = []
    full_chi_list = []
    full_chi_list = list(d)
    for chi in full_chi_list:
        partialphase = list(d[chi])
        full_phase_list+=list(d[chi])
        for phase in partialphase:
            full_fA_list+=list(d[chi][phase][:,0])
    full_phase_list = sorted(list(set(full_phase_list)))
    full_fA_list = sorted(list(set(full_fA_list)))
    return full_chi_list,full_phase_list,full_fA_list

def convert_dataframe(data_dict,full_chi_list,full_phase_list,full_fA_list):
    chi_len = len(full_chi_list)
    phase_len = len(full_phase_list)
    fa_len = len(full_fA_list)
    ALL_Data_Array = np.zeros((fa_len,(chi_len*phase_len)+1))
    ALL_Data_Array[:,0] = np.array(full_fA_list)
    header = []
    for i in range(0,chi_len,1):
    
        local_phase_list = sorted(list(data_dict[full_chi_list[i]]))
    
        for j in range(0,len(local_phase_list),1):
            phase_loc = full_phase_list.index(local_phase_list[j])
    
            for k in range(0,len(full_fA_list),1):
                # print(phase_len*i+j+1)
                fA_logic = np.where(full_fA_list[k]==\
                                    data_dict[full_chi_list[i]][full_phase_list[phase_loc]][:,0])[0]
                fA_length = len(fA_logic)
                
                if fA_length>1:
                    print('Error Issues')
                    continue
                elif fA_length==1:
                    ALL_Data_Array[k,(phase_len*i)+j+1] =\
                        data_dict[full_chi_list[i]][full_phase_list[phase_loc]][fA_logic,1]
    header = []
    header.append('F0')
    for i in range(0,chi_len,1):
        for j in range(0,phase_len,1):
            phase_loc = full_phase_list.index(local_phase_list[j])
            chiN = np.round(full_chi_list[i]*2000,3)
            header.append(f'chiN{chiN: 0.3f} {full_phase_list[phase_loc]}Phase')
    
    df = pd.DataFrame(ALL_Data_Array,columns = header)
    return df
def extend_data_cubicspline(df,dfA = 0.001):

    header = list(df)
    
    fA_array = np.array(df[header[0]])
    
    fA_interp = np.arange(np.min(fA_array),np.max(fA_array),dfA)
    
    large_data_array = np.zeros((len(fA_interp),len(header)))
    large_data_array[:,0] = fA_interp
    for i in range(1,len(header),1):
        nonzero_loc = np.nonzero(np.array(df[header[i]]))
        fA_array = np.array(df[header[0]])[nonzero_loc]
        if len(fA_array)>2:
            data_array = np.array(df[header[i]])[nonzero_loc]
            lower_bound = np.min(fA_array)
            upper_bound = np.max(fA_array)
            cs = CubicSpline(fA_array,data_array)
            eval_loc = np.where((fA_interp>lower_bound)&(fA_interp<upper_bound))
            large_data_array[eval_loc,i] = cs(fA_interp[eval_loc])
    
            
        if len(fA_array)<2 and len(fA_array)>0:
            data_array = np.array(df[header[i]])[nonzero_loc]
            large_data_array[nonzero_loc,i] = data_array
    new_df = pd.DataFrame(large_data_array,columns = header)
    return new_df


#file_path_full = 'asym-10-23-2020.dat'
file_path_full = '/home/tquah/BottlebrushPaper/PhaseBoundaries/asymbottlebrush.dat'
fo = open(file_path_full, 'r')
data_from_file = fo.read().splitlines()
fo.close()

def func(x, a, b, c):
    return a * np.exp(-b * x) + c

data_from_file = whitespace_remover(data_from_file)
phase_list, phase_loc = phase_list_parser(data_from_file)
data_dict = extract_data(phase_list, data_from_file, phase_loc, delimiter=' ')
from scipy.optimize import curve_fit
plt.figure(figsize=(5,5))
for boundary in data_dict:
    x = data_dict[boundary][:,0]
    y = data_dict[boundary][:,1]
    
    ysort = np.argsort(y)
    # popt, pcov = curve_fit(func, x, y)
    # xplot = np.linspace(np.min(x),np.max(x),100)
    # yplot = func(xplot,popt[0],popt[1],popt[2])
    plt.plot(x[ysort],y[ysort],'k')
    # plt.plot(x,y,'-ok',marker='s',markersize=5)

ypos = 6
plt.xlim(.0,1.0)
# plt.ylim(1.,19)
textsize = 13
plt.text(0.65,2.5,'L',color = 'k',size = textsize,ha = 'center')
plt.text(0.5,2.5,'G',color = 'k',size = textsize,ha = 'center')
plt.text(0.71,2.5,'G',color = 'k',size = textsize,ha = 'center')
plt.text(0.83,2.5,'C',color = 'k',size = textsize,ha = 'center')
plt.text(0.83,2.0,'S',color = 'k',size = textsize,ha = 'center')

# plt.arrow(0.37,2.05,.025,0,color = 'k',head_length=0.01,width = 0.05)
# plt.text(0.2,3.8,'A15',color = 'k',size = textsize)
plt.text(0.295,3.5,'A15',color = 'k',size = textsize,ha = 'center')
plt.text(0.19,2.5,'$\sigma$',color = 'k',size = textsize,ha = 'center')

plt.text(0.38,2.5,'C',color = 'k',size = textsize,ha = 'center')
plt.text(0.155,1.2,'S',color = 'k',size = textsize,ha = 'center')
plt.text(0.07,2.5,'D',color = 'k',size = textsize,ha = 'center')
plt.text(0.95,2.5,'D',color = 'k',size = textsize,ha = 'center')

plt.xlabel(r'$f_A$')
plt.ylabel(r'$\epsilon$')

plt.tight_layout()
plt.savefig('/home/tquah/Figures/asym_phasediagram.png',dpi=300)
