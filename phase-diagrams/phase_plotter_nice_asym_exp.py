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
plt.rc('font', family='sans-serif')
plt.rc('text', usetex=True)

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


file_path_full = 'asym-10-23-2020.dat'

fo = open(file_path_full, 'r')
data_from_file = fo.read().splitlines()
fo.close()

Neff = 1
def func(x, a, b, c):
    return a * np.exp(-b * x) + c

NscA = np.arange(2,18+1e-6,2)
NscB = 40-NscA
ratio = np.around(np.sqrt((NscB+1)/(NscA+1)),3)
ratioreplace = NscB/NscA

data_from_file = whitespace_remover(data_from_file)
phase_list, phase_loc = phase_list_parser(data_from_file)
data_dict = extract_data(phase_list, data_from_file, phase_loc, delimiter=' ')
from scipy.optimize import curve_fit
fig = plt.figure(figsize=(6.5,5))

ax = plt.subplot(111)
for boundary in data_dict:
    x = data_dict[boundary][:,0]
    y = data_dict[boundary][:,1]*Neff
    # count=0
    # for r in ratio:
    #     loc = np.where(y==r)[0]
    #     if len(loc)>0:
            
    #         y[loc] = ratioreplace[count]
    #     count+=1
    
    
    # popt, pcov = curve_fit(func, x, y)
    # xplot = np.linspace(np.min(x),np.max(x),100)
    # yplot = func(xplot,popt[0],popt[1],popt[2])
    ax.plot(x,y,'k')
    # plt.plot(x,y,'-ok',marker='s',markersize=5)

# path1_exp = '/home/tquah/Projects/PhaseDiagram/lam.csv'
# path2_exp = '/home/tquah/Projects/PhaseDiagram/cylinder.csv'
# path3_exp = '/home/tquah/Projects/PhaseDiagram/dis.csv'
# path4_exp = '/home/tquah/Projects/PhaseDiagram/spherical.csv'

path3_dob = '/media/tquah/TOSHIBA EXT/Projects/PhaseDiagram/Dobreninpts.dat'
# lamarray = np.loadtxt(path1_exp,delimiter=',')
# cylinderarray = np.loadtxt(path2_exp,delimiter=',')
# disarray = np.loadtxt(path3_exp,delimiter=',')
# sphericalarray = np.loadtxt(path4_exp,delimiter=',')
path_exp = '/media/tquah/TOSHIBA EXT/Projects/PhaseDiagram/Experimental_Dataset_v2.csv'

df = pd.read_csv(path_exp)
header = list(df)

references = sorted(list(set(list(df['Reference']))))
update_references = []
for ref in references:
    newref = ref.replace('-',' ')
    update_references.append(newref)
               
fullreference = list(df['Reference'])
fullphases =  list(df['Phase'])
                  
                  
phases = sorted(list(set(list(df['Phase']))))
symbol = ['o','+','+','h']
color = ['k','r','g','b']
ab_ratio = np.sqrt(np.array(df['NscB']/df['NscA']))*np.array((df['bB']/df['bA']))

NscAverage =  np.array((df['NscB']+df['NscA'])/2)

Nsc_Backbone = NscAverage/np.array(df['Nbb(Total)'])

NaNloc = np.where(np.isnan(Nsc_Backbone)==True)[0]
Nsc_Backbone[NaNloc]=0


count1 = 0
x = np.arange(10)
condlist = [x<3, x>5]
choicelist = [x, x**2]
# print(np.select(condlist, choicelist))





for ref in references:
    count2 = 0 
    for phase in phases:
        refindex= [i for i, x in enumerate(fullreference) if x ==ref]
        phaseindex= [i for i, x in enumerate(fullphases) if x ==phase]
        
        
        loc = list(set(refindex) & set(phaseindex))
        if len(loc)>0:
            
            # print(Nsc_Backbone[loc])
            
            

            mloc = np.where(Nsc_Backbone[loc]<0.8)[0]
            floc = []
            for m in mloc:
                floc.append(np.where(Nsc_Backbone[loc][m]==Nsc_Backbone)[0][0])
            if len(floc)>0:
                
                yraw = ab_ratio[floc]
                # print(yraw)
                xraw = np.array(df['fA'][floc])
    
                yloc = np.where(yraw>1)[0]
                xfilter = xraw[yloc]
                yfilter = yraw[yloc]
                
                print(ref)
                print(floc)

                
                
                # print(yfilter)
                # nscnbbloc = np.where(Nsc_Backbone[loc]<1.0)[0]
                # print(nscnbbloc)
                yloc = np.where(1<y)[0]
            
                # print(xfilter)
            
                name = phase+'-'+ref
                if len(xfilter)>0:
                    ax.scatter(xfilter,np.sqrt(yfilter),c=color[count1],marker = symbol[count2],label = name)      
        count2+=1
    count1+=1

dobarray = np.loadtxt(path3_dob)
ypos = 6
# plt.xlim(-.05,1.05)
# plt.xlim(0,0.5)

# plt.ylim(1.22,4)
# textsize = 15
# plt.text(0.40,1.3,'LAM',color = 'k',size = textsize)
# plt.text(0.3,2,'GYR',color = 'k',size = textsize)
# plt.arrow(0.37,2.05,.025,0,color = 'k',head_length=0.01,width = 0.05)
# plt.text(0.17,3.5,'A15',color = 'k',size = textsize)

# plt.text(0.28,3,'HEX',color = 'k',size = textsize)
# plt.text(0.12,2,'BCC',color = 'k',size = textsize)
# plt.text(0.025,1.5,'DIS',color = 'k',size = textsize)

# plt.xlabel(r'$f_A$')
# plt.ylabel(r'$N_{sc,B}/N_{sc,A}$')

# plt.plot(lamarray[:,0],lamarray[:,1],'+',label = 'LAM')
# plt.plot(cylinderarray[:,0],cylinderarray[:,1],'o',label = 'Cylinders')
# plt.plot(disarray[0],disarray[1],'<',label = 'DIS')
# plt.plot(sphericalarray[0],sphericalarray[1],'^',label = 'Spherical')

# plt.legend()
alpha = 0.5
# y = np.sqrt(np.power(dobarray[:,1],4))
# plt.plot(1-dobarray[:35,0],y[:35],'--r',alpha = 0.5)
# plt.plot(1-dobarray[36:77,0],y[36:77],'--r',alpha = 0.5)
# plt.plot(1-dobarray[78:134,0],y[78:134],'--r',alpha = 0.5)
# plt.plot(1-dobarray[135:,0],y[135:],'--r',alpha = 0.5)

# pos1 = (np.max(1-dobarray[:35,0])+0)/2


# plt.text(0.1,3,'S',color = 'r',size = textsize)
# plt.text(0.25,3,'C',color = 'r',size = textsize)
# plt.text(0.48,3,'L',color = 'r',size = textsize)

# plt.text(0.2,8,'S',color = 'r',size = textsize)
# plt.text(0.37,8,'C',color = 'r',size = textsize)
# plt.text(0.48,8,'L',color = 'r',size = textsize)
# plt.text(0.9,8,'C',color = 'r',size = textsize)
# plt.text(1.0,88,'S',color = 'r',size = textsize)



# plt.text(0.1,18,'S',color = 'r',size = textsize)
# plt.text(0.37,18,'C',color = 'r',size = textsize)
# plt.text(0.6,18,'L',color = 'r',size = textsize)
# plt.text(0.9,18,'C',color = 'r',size = textsize)
# plt.text(1.0,18,'S',color = 'r',size = textsize)

# plt.text(0.2,12,'S',color = 'k',size = textsize)
# plt.text(0.37,12,'C',color = 'k',size = textsize)
# plt.text(0.6,12,'L',color = 'k',size = textsize)
# plt.text(0.82,12,'C',color = 'k',size = textsize)
# plt.text(0.9,12,'S',color = 'k',size = textsize)





ypos = 6
plt.xlim(-0.05,1.05)
plt.ylim(0.9,3.7)
# textsize = 15
# plt.text(0.40,1.5,'LAM',color = 'k',size = textsize)
# plt.text(0.3,2,'GYR',color = 'k',size = textsize)
# plt.arrow(0.37,2.05,.025,0,color = 'k',head_length=0.01,width = 0.05)
# plt.text(0.185,3.5,'A15',color = 'k',size = textsize)

# plt.text(0.28,3,'HEX',color = 'k',size = textsize)
# plt.text(0.135,2,'BCC',color = 'k',size = textsize)
# plt.text(0.05,1.5,'DIS',color = 'k',size = textsize)
ax.legend(bbox_to_anchor=(1.1, 1.05))
plt.tight_layout()

plt.xlabel(r'$f_A$')
plt.ylabel(r'$\epsilon$')
plt.tight_layout()
# plt.savefig('/home/tquah/Presentations/FirstYearTalkQuah/images/phasediagramasym_1.png',dpi=300)
plt.savefig('/home/tquah/Figures/10-23-2020-GFKD-Meeting/asymexp_phasediagram.png',dpi=300)

# plt.figure()
# plt.plot(1-dobarray[:35,0],y[:35],'--r')
# plt.plot(1-dobarray[36:77,0],y[36:77],'--r')
# plt.plot(1-dobarray[78:134,0],y[78:134],'--r')
# plt.plot(1-dobarray[135:,0],y[135:],'--r')
# plt.ylim(0,19)
