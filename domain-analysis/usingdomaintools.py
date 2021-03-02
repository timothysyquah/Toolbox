#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Feb 28 15:59:16 2021

@author: tquah
"""

import numpy as np
import regex as re
import os

IDIR = os.getcwd()
datapath = '/media/tquah/TOSHIBA EXT/Projects/DMREF/sweep-asym-armlength_corrected'

path = '/home/tquah/toolbox_github/Daniel_Phase_Diagram_Tool/Implimentation_Dataset/under.dat'
op = open(path,'r')
startlist = op.read().split('\n')
op.close()


#first we need to only get what we want
def return_likely_point(pt,data):
    r = np.sum((data[:,0:2]-pt)**2,axis=1)
    loc = np.where(r==np.min(r))[0]
    return data[loc],loc




phases = ['BCC','SIGMA','A15']
desired_fA  =[[[0.14,1.0],[0.14,1.5],[0.13,1.9],[0.13,2.2],[0.13,2.7]],\
              [[0.19,1.5],[0.21,1.9],[0.20,2.2],[0.20,2.7]],\
              [[0.25,1.5],[0.26,1.9],[0.26,2.2],[0.26,2.7]]
              ]
    
dirlist =[ele for ele in startlist if ele.strip()] 


#only choose 2
key_locations = [3,[1,2]]
#rules
tuplerule = lambda x: np.sqrt((x[1]+1)/(x[0]+1))
fA = lambda x : x*1

rules = [fA,tuplerule]




#not going to make it too general take 2nd column and 3rd column
def make_grid(dirlist,rules,key_locations):
    xlist = []
    ylist = []
    for i in range(len(dirlist)):
        listofstr = re.findall("\d+\.\d+", dirlist[i])
        listofvalue = np.array([float(i) for i in listofstr])
        x = rules[0](listofvalue[key_locations[0]])
        y = rules[1](listofvalue[key_locations[1]])
        xlist.append(x)
        ylist.append(y)
    return np.vstack((xlist,ylist)).transpose()
    
    
datalist = make_grid(dirlist,rules,key_locations)
for i in range(len(phases)):
    newsubmit_directory_list = []
    os.chdir(datapath)

    for j in range(len(desired_fA[i])):
        datapt,location = return_likely_point(desired_fA[i][j],datalist)
        
        op = open(dirlist[location[0]].strip())
        phasedatatemp = op.read().split('\n')
        phasedata =[ele for ele in phasedatatemp if ele.strip()] 

        op.close()
        FE = []
        Cname = []
        for k in range(len(phasedata)):
            component = phasedata[k].split(' ')
            if phases[i] in component[0]:
                Cname.append(component[0])
                FE.append(float(component[1]))
    
        if len(Cname)==1:
            newsubmit_directory_list.append(dirlist[location[0]][:-15]+Cname[0])
        else:
            FEloc = np.where(np.min(FE)==FE)[0][0]
            newsubmit_directory_list.append(dirlist[location[0]][:-15]+Cname[FEloc])
    
        
    os.chdir(IDIR)
    op = open(phases[i]+'.imp','w+')
    for j in range(len(newsubmit_directory_list)):
        op.write(newsubmit_directory_list[j]+'\n')
    op.close()

