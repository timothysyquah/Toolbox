#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 19:11:33 2020

@author: tquah
"""

import os 
import glob
import numpy as np

file = 'LAM.out'
dir = os.getcwd()
dir_list = glob.glob(dir+'/L_*')
D_0 = []
L = []



for i in range(0,len(dir_list),1):
    fullpath = os.path.join(dir_list[i],file)
    segment = fullpath.split('/')
    if os.path.exists==False:
        print(file+' DOES NOT EXIST...SKIP')
        continue
    convergepath = os.path.join(dir_list[i],'STATUS')
    co = open(convergepath,'r')
    statusread=co.read()
    co.close()
    
    if int(statusread)!=2:
        print(segment[-2]+' NOT CONVERGED...SKIP')
        continue
    
    
    
    
    fo = open(fullpath,'r')
    content = fo.read().split('\n')
    fo.close()

    r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
    result = float(r[r.find("(")+1:r.find(")")])
    D_0.append(result)
    
    L_split = segment[-2].split('_')
    L.append(float(L_split[1]))
    
    
Lnp = np.array(L).T
Dnp = np.array(D_0).T

data = np.vstack((L,D_0)).T
export_path = os.path.join(dir,'LD_data.dat')
np.savetxt(export_path,data)
