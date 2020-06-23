#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 19:11:33 2020

@author: tquah
"""
import os
import glob
import datetime
import shutil
import numpy as np
from collections import Counter


filetype = 'LAM.out'
EXPORT_PATH = '/home/tquah/toolbox/SandBox/extractlengthdomain/EXPORT_LOCAL'
export_base_name = 'LD.dat'
density_file = 'density.dat'
status_file = 'STATUS'
Today  = str(datetime.date.today())
today_dir = os.path.join (EXPORT_PATH,Today)

if os.path.exists(today_dir)==False:
    print('making today directory')
    os.mkdir(today_dir)

IDIR = os.getcwd()
IDIR_pieces = IDIR.split('/')
path_extraction = os.path.join(IDIR,'**/'+filetype)
save_other_list = []
tol = 1e-3
for file in glob.iglob(path_extraction,recursive=True):
    file_pieces = file.split('/')
    overall_list = file_pieces+IDIR_pieces
    counts = Counter(overall_list)
    category_list = [k for k in overall_list if counts[k] == 1]
    save_other_list.append(category_list[:-2])
    
b_set = list(set(map(tuple,save_other_list)))  #need to convert the inner lists to tuples so they are hashable



for i in range(0,len(b_set),1):
    for j in range(0,len(b_set[i]),1):
        if j==0:
            dir_path = os.path.join(IDIR,b_set[i][j])    
            name = b_set[i][j]+'_'+export_base_name
        else:
            dir_path = os.path.join(dir_path,b_set[i][j])
            name =  b_set[i][j]+'_'+name
        path_L = os.path.join(dir_path,'L_**/'+filetype)
    L = []
    D = []
    for file in glob.iglob(path_L,recursive=True):
        file_pieces = file.split('/')
        overall_list = file_pieces+IDIR_pieces
        counts = Counter(overall_list)
        category_list = [k for k in overall_list if counts[k] == 1]
        WDIR = os.path.join(dir_path,category_list[-2])
        # print(WDIR)
        
        density_path = os.path.join(WDIR,density_file)
        status_path = os.path.join(WDIR,status_file)
        co = open(status_path,'r')
        statusread=int(co.read())
        co.close()
        
    

        if int(statusread)!=2:
            print("Convergance Fail")
            continue
        if statusread==2:
            density_data = np.loadtxt(density_path)
            amp_check = abs(np.min(density_data[:,1])-np.max(density_data[:,1]))
            if amp_check<tol:
                print("Disorder Fail")
                continue

        fo = open(file,'r')
        content = fo.read().split('\n')
        fo.close()

        L.append(float(category_list[-2][2:]))

        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        result = float(r[r.find("(")+1:r.find(")")])
        D.append(result)
    print(name)
    data = np.vstack((L,D)).T
    
    export_path = os.path.join(today_dir,name)
    np.savetxt(export_path,data)