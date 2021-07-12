#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 12 15:04:54 2021

@author: tquah
"""

import numpy as np
from data_to_pif_json import *
import glob
import shutil
import os
import re


def check_status(path):
    status_path = os.path.join(path,'STATUS')
    return np.loadtxt(status_path)
    
##################################################################################
# Parameter I need to change
##################################################################################
#Directories I need to run
path = '/home/tquah/Photonic_Band_Gap_Dataset/DGC'

# path = '/home/tquah/Photonic_Band_Gap_Dataset/FJC'

chainstatistic = 'discrete_gaussian'

# chainstatistic = 'freely_jointed'


DOI = '10.1021/acsmacrolett.0c00380'
SOFTWARE = 'PolyFTS'
METHOD = 'SCFT'
NAME = 'Timothy Quah'
EMAIL = 'timothy_quah@ucsb.edu'
ExportDirectory = 'PIF_JSON'
Overwrite = True
#################################################################################
#################################################################################

IDIR = os.getcwd()
os.chdir(path)

if os.path.exists(ExportDirectory):
    if Overwrite:
        print('Overwriting...')
        shutil.rmtree(ExportDirectory)
        os.mkdir(ExportDirectory)
    else:
        print('Directory exists-since overwrite is true only updating directory')
else:
    os.mkdir(ExportDirectory)





list_of_dir = glob.glob('./**/**/*Phase/')
for i in range(len(list_of_dir)):
    local_path = list_of_dir[i]
    split_name = local_path.split('/')
    Nsc = float(re.findall('[\d]*[.][\d]+',split_name[1])[0])
    fA = float(re.findall('[\d]*[.][\d]+',split_name[2])[0])
    Phase = split_name[3][:-5]
    outname = f'{chainstatistic}-Nsc_{Nsc}-fA_{fA}-{Phase}Phase'
    export_path = os.path.join(ExportDirectory,outname)
    infile_path = os.path.join(local_path,f'{Phase}.out')
    if os.path.exists(export_path):
        if Overwrite==False:
            continue
        else:
            pass
    
    
    if check_status(local_path)==2:
        print(infile_path)

        LAZY_maker(infile_path,export_path,DOI,SOFTWARE,METHOD,NAME,EMAIL)
        
        
    else:
        print(f'Skipping...returned STATUS {check_status(local_path)}')
    
    
    




os.chdir(IDIR)