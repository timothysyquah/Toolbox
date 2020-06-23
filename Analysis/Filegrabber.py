#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 14:18:31 2020

@author: tquah
"""

import os
import glob
import datetime
import shutil

n_cat = -3

EXPORT_PATH = '/home/tquah/EXPORTLOCAL/'
Today  = str(datetime.date.today())
today_dir = os.path.join (EXPORT_PATH,Today)

if os.path.exists(today_dir)==False:
    print('making today directory')
    os.mkdir(today_dir)
    
filetype = 'test.dat'



file_dir = os.path.join (today_dir,filetype)
if os.path.exists(file_dir)==False:
    os.mkdir(file_dir)
    print(f'making {filetype} directory')

IDIR = os.getcwd()
path_extraction = os.path.join(IDIR,'**/'+filetype)

for file in glob.iglob(path_extraction,recursive=True):
    pieces = file.split('/')
    category = pieces[n_cat:-1]
    name = ''
    for i in range(0,len(category),1):
        name+=category[i]+'_'
    name+=pieces[-1]
    export_file_path = os.path.join(file_dir,name)
    shutil.copyfile(file,export_file_path)
