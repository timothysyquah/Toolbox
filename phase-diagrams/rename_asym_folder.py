#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug  6 15:44:57 2020

@author: tquah
"""
import os
import re


Import_Path ="/home/tquah/IMPORT_BRAID/NSCASYM/chiAB_0.0289/"
newname = 'ABratio_'

list_of_directories = os.listdir(Import_Path)

for direct in list_of_directories:
    list_values = re.findall("\d+\.\d+", direct)
    print(list_values)
    ratio = float(list_values[1])/float(list_values[0])
    newname_value = newname+('%0.3f'%ratio)
    print(newname_value)
    oldpath = os.path.join(Import_Path,direct)
    newpath = os.path.join(Import_Path,newname_value)
    os.rename(oldpath,newpath)