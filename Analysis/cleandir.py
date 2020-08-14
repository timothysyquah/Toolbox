#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 13 23:50:41 2020

@author: tquah
"""
import os
import shutil
import re
path = os.getcwd()
# path = '/home/tquah/IMPORT_BRAID/NSCASYM_02_other'
chilist = os.listdir(path)
for chi in chilist:
    if chi.find('chi')!=-1:
        path2 = os.path.join(path,chi)
        # os.chdir(path2)
        Nsclist = os.listdir(path2)
        for nsc in Nsclist:
            if nsc.find('Nsc')!=-1:
                path3 = os.path.join(path2,nsc)
                flist = os.listdir(path3)
                for f in flist:
                    if f.find('f')!=-1:
                        path4 = os.path.join(path3,f)
                
                        fval = float(re.findall("\d+\.\d+", f)[0])
                        if fval>0.5:
                            print('Delete %0.3f'%fval)
                            shutil.rmtree(path4)

        
        


