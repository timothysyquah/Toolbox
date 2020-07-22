#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:44:05 2020

@author: tquah
"""
import os
import numpy as np
dir_path = "/home/tquah/IMPORT_BRAID/fA0.30000/"
def domain_size_extractor(content):
    list_of_cell = list(filter(lambda x: 'a_1' in x, content))
    if len(list_of_cell)>2:
        print('High Dim')
        r = list_of_cell[0]
        Box_initial = abs(float(r[r.find("(")+1:r.find(",")]))
        r = list_of_cell[1]
        Box_change = abs(float(r[r.find("(")+1:r.find(",")]))
        conversion = Box_initial/Box_change
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        Box_Final = np.around(abs(float(r[r.find("(")+1:r.find(",")]))*conversion,5)
    else:
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        Box_Final = np.around(float(r[r.find("(")+1:r.find(")")]),5)
    return Box_Final


def volume_fraction(content):
    r = list(filter(lambda x: 'volfrac' in x, content))[0]
    return np.around(float(r[r.find("=")+1:r.find(";")]),5)

def stress_tensor_check(content):
    r = list(filter(lambda x: 'L2-norm of stress tensor' in x, content))[-1]
    return float(r[r.find("+")+1:])



dirlist = os.listdir(dir_path)
Phaselist = [x for x in dirlist if "Phase"  in x]

os.chdir(dir_path)
IDIR = os.getcwd()

for phase in Phaselist:
    newpath = os.path.join(IDIR,phase)
    os.chdir(newpath)
    Phase = phase[:-5] 
    
    so = open('STATUS','r')
    status = int(so.read())
    so.close()
    if status==2:
        print('Simulation Converged!')
        
        phaseout = open(f'{Phase}.out')
        content = phaseout.read().splitlines()
        phaseout.close()
        
        finalcell = domain_size_extractor(content)
        fA = volume_fraction(content)
        stress_tensor = stress_tensor_check(content)
        if stress_tensor<1e-10:
            print(phase)
            print(stress_tensor)
            print('DIS')
    

        # if Phase=='DIS':
        #     break
    else:
        print('Simulation Not Converged...Quitting!')
    
    
    os.chdir(IDIR)    