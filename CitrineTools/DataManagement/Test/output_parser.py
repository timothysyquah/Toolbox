#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  7 14:34:42 2021

@author: tquah
"""

import numpy as np
import os
import re
import math as m
import time


os.chdir('/home/tquah/Projects/TESTPHASES')
#there are duplicates of 3D only true one is labeled by phase
listoffiles = ['AGYRPhase/AGYR.out','ACYLPhase/ACYL.out','DISPhase/DIS.out','LAMPhase/LAM.out']
def nCr(n,r):
    return int(m.factorial(n)/m.factorial(r)/(m.factorial(n-r)))

def OutputParser(text):
    DataDictionary = dict()
    for line in text:
        if 'initialized in' in line:
            d = int(re.findall(r'\d+', line)[0])
            break
    arch_list = []

    for i in range(len(text)):
        spacegroup_name = 'None'
        if 'Space group name' in text[i]: 
            spacegroup_name = text[i][text[i].index('=')+1:].strip()
        
        if 'Real-space collocation mesh' in text[i]: 
            npwarray = np.zeros(d)
            npw = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i])
            for k in range(d):
                npwarray[k] = int(npw[k])

        if "Number of monomer species" in text[i]:
            specs = int(re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i])[0])

        if 'Chi parameters (original):' in text[i]: 
            combo = nCr(specs,2)
            chimatrix = np.zeros([combo,combo])
            for j in range(combo):
                values = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+j+1])
                for k in range(combo):
                    chimatrix[j,k] = float(values[k])
                    
        
        if 'Final simulation cell' in text[i]:
            celltensor = np.zeros([d,d])
            for j in range(d):
                values = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+j])
                for k in range(d):
                    a = float(values[2*k+1])
                    b = float(values[2*k+2])
                    finval = a*10**b
                    celltensor[j,k] = finval
            hamval = float(re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i-2])[0])
        
        
        
        if 'Discrete Backbone' in text[i]:
            stat = text[i+1]
            loc = stat.index('+')
            stat = stat[loc+1:].strip()
            
            speclist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+4])
            specarray = np.zeros(len(speclist))
            for j in range(len(speclist)):
                specarray[j] = float(speclist[j])
            
            blist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+5])
            barray = np.zeros(len(blist))
            for j in range(len(blist)):
                barray[j] = float(blist[j])
            beadlist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+6])
            beadarray = np.zeros(len(beadlist))
            for j in range(len(beadlist)):
                beadarray[j] = float(beadlist[j])

            mainback = [stat,specarray,barray,beadarray]
            arch_list.append(mainback)
            
            
            
        if 'Side arm index' in text[i]:
            stat = text[i+1]
            loc = stat.index('+')
            stat = stat[loc+1:].strip()
            
            speclist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+4])
            specarray = np.zeros(len(speclist))
            for j in range(len(speclist)):
                specarray[j] = float(speclist[j])
            
            blist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+5])
            barray = np.zeros(len(blist))
            for j in range(len(blist)):
                barray[j] = float(blist[j])
            beadlist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+6])
            beadarray = np.zeros(len(beadlist))
            for j in range(len(beadlist)):
                beadarray[j] = float(beadlist[j])
            graftlist = re.findall('-?\ *[0-9]+\.?[0-9]*(?:[Ee]\ *-?\ *[0-9]+)?', text[i+9])
            graftarray = np.zeros(len(graftlist))
            for j in range(len(graftlist)):
                graftarray[j] = float(graftlist[j])

            mainback = [stat,specarray,barray,beadarray,graftarray]
            arch_list.append(mainback)

    DataDictionary['Space_Group'] = spacegroup_name
    DataDictionary['NPW'] = npwarray
    DataDictionary['CHI'] = chimatrix
    DataDictionary['ARCH'] = arch_list
    DataDictionary['CELL'] = celltensor
    DataDictionary['FREE_ENERGY'] = hamval
    return DataDictionary






from pypif import pif
from pypif.obj import *


# chemical_system = ChemicalSystem()
# chemical_system.chemical_formula = 'A_30(A\')_2 B_70(B\')_2 C_30(C\')_2'

# resolution = Property()
# resolution.name = 'NPW'
# resolution.vectors = [32,32,32]
# FreeEnergy = Property()
# FreeEnergy.name = 'Free Energy'
# FreeEnergy.scaler = 3.2



chemical_system.properties = [resolution,FreeEnergy]

test = pif.dumps(chemical_system, indent=4)





for file in listoffiles:
    start = time.time()

    op = open(file,'r')
    text = op.read().splitlines()
    op.close()
    DataDictionary = OutputParser(text)
    chemical_system = ChemicalSystem()


    end = time.time()
    print(end - start)


    
    # break
    # if i==3:
        

    # * A simulation cell has been initialized in 3 dimensions
    






