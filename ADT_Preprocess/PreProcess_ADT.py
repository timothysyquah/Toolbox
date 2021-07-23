#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 14:41:22 2021

@author: tquah
"""

import numpy as np
import matplotlib.pyplot as plt
import sys
import pandas as pd
import os
IDIR = os.getcwd()

#variables that have to be replaced by args
path = '/home/tquah/IMPORT_POD/L_116_EMADT_finitesize_larger/'
path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/rho_1_z_2_point_invest/CL/f0.5/chi0.2/nsc10.0/nbb99.0/LAM3DPhase/L_60'
operatorsfile = 'operators.dat'
dtfile = 'ADT_tsteps.dat'


os.chdir(path)


#load data in 
f = open(operatorsfile, 'r')
line1 = f.readline()[1:]
df = pd.read_csv(f, sep='\s+', names=line1.replace(' #', '').split(), dtype=np.float)
dtarray = np.loadtxt(dtfile)
operator_array = np.array(df)

#figure out how big blocks are
blockssize = int(df['step'][1])


#determine number of blocks
shape = np.shape(dtarray)
numofblocks = int(shape[0]/blockssize)
cutoff = numofblocks*blockssize

#get min dt


reshape_array = dtarray[0:cutoff,0].reshape(numofblocks,blockssize)
blocktime = np.sum(reshape_array, axis=1).reshape(1,numofblocks)
totaltime = np.sum(dtarray[:cutoff,0])

reweight = blocktime.transpose()*operator_array[1:,1:]
operator_array[1:,1:] = reweight

header = list(df)
header.append('BlockTime')
string = ''
for i in range(len(header)):
    string+=header[i]+' '

# line1+=' TotalTime'
blocktime_include = np.zeros((np.shape(blocktime)[1]+1,1))
blocktime_include[1:,0] = blocktime[0,:]
newoperators = np.hstack((operator_array,blocktime_include))
np.savetxt('weightedoperators.dat',newoperators,header = string)
#check


test = np.mean(df['StressXX.Real'])
test1 = np.mean(operator_array[:,3])*len(blocktime[0])/np.sum(blocktime)



loadweight = np.loadtxt('RW_operators.dat')
blockload = np.loadtxt('BlockSteps.dat')


test2 = np.mean(loadweight[:,3])*len(blockload)/np.sum(blockload)
test3 = np.mean(loadweight[:,3])/np.mean(blockload)

os.chdir(IDIR)

# timeinblocks = 200




# data = np.loadtxt(path)
# shape = np.shape(data)



# numofblocks = int(shape[0]/timeinblocks)
# cutoff = numofblocks*timeinblocks


# fields = int(shape[1]/2)

# # for i in range(fields):
#     # data[2*i+1,:]

# see = data[0:cutoff,0].reshape(timeinblocks,numofblocks)
# means = np.sum(see, axis=0)