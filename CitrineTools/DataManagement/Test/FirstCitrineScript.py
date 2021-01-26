#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 12 00:11:03 2020

@author: tquah
"""


from citrine import Citrine
import numpy as np
import glob,re
import os,sys
import matplotlib.pyplot as plt
from collections import defaultdict
import pickle
from scipy.interpolate import CubicSpline,interp1d,UnivariateSpline,barycentric_interpolate
from copy import deepcopy
import scipy as sp
import scipy.stats
import itertools
from citrination_client import CitrinationClient
from os import environ

import json
apikey_path = os.path.join('/home/tquah/.citrine_api_key')
op = open(apikey_path,'r')
apikey = op.read().split('\n')
op.close()

client = CitrinationClient(apikey[0])


data_client = client.data



from pypif import pif
from pypif.obj import *


chemical_system = ChemicalSystem()
chemical_system.chemical_formula = 'A_30(A\')_2 B_70(B\')_2 C_30(C\')_2'

resolution = Property()
resolution.name = 'NPW'
resolution.vectors = [32,32,32]
FreeEnergy = Property()
FreeEnergy.name = 'Free Energy'
FreeEnergy.scaler = 3.2



chemical_system.properties = [resolution,FreeEnergy]

test = pif.dumps(chemical_system, indent=4)

op = open('test.json','w+')

op.write(test)
op.close()

# ... client initialization left out
data_client = client.data

file_path = "./test.json"
dataset_id = 195077
data_client.upload(dataset_id, file_path)

# def data_extraction(path,filename,keyword):
#     IDIR_og = os.getcwd()
#     os.chdir(path)
#     IDIR = os.getcwd()
#     dir_dict = dict()
#     # print(glob.glob(f'./Nsc*/fAfC*/{filename}'))
#     # print(glob.glob('./Nsc*'))
    
#     for file in glob.glob(f'./Nsc*/fA**/{filename}', recursive = True):
#         dir_split = file.split('/')
#         full_path = os.path.join(IDIR,file)
#         op = open(full_path,'r')
#         dataread = op.read().splitlines()
#         op.close()
#         dictionary_list = []
#         for i in range(0,len(keyword)-1):
#             dir1loc = dir_split.index([s for s in dir_split if keyword[i] in s][0])
s#             dir1=dir_split[dir1loc]
#             dictionary_list.append(float(re.findall("\d+\.\d+", dir1)[0]))
    
#         dir2loc = dir_split.index([s for s in dir_split if keyword[-1] in s][0])
#         dir2=dir_split[dir2loc]
#         innerstore = float(re.findall("\d+\.\d+", dir2)[0])
        
#         for i in range(0,len(dataread),1):
#             newlist = []
#             splitdata = dataread[i].split(' ')
#             phase = splitdata[0][:-5]
#             newlist+=dictionary_list
#             newlist.append(phase)      
#             if int(splitdata[2])==2:
                
                
#                 if tuple(newlist) in dir_dict:
#                     temparray = np.array([innerstore,float(splitdata[1])]) 
#                     dir_dict[tuple(newlist)] = np.vstack((dir_dict[tuple(newlist)],temparray))
#                 else:
#                     dir_dict[tuple(newlist)] = np.array([innerstore,float(splitdata[1])])
#             else:
#                 print('Simulations have STATUS 0 1 3...revist extractF0.dat...')
#                 breakuntitled1
#     os.chdir(IDIR_og)
#     return dir_dict


# def input_extraction(path,keyword,phases):
#     IDIR_og = os.getcwd()
#     os.chdir(path)
#     IDIR = os.getcwd()
#     dir_dict = dict()
#     # print(glob.glob(f'./Nsc*/fAfC*/{filename}'))
#     # print(glob.glob('./Nsc*'))
#     test = []
#     for phase in phases:
#         # print(glob.glob(f'./Nsc*/fA*/**/{phase}.out', recursive = True))
#         for file in glob.glob(f'./Nsc*/fA*/**/{phase}.out', recursive = True):
#             dir_split = file.split('/')
#             full_path = os.path.join(IDIR,file)
#             print(full_path)
#             op = open(full_path,'r')
#             dataread = op.read().splitlines()
#             op.close()
            
#             data_select = dataread[0:500]
#             del dataread
#             test.append(data_select)
            
#             break
#         break
#     os.chdir(IDIR_og)
#     return testuntitled1
# def extract_xy(data_dict):
#     l1 = []
#     l2 = []
#     for tup in list(data_dict):
#         l1.append(tup[0])
#         l2.append(tup[1])
#     return sorted(list(set(l1))),sorted(list(set(l2)))

# path_to_free_energies = '/home/tquah/Projects/DMREF/Discrete_Chains/sweep-fAfC-armlength-chiNAB13-k2.7/sweep-fAfC-armlength-chiNAB13-k2.7/PhaseDiagram/'



# keyword = ['Nsc','f']

# filename = 'F0_phases.dat'

# data_dict = data_extraction(path_to_free_energies,filename,keyword)


# Nsc,phase = extract_xy(data_dict)


# input_dict = input_extraction(path_to_free_energies,keyword,phase)







# chemical_system = ChemicalSystem()
# chemical_system.chemical_formula = 'MgO2'

# band_gap = Property()
# band_gap.name = 'Band gap'
# band_gap.scalars = 7.8
# band_gap.units = 'eV'

# chemical_system.properties = band_gap

# save = pif.dumps(chemical_system, indent=4)

# data_dict_test = dict()

# data_dict_test['Test'] =dict()

# data_dict_test['Test']['hello'] = 1
