#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 23:53:05 2020

@author: tquah
"""

import sys
path = '/home/tquah/tools/geninputfile/' 
sys.path.append(path)
from Chains import Lazy_Input_Generator
#from string import 
    

#Bottlebrush Build

chain_label = 'AB-Bottlebrush'
backbone_composition = [0.25,0.5,0.25]
backbone_species = [1,2,3]
backbone_length = 99.0
backbone_statistics = 'DGC'

backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                 backbone_composition]

sidearm_composition_1 = [1.0]
sidearm_species_1 = [1]
sidearm_length_1 = 4.0
sidearm_statistics_1 = 'DGC'
sidearm_coverage_1 = 0.25
spacing_1 = 1
sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
                 sidearm_composition_1,sidearm_coverage_1,spacing_1]

sidearm_composition_2 = [1.0]
sidearm_species_2 = [2]
sidearm_length_2 = 4.0
sidearm_statistics_2 = 'DGC'
sidearm_coverage_2 = 0.5
spacing_2 = 1
sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
                 sidearm_composition_2,sidearm_coverage_2,spacing_2]

sidearm_composition_3 = [1.0]
sidearm_species_3 = [3]
sidearm_length_3 = 4.0
sidearm_statistics_3 = 'DGC'
sidearm_coverage_3 = 0.25
spacing_3 = 1
sidechain_3_list = [sidearm_statistics_3,sidearm_length_3,sidearm_species_3,\
                 sidearm_composition_3,sidearm_coverage_3,spacing_3]


    
chain_list = [backbone_list,sidechain_1_list,sidechain_2_list,sidechain_3_list]

#dynamic parameters
dt = 0.1
stress_scale = 25/dt
#parameters 
field = './fields.in'
test_file = './LAMtest.in'
dS = 1
d = 1
initial_box = 20
npw = 128
k = 127
chiN = [0,0,0.1]
force_scale = [0.1]
ends = True
# Make Input File
Lazy_Input_Generator(test_file,field,chain_list,chiN,dS,npw,dt,\
                         initial_box,stress_scale,force_scale,d,chain_label,\
                         ends,diffuser_method='SOS',Nref=1,invzeta=0.1,\
                         kuhn_length=[1.0],stress_tol=1e-4,\
                         force_tol=1e-5,CellUpdater = 'Broyden',add_phase = True,\
                         space_group = 'Fd-3m:1',non_primitive_centering = 'True',\
                         symmetrize = None,parallel_cuda = 6,cuda_thread_block_size = 64,\
                         nThreads = 1)