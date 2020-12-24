#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 23:53:05 2020

@author: tquah
"""

import sys
path = '/home/tquah/tools/geninputfile/' 
sys.path.append(path)
from Chains import Lazy_Input_Generator_General
#from string import 
    

#Bottlebrush Build

chain_label = 'AB-Bottlebrush'
backbone_composition = [1.0]
backbone_species = [1]
backbone_length = 19.0
backbone_statistics = 'DGC'

backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                 backbone_composition]

sidearm_composition = [1.0]
sidearm_species = [2]
sidearm_length = 5.0
sidearm_statistics = 'DGC'
sidearm_coverage1 = 0.5
spacing = 1
sidechain_1_list = [sidearm_statistics,sidearm_length,sidearm_species,\
                 sidearm_composition,sidearm_coverage1,spacing]

sidearm_composition = [1.0]
sidearm_species = [3]
sidearm_length = 5.0
sidearm_statistics = 'DGC'
sidearm_coverage = 1-sidearm_coverage1
sidechain_2_list = [sidearm_statistics,sidearm_length,sidearm_species,\
                 sidearm_composition,sidearm_coverage,spacing]

chain_list = [backbone_list,sidechain_1_list,sidechain_2_list]

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
chiN = [0,0,0.1]
force_scale = [0.1]
# Make Input File
Lazy_Input_Generator_General(test_file,field,chain_list,chiN,dS,npw,dt,\
                         initial_box,stress_scale,force_scale,d,chain_label,\
                         diffuser_method='SOS',Nref=1,invzeta=0.1,\
                         kuhn_length=[1.0],stress_tol=1e-4,\
                         force_tol=1e-5)
