#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Feb 20 10:39:43 2020

@author: tquah
"""
#subprocess.call('qstat -u tquah'.split())

import numpy as np
import os
from make_input_file import make_input
#Functions


###############################################################################
Sstart = 49
Send = 300
n = 50
chi12 = 0
chi23 = 0.1
chi13 = 0
chiN = [chi12,chi13,chi23]
stats = 'DGC'
backbone_composition = [0.5,0.5]
backbone_species = [1,2]
backbone_statistics = stats
sidearm_composition_1 = [1.0]
sidearm_species_1 = [1]
sidearm_length_1 = 5.0
sidearm_statistics_1 = stats
sidearm_coverage_1 = 0.5
spacing_1 = 1
sidearm_composition_2 = [1.0]
sidearm_species_2 = [2]
sidearm_length_2 = 5.0
sidearm_statistics_2 = stats
sidearm_coverage_2 = 1-sidearm_coverage_1
spacing_2 = 1
###############################################################################
#alter params
dS = 1.0  #countour length
DT = 0.01  #time length
lambda_force = [0.1,0.1]
lambda_stress = 0.001
###############################################################################
#other settings
base = 27
user='tquah'
npw = 128 #number of plane waves
fp =10 #floating point
d = 1
phase = 'LAM'
submitfile='./submit.sh'
jobname = 'AB-Bottlebrush'
field_file = 'fields.dat'
density_file = 'density.dat'
status_file = 'STATUS'
initial_field_dir = 'seed_L_39_chi0.1_alpha_5'
start_wait_interval = 15 #check existance in s
job_wait = 15 #checks if job is complete 
max_reparam_itter = 20
tol = 1e-3 #tolerance for amplitude test
tune=0.5 #how aggressive is the change
###############################################################################

#chisweep = np.round(np.linspace(chistart,chiend,n),2)
Ssweep = np.round(np.arange(Sstart,Send,n),5)

#need to populate

param_full_order = []

param_full_order = [1,1,1,1,1,1,1,1]

IDIR = os.getcwd()
WDIR=None
past_field_dir = []
jobid = None
reparam_itter_store_count = np.zeros(len(Ssweep),dtype=int)
reparam_overal_count= 0
replace_list = []
outfile = phase+'.out'


if os.path.exists('./L_49'):
    os.rmdir('./L_49')

i = 0
while i!=len(Ssweep):
    #chi_N = chi23*Ssweep[i]
    backbone_list = [backbone_statistics,Ssweep[i],backbone_species,backbone_composition]

    sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
         sidearm_composition_1,sidearm_coverage_1,spacing_1]

    sidechain_2_list = [sidearm_statistics_2,sidearm_length_1,sidearm_species_2,\
             sidearm_composition_2,sidearm_coverage_2,spacing_2]

    chain_list = [backbone_list,sidechain_1_list,sidechain_2_list]

    paramlist = [dS,DT,lambda_stress,lambda_force,True] #number of changable parameters
    replace_list = ['field.in',chain_list,chiN,npw,20,d,jobname]+paramlist
    field_wait_count = 0
    statusread=0
    WDIR_NAME = 'L_'+str(Ssweep[i])
    WDIR = os.path.join(IDIR,WDIR_NAME)
    os.makedirs(WDIR)
    make_input(phase,replace_list,WDIR)
    i+=1
    break
