#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:17:11 2020

@author: tquah
"""
#subprocess.call('qstat -u tquah'.split())

import numpy as np
import os
import subprocess
import time
import shutil
from make_input_file import make_input
#Functions

def qsub_py(JOB,WDIR,IDIR):
    os.chdir(WDIR)
    job = subprocess.run(['qsub',JOB],capture_output=True)
    os.chdir(IDIR)
    job_id = str(job.stdout)
    return job_id[job_id.find("b'")+1:job_id.find("si")][1:]

def make_submit(SUBMITFILE,PHASE,JOBNAME,WDIR):
    with open('%s/submit.sh' % WDIR, 'w') as fout:
        with open(SUBMITFILE,'r') as f:
            for line in f:
                line = line.replace('__JOBNAME__',JOBNAME)
                line = line.replace('__INPUT__',f'{PHASE}.in')
                line = line.replace('__OUTPUT__',f'{PHASE}.out')
                fout.write(line)
    fout.close
    f.close()
    
    
def job_status(USR,JOBID):
    job_call = subprocess.run(['qstat','-u',USR],capture_output=True)
    job_list = str(job_call.stdout)
    index = job_list.find(JOBID)
    if index==-1:
        return 1
    else:
        return 0

###############################################################################
Sstart = 49
Send = 300
n = 50
chi12 = [0.1]
stats = 'DGC'
chiN = [chi12]
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
DT = 1  #time length
lambda_force = [0.01,0.01]
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
i = 0
jobid = None
reparam_itter_store_count = np.zeros(len(Ssweep),dtype=int)
reparam_overal_count= 0
replace_list = []
outfile = phase+'.out'





field_path = './fields.in'
itter =0
while i!=len(Ssweep):
    #chi_N = chi23*Ssweep[i]
    backbone_list = [backbone_statistics,Ssweep[i],backbone_species,\
         backbone_composition]

    sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
             sidearm_composition_1,sidearm_coverage_1,spacing_1]

    sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
                 sidearm_composition_2,sidearm_coverage_2,spacing_2]
    
    chain_list = [backbone_list,sidechain_1_list,sidechain_2_list]
    lambda_stress = 5/DT

    field_wait_count = 0
    statusread=0

    if i==0:
        statusread = 2
        
    if reparam_itter_store_count[i-1]>2:
        statusread=2

    if int(statusread)!=2:
        print(statusread)
        shutil.rmtree(past_field_dir[-1])

        del past_field_dir[-1]
        i-=1
        reparam_itter_store_count[i]+=1
        print('Delete!')
        continue
         
    if reparam_itter_store_count[i]==len(param_full_order):
        break

    Box_initial = 10
    
    if reparam_itter_store_count[i]==0:
        paramlist = [dS,DT,lambda_stress,lambda_force]
        replace_list.append([field_path,chain_list,chiN,npw,\
                             Box_initial,d,jobname,paramlist[0],\
                             paramlist[1],paramlist[2],paramlist[3]])    
 
    else:
        DT=DT*tune 
        paramlist = [dS,DT,lambda_stress,lambda_force]
        replace_list.append([field_path,chain_list,chiN,npw,\
                             Box_initial,d,jobname,paramlist[0],\
                             paramlist[1],paramlist[2],paramlist[3]])    
    
    WDIR_NAME = 'L_'+str(Ssweep[i])
    WDIR = os.path.join(IDIR,WDIR_NAME)
    jobstatus=0
    os.makedirs(WDIR)
    make_input(phase,replace_list[i],WDIR)
    past_field_dir.append(WDIR)
    print('pass')
    i+=1
    itter=+1
