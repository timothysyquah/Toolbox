#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 13:17:11 2020

@author: tquah
"""
#subprocess.call('qstat -u tquah'.split())

import numpy as np
import os
import time
import shutil
from make_input_file import make_input
#Functions


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
    
    

###############################################################################
Sstart = 49
Send = 300
n = 50
chi12 = 0
chi23 = 0.1
chi13 = 0
stats = 'CGC'
chiN = [chi12,chi13,chi23]
backbone_composition = [1.0]
backbone_species = [1]
backbone_statistics = stats
sidearm_composition_1 = [1.0]
sidearm_species_1 = [2]
sidearm_length_1 = 5.0
sidearm_statistics_1 = stats
sidearm_coverage_1 = 0.5
spacing_1 = 1
sidearm_composition_2 = [1.0]
sidearm_species_2 = [3]
sidearm_length_2 = 5.0
sidearm_statistics_2 = stats
sidearm_coverage_2 = 1-sidearm_coverage_1
spacing_2 = 1
###############################################################################
#alter params
dS = 1.0  #countour length
DT = 0.01  #time length
lambda_force = [0.05,0.01,0.5]
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
field_file = 'fields.in'
density_file = 'density.dat'
status_file = 'STATUS'
initial_field_dir = 'seed'
start_wait_interval = 1 #check existance in s
job_wait = 1 #checks if job is complete 
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



while i!=len(Ssweep):
    #chi_N = chi23*Ssweep[i]
    backbone_list = [backbone_statistics,Ssweep[i],backbone_species,\
         backbone_composition]

    sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
             sidearm_composition_1,sidearm_coverage_1,spacing_1]

    sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
                 sidearm_composition_2,sidearm_coverage_2,spacing_2]
    
    chain_list = [backbone_list,sidechain_1_list,sidechain_2_list]

     #number of changable parameters
    field_wait_count = 0
    statusread=0
    if i==0:
        field_DIR = os.path.join(IDIR,initial_field_dir)
        jobstatus = 1
    else:
        field_DIR = os.path.join(IDIR,past_field_dir[i-1])

    if jobstatus:
        if os.path.exists(field_DIR):
            status_path = os.path.join(field_DIR,status_file)
            density_path = os.path.join(field_DIR,density_file)
            field_path = os.path.join(field_DIR,field_file)
            print('Field Exist')
            
#            if os.path.exists(status_path):
#                co = open(status_path,'r')
#                statusread=int(co.read())
#                co.close()
#                print('Status Exist')
                
            statusread = 3
        
            if int(statusread)!=2:
                shutil.rmtree(field_DIR)
#                del past_field_dir[-1]
                i-=1
                reparam_itter_store_count[i]+=1
                print('Did Not Converge')
                continue
                            
#            else:
#                jobstatus=0
#                if past_field_dir[-1]!=WDIR:
#                    past_field_dir.append(WDIR)
#                continue
            
        else:
            continue
        
    
        
        
        if reparam_itter_store_count[i]==len(param_full_order):
            break
        out_DIR = os.path.join(field_DIR,outfile)
    
        Box_initial = 10

        if reparam_itter_store_count[i]==0:
            paramlist = [dS,DT,lambda_stress,lambda_force]
            replace_list.append([field_path,chain_list,chiN,npw,\
                                 Box_initial,d,jobname,paramlist[0],\
                                 paramlist[1],paramlist[2],paramlist[3]])    
 
        else:
            paramlist[1] = 0.5*paramlist[1]
            replace_list.append([field_path,chain_list,chiN,npw,\
                                 Box_initial,d,jobname,paramlist[0],\
                                 paramlist[1]*tune,paramlist[2],paramlist[3]])    
    else:
        continue

    WDIR_NAME = 'L_'+str(Ssweep[i])
    WDIR = os.path.join(IDIR,WDIR_NAME)

    if os.path.exists(WDIR):
        jobstatus=1
        print("{} exists...skipping.".format(WDIR_NAME)) 
        past_field_dir.append(WDIR)
        i+=1
        continue
    else:
        print('past')
        jobstatus=1
        os.makedirs(WDIR)
        make_input(phase,replace_list[i],WDIR)
        statuspath = os.path.join(WDIR,'STATUS')
        so = open(statuspath,'w+')
        so.write('0')
        so.close()
        fieldspath = os.path.join(WDIR,'fields.in')
        fo = open(fieldspath,'w+')
        fo.write('0')
        fo.close()
        i+=1

        
