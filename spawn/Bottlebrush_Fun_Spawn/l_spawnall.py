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
chi12 = 0.1
stats = 'CGC'
Sstart = 130
Send = 500
n = 10
chiN = [chi12]
backbone_composition = [0.5,0.5]
backbone_species = [1,2]
backbone_statistics = stats
###############################################################################
#alter params
dS = 1.0  #countour length
DT = 0.002  #time length
lambda_force = [1.0,0.002]
lambda_stress = 0.001
###############################################################################
#other settings
base = 27
user='tquah'
npw = 1024 #number of plane waves
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




while i!=len(Ssweep):
    #chi_N = chi23*Ssweep[i]
    backbone_list = [backbone_statistics,Ssweep[i],backbone_species,\
         backbone_composition]    
    chain_list = [backbone_list]

    lambda_stress = 5/DT
    field_wait_count = 0
    statusread=0
    if i==0:
        field_DIR = os.path.join(IDIR,initial_field_dir)
        jobstatus=1
    else:
        field_DIR = os.path.join(IDIR,past_field_dir[i-1])

    if jobstatus:
        if os.path.exists(field_DIR):
            status_path = os.path.join(field_DIR,status_file)
            density_path = os.path.join(field_DIR,density_file)
            field_path = os.path.join(field_DIR,field_file)
    
            if os.path.exists(status_path):
                co = open(status_path,'r')
                statusread=int(co.read())
                co.close()
                
            
        
                if int(statusread)!=2:
                    shutil.rmtree(field_DIR)
                    del past_field_dir[-1]
                    i-=1
                    reparam_itter_store_count[i]+=1
                    time.sleep(60)
                    
                    continue
                
                if statusread==2:
                    density_data = np.loadtxt(density_path)
                    amp_check = abs(np.min(density_data[:,2])-np.max(density_data[:,2]))
                    if amp_check<tol:
                        shutil.rmtree(field_DIR)
                        i-=1
                        reparam_itter_store_count[i]+=1
                        del past_field_dir[-1]
                        time.sleep(60)
                        
                        continue
            
            else:
                jobstatus=0
                jobid = qsub_py(submitfile,WDIR,IDIR)
                if past_field_dir[-1]!=WDIR:
                    past_field_dir.append(WDIR)

                continue
            
        else:
            time.sleep(start_wait_interval)
            continue
        
    
        
        
        if reparam_itter_store_count[i]==len(param_full_order):
            break
        out_DIR = os.path.join(field_DIR,outfile)
        fo = open(out_DIR,'r')
        content = fo.read().split('\n')
        fo.close()
    
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        Box_initial = float(r[r.find("(")+1:r.find(")")])

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
        time.sleep(job_wait)
        jobstatus = job_status(user,jobid)
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
        jobstatus=0
        os.makedirs(WDIR)
        make_input(phase,replace_list[i],WDIR)
        make_submit(submitfile,phase,jobname,WDIR)
        jobid = qsub_py(submitfile,WDIR,IDIR)
        past_field_dir.append(WDIR)
        i+=1
        print('Running--'+WDIR_NAME+'--'+jobid)
