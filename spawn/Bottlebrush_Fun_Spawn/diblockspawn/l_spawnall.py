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
    
def check_equality(a,b,message):
    assert a==b, message


###############################################################################
chi12 = 1.0 #chivalues
stats = 'CGC' #statistics model
chiN = [chi12] #increase for more complex systems
backbone_composition = [0.5,0.5] #maxbone composition
backbone_species = [1,2] #species type
backbone_statistics = stats
Ssweep = np.array([14,16,18,20,22,24,26,30,34,36,42,48,54,60,66,72,78,84,90])-1
chiN_check = chiN*Ssweep
###############################################################################
#alter params
dS = 1.0  #countour length
DT = 0.002  #time length
lambda_force = [1.0,0.002]
lambda_stress = 0.001
npw = 1024 #number of plane waves
d = 1
paramlist = [dS,DT,lambda_stress,lambda_force]
###############################################################################
#other settings
user='tquah'
fp =10 #floating point
phase = 'LAM'
submitfile='./submit.sh'
jobname = 'AB-Bottlebrush'
field_file = 'fields.dat'
density_file = 'density.dat'
status_file = 'STATUS'
initial_field_dir = 'seed_L'
start_wait_interval = 15 #check existance in s
job_wait = 15 #checks if job is complete 
max_reparam_itter = 10
tol = 1e-3 #tolerance for amplitude test
tune=0.5 #how aggressive is the change
sleepcount = 30
###############################################################################

#chisweep = np.round(np.linspace(chistart,chiend,n),2)

#need to populate
IDIR = os.getcwd()
WDIR=None
past_field_dir = []
jobid = None
reparam_itter_store_count = 0
reparam_overal_count= 0
replace_list = []
outfile = phase+'.out'


#Some good ol prechecks
list_of_checks = [[len(backbone_composition),len(backbone_species),'check backbone species and composition'],\
                   [np.isclose(np.sum(backbone_composition),1),True,'backbone composition does not sum to 1'],\
                   [len(lambda_force),len(backbone_species),'Force and Backbone should align']]

for i in range(0,len(list_of_checks),1):
    check_equality(list_of_checks[i][0],list_of_checks[i][1],list_of_checks[i][2])

i = 0

while i!=len(Ssweep):
    
    
    check_equality(Ssweep[i]%2,1,'S values must be odd')

    backbone_list = [backbone_statistics,Ssweep[i],backbone_species,\
         backbone_composition]    
    chain_list = [backbone_list]
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
                    time.sleep(sleepcount)
                    shutil.rmtree(field_DIR)
                    del past_field_dir[-1]
                    i-=1
                    reparam_itter_store_count+=1
                    print()
                    continue
                
                if statusread==2:
                    
                    density_data = np.loadtxt(density_path)
                    amp_check = abs(np.min(density_data[:,2])-np.max(density_data[:,2]))
                    if amp_check<tol:
                        time.sleep(sleepcount)
                        shutil.rmtree(field_DIR)
                        i-=1
                        reparam_itter_store_count+=1
                        del past_field_dir[-1]
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
        if reparam_itter_store_count==max_reparam_itter:
            break
        out_DIR = os.path.join(field_DIR,outfile)
        fo = open(out_DIR,'r')
        content = fo.read().split('\n')
        fo.close()
    
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        Box_initial = float(r[r.find("(")+1:r.find(")")])

        if reparam_itter_store_count==0:
            replace_list=[field_path,chain_list,chiN,npw,\
                                 Box_initial,d,jobname,paramlist[0],\
                                 paramlist[1],paramlist[2],paramlist[3]]
 
        else:
            print('Reparameterization has been initiated now %d out of %d'%(reparam_itter_store_count,max_reparam_itter))
            paramlist[1] = tune*paramlist[1]
            replace_list = []
            replace_list=[field_path,chain_list,chiN,npw,\
                                 Box_initial,d,jobname,paramlist[0],\
                                 paramlist[1],paramlist[2],paramlist[3]]  
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
        make_input(phase,replace_list,WDIR)
        make_submit(submitfile,phase,jobname,WDIR)
        jobid = qsub_py(submitfile,WDIR,IDIR)
        past_field_dir.append(WDIR)
        i+=1
        print('Running--'+WDIR_NAME+'--'+jobid)
