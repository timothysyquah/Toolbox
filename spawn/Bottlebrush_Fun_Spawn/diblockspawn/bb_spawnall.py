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
chi12 = 1.0
stats = 'DGC'
chiN = [chi12]
backbone_composition = [0.5,0.5]
backbone_species = [1,2]
backbone_statistics = stats
sidearm_composition_1 = [1.0]
side_arms_length_sym = 15.0
sidearm_species_1 = [1]
sidearm_length_1 = side_arms_length_sym
sidearm_statistics_1 = stats
sidearm_coverage_1 = 0.5
spacing_1 = 1
sidearm_composition_2 = [1.0]
sidearm_species_2 = [2]
sidearm_length_2 = side_arms_length_sym
sidearm_statistics_2 = stats
sidearm_coverage_2 = 1-sidearm_coverage_1
spacing_2 = 1
total_species = list(set(backbone_species+sidearm_species_1+sidearm_species_2))
###############################################################################
#alter params
Ssweep = np.array([6,8,10,12,14,16,18,20,22,24,26,30,34,42,48,54,60,66,72,78,84,90])-1
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
max_reparam_itter = 20
tol = 1e-3 #tolerance for amplitude test
tune=0.5 #how aggressive is the change
sleepcount = 30
###############################################################################

#need to populate
IDIR = os.getcwd()
WDIR=None
past_field_dir = []
jobid = None
reparam_itter_store_count = 0
reparam_overal_count= 0
replace_list = []
outfile = phase+'.out'


i = 0 

while i!=len(Ssweep):
    check_equality(Ssweep[i]%2,1,'S values must be odd')

    backbone_list = [backbone_statistics,Ssweep[i],backbone_species,\
         backbone_composition]

    sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
             sidearm_composition_1,sidearm_coverage_1,spacing_1]

    sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
                 sidearm_composition_2,sidearm_coverage_2,spacing_2]
    
    chain_list = [backbone_list,sidechain_1_list,sidechain_2_list]
    check_equality(len(chain_list)==0,False,"Need more than 0 chains")
    check_equality(len(total_species),len(lambda_force),'Force and species should align')
    for j in range(0,len(chain_list),1):
        check_equality(len(chain_list[j])>6,False,"Too many side chain inputs")
        check_equality(len(chain_list)!=0,True,"Need more than 0 chains")
        list_of_checks = [[len(chain_list[j][2]),len(chain_list[j][3]),'check chain species and composition'],\
                            [np.isclose(np.sum(chain_list[j][3]),1),True,'backbone composition does not sum to 1']]

        for checkchecks in list_of_checks:
            check_equality(checkchecks[0],checkchecks[1],checkchecks[2])                 
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
