#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 15 22:02:01 2020

@author: tquah
"""

import numpy as np

import numpy as np
import os
import subprocess
import time
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

stats = 'DGC'
backbone_composition = [0.5,0.5]
backbone_species = [1,2]
backbone_statistics = stats
sidearm_composition_1 = [1.0]
sidearm_species_1 = [1]
sidearm_statistics_1 = stats
sidearm_coverage_1 = 0.5
spacing_1 = 1
sidearm_composition_2 = [1.0]
sidearm_species_2 = [2]
sidearm_statistics_2 = stats
sidearm_coverage_2 = 1-sidearm_coverage_1
spacing_2 = 1
total_species = list(set(backbone_species+sidearm_species_1+sidearm_species_2))
###############################################################################
#alter params
#Ssweep = np.array(list(np.arange(110,150+1e-6,2))+\
#                  list(np.arange(154,200+1e-6,4))+\
#                  list(np.arange(200,250+1e-6,6))+\
#                  list(np.arange(250,300+1e-6,8))+\
#                  list(np.arange(300,500+1e-6,10)))-1




dS = 0.1  #countour length
DT = 0.001  #time length
lambda_force = [1.0,0.01]
lambda_stress = 0.001
npw = 512 #number of plane waves
d = 1
paramlist = [dS,DT,lambda_stress,lambda_force]
###############################################################################
#other settings
user='tquah'
fp =10 #floating point
phase = 'LAM'
submitfile='./submit.sh'
jobname = 'AB-Bottlebrush'
field_file = 'fields_k.dat'
density_file = 'density.dat'
status_file = 'STATUS'
initial_field_dir = 'seed_L'
start_wait_interval = 15 #check existance in s
job_wait = 15 #checks if job is complete 
max_reparam_itter = 20
tol = 1e-3 #tolerance for amplitude test
tune=0.5 #how aggressive is the change
sleepcount = 30
chain_seed = True
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
chi = [0.05,0.1,0.15]
backbone = [50,100,150]

chain_seed_dict = dict()

chiN_initial = 25


for m in range(0,len(chi),1):
    chain_seed_dict[m]= dict()
    dir_chi = f'chi_{chi[m]}'
    chipath = os.path.join(IDIR,dir_chi)
    if os.path.isdir(chipath)==False:
        os.mkdir(chipath)

    chiN = [chi[m]]
    for n in range(0,len(backbone),1):
        dir_back = f'L_{backbone[n]}'
        backpath = os.path.join(chipath,dir_back)
        if os.path.isdir(backpath)==False:    
            os.mkdir(backpath)
        i = 0   
        
        N_initial = int(np.around(chiN_initial/(chi[m]*(backbone[n]))))
        
        
        Ssweep = np.arange(N_initial,0,-1,dtype = int)

        while i!=len(Ssweep):
            armlength = Ssweep[i]
            # check_equality(Ssweep[i]%2,1,'S values must be odd')
        
            backbone_list = [backbone_statistics,backbone[n],backbone_species,\
                  backbone_composition]
            chain_list = [backbone_list]
            sidearmlist = []
            if armlength>0:
                sidechain_1_list = [sidearm_statistics_1,armlength,sidearm_species_1,\
                          sidearm_composition_1,sidearm_coverage_1,spacing_1]
            
                sidechain_2_list = [sidearm_statistics_2,armlength,sidearm_species_2,\
                              sidearm_composition_2,sidearm_coverage_2,spacing_2]
                sidearmlist = [sidechain_1_list,sidechain_2_list]
                
            chain_list+=sidearmlist
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
                seed_DIR = os.path.join(IDIR,initial_field_dir)
                field_DIR = os.path.join(seed_DIR,str(chi[m]))
                jobstatus=1                
                if chain_seed:
                    if m!=0 or n!=0:
                        print(m)
                        print(n)
                        print(chain_seed_dict[m])
                        if m==0:
                            field_DIR = chain_seed_dict[m][n-1] 
                        if m>0 and n==0:
                            field_DIR = chain_seed_dict[m-1][n] 
                        
                        if m>0 and n>0:
                            field_DIR = chain_seed_dict[m][n-1] 

                        
                
                
            else:
                field_DIR = past_field_dir[-1]
        
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
                            print("Convergance Fail")
                            i-=1
                            break
                        
                        if statusread==2:
                            
                            density_data = np.loadtxt(density_path)
                            amp_check = abs(np.min(density_data[:,2])-np.max(density_data[:,2]))
                            if amp_check<tol:
                                i-=1
                                print("Disorder Fail")
                                break
                    else:
                        jobstatus=0
                        jobid = qsub_py(submitfile,WDIR,backpath)
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
                time.sleep(job_wait)
                jobstatus = job_status(user,jobid)
                continue
        
            WDIR_NAME = 'alpha_'+str(Ssweep[i])
            WDIR = os.path.join(backpath,WDIR_NAME)
        
            if os.path.exists(WDIR):
                jobstatus=1
                print("{} exists...skipping.".format(WDIR_NAME)) 
                past_field_dir.append(WDIR)

                if chain_seed:    
                    if i==0:
                        chain_seed_dict[m][n] = WDIR
                i+=1
                continue
            else:
                jobstatus=0
                os.makedirs(WDIR)
                make_input(phase,replace_list,WDIR)
                make_submit(submitfile,phase,jobname,WDIR)
                jobid = qsub_py(submitfile,WDIR,IDIR)
                past_field_dir.append(WDIR)
                if chain_seed:    
                    if i==0:
                        chain_seed_dict[m][n] = WDIR
                i+=1
                print('Running--'+WDIR_NAME+'--'+jobid)



