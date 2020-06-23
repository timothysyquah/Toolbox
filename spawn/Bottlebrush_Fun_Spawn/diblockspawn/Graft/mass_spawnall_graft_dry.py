#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 24 16:10:03 2020

@author: tquah
"""

import numpy as np
import os
import subprocess
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
# Ssweep = np.array(list(np.arange(50,100+1e-6,2))+\
#                   list(np.arange(104,150+1e-6,4))+\
#                   list(np.arange(154,196+1e-6,6))+\
#                   list(np.arange(204,250+1e-6,8))+\
#                   list(np.arange(250,400+1e-6,10)))-1

Ssweep  = np.array([10,20,30])-1
    
dS = 1.0  #countour length
DT = 0.002  #time length
lambda_force = [1.0,0.02]
lambda_stress = 0.001
npw = 1024 #number of plane waves
d = 1
force_tol = 1e-5
stress_tol = 5e-5
CellUpdater = 'Euler'
paramlist = [dS,DT,lambda_stress,lambda_force,[force_tol,stress_tol],CellUpdater]

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
outfile = phase+'.out'
armlength = [1,2,3]
chain_seed_dict = dict()
graft = [2,5,10]
chi = 0.1
chiN = [chi]


for m in range(0,len(armlength),1):
    chain_seed_dict[m] = dict()
    dir_arm = f'alpha_{armlength[m]}'
    armpath = os.path.join(IDIR,dir_arm)
    for n in range(0,len(graft),1):
        dir_graft = f'graft_{graft[n]}'
        graftpath = os.path.join(armpath,dir_graft)
    
    
        i = 0   
         
        while i!=len(Ssweep):
            check_equality(Ssweep[i]%2,1,'S values must be odd')
        
            backbone_list = [backbone_statistics,Ssweep[i],backbone_species,\
                  backbone_composition]
            chain_list = [backbone_list]
            sidearmlist = []
            if armlength[m]>0:
                sidechain_1_list = [sidearm_statistics_1,armlength[m],sidearm_species_1,\
                          sidearm_composition_1,sidearm_coverage_1,spacing_1]
            
                sidechain_2_list = [sidearm_statistics_2,armlength[m],sidearm_species_2,\
                              sidearm_composition_2,sidearm_coverage_2,spacing_2]
                sidearmlist = [sidechain_1_list,sidechain_2_list]
                
            chain_list+=(graft[n]*sidearmlist)
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
                field_DIR = os.path.join(seed_DIR,str(chiN[0]))
                jobstatus=1
                if chain_seed:
                    if m!=0 or n!=0:
                        if m==0:
                            field_DIR = chain_seed_dict[m][n-1] 
                        if m>0 and n==0:
                            field_DIR = chain_seed_dict[m-1][n] 
                        
                        if m>0 and n>0:
                            field_DIR = chain_seed_dict[m][n-1] 
                
                print(field_DIR)


                        
                        
                
                
            else:
                field_DIR = past_field_dir[-1]
            
            status_path = os.path.join(field_DIR,status_file)
            density_path = os.path.join(field_DIR,density_file)
            field_path = os.path.join(field_DIR,field_file)
        
                
            out_DIR = os.path.join(field_DIR,outfile)
        
            Box_initial = 10 #float(r[r.find("(")+1:r.find(")")])
    
            if reparam_itter_store_count==0:
                    replace_list=[field_path,chain_list,chiN,npw,\
                                          Box_initial,d,jobname]+paramlist
                                 
                    
            WDIR_NAME = 'L_'+str(Ssweep[i])
            WDIR = os.path.join(graftpath,WDIR_NAME)
            if chain_seed:    
                if i==0:
                    chain_seed_dict[m][n] = WDIR
                # print(WDIR)

                
            jobstatus=0
            past_field_dir.append(WDIR)
            i+=1