#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:08:53 2020

@author: tquah
"""
from make_input_file import make_input
import numpy as np
import os
import pdb
import shutil
import time
import subprocess
def qsub_py(JOB,WDIR,IDIR):
    os.chdir(WDIR)
    job = subprocess.run(['qsub',JOB],capture_output=True)
    os.chdir(IDIR)
    job_id = str(job.stdout)
    return job_id[job_id.find("b'")+1:job_id.find("si")][1:]
def job_status(USR,JOBID):
    job_call = subprocess.run(['qstat','-u',USR],capture_output=True)
    job_list = str(job_call.stdout)
    index = job_list.find(JOBID)
    if index==-1:
        return 1
    else:
        return 0

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
    
    

def check_equality(a,b,message):
    assert a==b, message

def reorder_seed(seedpath,header):
    if os.path.exists(seedpath):
        so = open(seedpath,'r')
        seed_locs = so.read()
        loc_list = list(set(seed_locs.splitlines()[1:]))
        loc_list.sort()
        so.close()
        so = open(seedpath,'w')
        so.write(header)
        so.write('\n')
        for i in range(0,len(loc_list),1):
            so.write(loc_list[i])
            so.write('\n')
        so.close()
    else:
        so = open(seedpath,'w')
        so.write(header)
        so.write('\n')
        so.close()


if __name__=="__main__":
    user = 'tquah'
    PolyFTS_Path = '/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release'
    field_file = 'fields_k.bin'
    status_file = 'STATUS'
    seedin_file = 'SEED_INFO.in'
    seed_header = 'Phase Nsc fAfC Path'
    PhaseList      = ['AGYR', 'ADIA', 'ACYL'] 
    InitialL0Guess = [5.0, 5.0,3.5]
    NThreads       = [ -1, -1,1]
    SpaceGroup     = ['I4_1.32','Fd-3m:1','p4mm']



    Nscstart = 0
    Nscend = 4
    Nscdelta = 1
    fAfC = 0.3
    
    chain_label = 'ABC-Bottlebrush'
    dS = 0.01
    dt = 0.002
    stress_scale = 1/dt
    force_scale = [1.0,1.0,1.0]


    #==========================================================================
    #Fixed Parameters
    spacing = 1
    ds = 0.1
    narmtotal = 100
    Nref = 1
    chiAB = 13.0/(spacing*(narmtotal-1))*Nref
    k = 2.7
    chiBC = chiAB
    chiAC = chiAB*k
    chiN = [np.round(chiAB,3),np.round(chiBC,3),np.round(chiAC,3)]

    #==========================================================================

    #chain default parameters
    backbone_statistics = 'DGC'
    sidearm_statistics_1 = 'DGC'
    sidearm_statistics_2 = 'DGC'
    sidearm_statistics_3 = 'DGC'

    backbone_species = [1,2,3]
    sidearm_species_1 = [1]
    sidearm_species_2 = [2]
    sidearm_species_3 = [3]
    backbone_length = narmtotal-1
    
    sidearm_composition_1 = [1.0]
    sidearm_composition_2 = [1.0]
    sidearm_composition_3 = [1.0]

    spacing_1 = 1
    spacing_2 = 1
    spacing_3 = 1
    #==========================================================================
    total_species = \
        list(set(backbone_species+sidearm_species_1+sidearm_species_2+sidearm_species_3))

    Nscs = np.arange(Nscstart,Nscend+1e-6,Nscdelta,dtype = int)
    # fAfCs = np.arange(fAfCstart,fAfCend+1e-6,fAfCdelta)
    try: 
        assert(len(PhaseList) == len(InitialL0Guess) == len(NThreads))
    except:
        if type(NThreads) == int and type(InitialL0Guess) == float and type(PhaseList) == str:
            pass
        else:
            raise RuntimeError("Lengths of PhasesList, InitialL0Guess and NThreads not equal! (%d != %d = %d)" % (len(PhaseList),len(InitialL0Guess),len(NThreads)) )

    IDIR = os.getcwd()
    seedpath = os.path.join(IDIR,'SEEDS')
    seedin_path = os.path.join(IDIR,seedin_file)
    if os.path.exists(seedpath)!=True:
        raise RuntimeError('Seed Path Does not Exist')
    reorder_seed(seedin_path,seed_header)
    for i in range(0,len(PhaseList),1):
        field_path_list = []
        for j in range(0,len(Nscs),1):
            if j ==0:
                fieldname = PhaseList[i]+'_fields.in'
                fieldpath = os.path.join(seedpath,fieldname)
            if j>0:
                fieldpath = field_path_list[-1]
                
            if os.path.exists(fieldpath)!=True:
                print('Field Path Does not Exist')
                break
                #raise RuntimeError('Field Path Does not Exist')

            if PhaseList[i]=='LAM':
                d = 1
            elif PhaseList[i]=='ACYL':
                d = 2
            else:
                d = 3
            

            
            r_WDIR=f'Nsc{Nscs[j]:0.3f}/fAfC{fAfC:0.2f}/{PhaseList[i]}Phase/'
            WDIR = os.path.join(IDIR,r_WDIR)
            if os.path.exists(WDIR):
                print("{} exists...skipping.".format(WDIR)) 
                status_path = os.path.join(WDIR,status_file)
                field_path = os.path.join(WDIR,field_file)
                if os.path.exists(status_path):
                    co = open(status_path,'r')
                    statusread=int(co.read())
                    co.close()
                statusread=2
                if statusread!=2:
                    print('Simulation not converged quitting...')
                    break
                if statusread==2:
                    field_path_list.append(field_path)
                    so = open(seedin_path,'a')
                    so.write(f'{PhaseList[i]} {Nscs[j]:0.3f} {field_path} \n')
                    so.close()
                continue
            else:
                os.makedirs(WDIR)

        
            sidearm_length_1 = Nscs[j]
            sidearm_length_2 = Nscs[j]
            sidearm_length_3 = Nscs[j]

            sidearm_coverage_1 = np.round(fAfC,5)
            sidearm_coverage_2 = np.round(1-(2*fAfC),5)
            sidearm_coverage_3 = np.round(fAfC,5)
            
            backbone_composition = [sidearm_coverage_1,sidearm_coverage_2,sidearm_coverage_3]


            backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                              backbone_composition]
            chain_list = [backbone_list]

            if Nscs[j]>0:
                sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
                                  sidearm_composition_1,sidearm_coverage_1,spacing_1]
                
                sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
                                  sidearm_composition_2,sidearm_coverage_2,spacing_2]
                
                sidechain_3_list = [sidearm_statistics_3,sidearm_length_3,sidearm_species_3,\
                                  sidearm_composition_3,sidearm_coverage_3,spacing_3]
                sidechain_list = [sidechain_1_list,sidechain_2_list,sidechain_3_list]
                chain_list+=sidechain_list
            #pass internal checks
            check_equality(len(chain_list)==0,False,"Need more than 0 chains")
            for k in range(0,len(chain_list),1):
                check_equality(len(chain_list[k])>6,False,"Too many side chain inputs")
                check_equality(len(chain_list)!=0,True,"Need more than 0 chains")
                list_of_checks = [[len(chain_list[k][2]),len(chain_list[k][3]),'check chain species and composition'],\
                                    [np.isclose(np.sum(chain_list[k][3]),1),True,'backbone composition does not sum to 1']]
        
                for checkchecks in list_of_checks:
                    check_equality(checkchecks[0],checkchecks[1],checkchecks[2])                 

            REPLACE_lIST = [fieldpath,chain_list,chiN,d,chain_label,dS,dt,stress_scale,\
                            force_scale,InitialL0Guess[i],SpaceGroup[i]]
            
            make_input(PhaseList[i],REPLACE_lIST,WDIR)
            
            if NThreads[i] == 1:
                submitfile='SEEDS/submit.sh'
            elif NThreads[i] == -1:
                if PhaseList[i] != 'SIGMA':
                    submitfile='SEEDS/submitGPU.sh'
                else:
                    submitfile='SEEDS/submitGPU-P100.sh'
            else:
                submitfile='SEEDS/submitPLL.sh'
            
            
            
            with open('%s/submit.sh' % WDIR, 'w') as fout:
                with open(submitfile,'r') as f:
                    for line in f:
                        line = line.replace('__JOBNAME__',WDIR)
                        line = line.replace('__INFILE__',f'{PhaseList[i]}.in')
                        line = line.replace('__OUTFILE__',f'{PhaseList[i]}.out')
                        line = line.replace('__PFTPATH__',PolyFTS_Path)

                        if NThreads[i] > 1:
                            line = line.replace('__NTHREADS__',WDIR)
                        fout.write(line)
                        
            # import subprocess
            # cmd="qsub submit.sh"
            # cmd="sbatch submit.sh"
            # subprocess.call(cmd.split())
            jobid = qsub_py('submit.sh',WDIR,IDIR)
            time.sleep(5)
            jobstatus=0
            jobstatus = job_status(user,jobid)

            
            while jobstatus==0:
                time.sleep(30)
                jobstatus = job_status(user,jobid)
            
            status_path = os.path.join(WDIR,status_file)
            field_path = os.path.join(WDIR,field_file)
            if os.path.exists(status_path):
                co = open(status_path,'r')
                statusread=int(co.read())
                co.close()
                
            if statusread!=2:
                print('Simulation not converged quitting...')
                break
            if statusread==2:
                field_path_list.append(field_path)
                so = open(seedin_path,'a')
                so.write(f'{PhaseList[i]} {Nscs[j]:0.3f} {fAfC} {field_path} \n')
                so.close()
    
    reorder_seed(seedin_path,seed_header)
