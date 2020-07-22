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




if __name__=="__main__":
    PolyFTS_Path = '/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release'
    PhaseList      = ['AGYR32', 'ADIA32', 'ACYL'] 
    InitialL0Guess = [5.0, 5.0,3.5]
    NThreads       = [ -1, -1,1]
    SpaceGroup     = ['I4_1.32','Fd-3m:1','p4mm']



    Nscstart = 1
    Nscend = 1
    Nscdelta = 1
    fAfCstart = 0.19 # 0.2
    fAfCend = 0.30 # 0.2
    fAfCdelta = 0.01 # 0.01
    
    
    chain_label = 'ABC-Bottlebrush'
    dS = 0.01
    dt = 0.01
    stress_scale = 0.1
    force_scale = [1.0,1.0]


    #==========================================================================
    #Fixed Parameters
    spacing = 1
    ds = 0.1
    narmtotal = 100
    Nref = 100
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
    fAfCs = np.arange(fAfCstart,fAfCend+1e-6,fAfCdelta)
    
    try: 
        assert(len(PhaseList) == len(InitialL0Guess) == len(NThreads))
    except:
        if type(NThreads) == int and type(InitialL0Guess) == float and type(PhaseList) == str:
            pass
        else:
            raise RuntimeError("Lengths of PhasesList, InitialL0Guess and NThreads not equal! (%d != %d = %d)" % (len(PhaseList),len(InitialL0Guess),len(NThreads)) )

    IDIR = os.getcwd()
    seedpath = os.path.join(IDIR,'SEEDS')
    
    if os.path.exists(seedpath)!=True:
        raise RuntimeError('Seed Path Does not Exist')
    

    for i in range(0,len(Nscs),1):
        for j in range(0,len(fAfCs),1):
            for k in range(0,len(PhaseList),1):
                fieldname = PhaseList[k]+'_fields.in'
                fieldpath = os.path.join(seedpath,fieldname)
                if os.path.exists(fieldpath)!=True:
                    raise RuntimeError('Field Path Does not Exist')

                if PhaseList[k]=='LAM':
                    d = 1
                elif PhaseList[k]=='ACYL':
                    d = 2
                else:
                    d = 3
                

                
                WDIR=f'Nsc{Nscs[i]:0.3f}/fAfC{fAfCs[j]:0.2f}/{PhaseList[k]}Phase/'
                if os.path.exists(WDIR):
                    print("{} exists...skipping.".format(WDIR)) 
                    continue
                else:
                    os.makedirs(WDIR)

            
                sidearm_length_1 = Nscs[i]
                sidearm_length_2 = Nscs[i]
                sidearm_length_3 = Nscs[i]

                sidearm_coverage_1 = np.round(fAfCs[j],5)
                sidearm_coverage_2 = np.round(1-(2*fAfCs[j]),5)
                sidearm_coverage_3 = np.round(fAfCs[j],5)
                
                backbone_composition = [sidearm_coverage_1,sidearm_coverage_2,sidearm_coverage_3]


                backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                                  backbone_composition]
                chain_list = [backbone_list]

                if Nscs[i]>0:
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
                for l in range(0,len(chain_list),1):
                    check_equality(len(chain_list[l])>6,False,"Too many side chain inputs")
                    check_equality(len(chain_list)!=0,True,"Need more than 0 chains")
                    list_of_checks = [[len(chain_list[l][2]),len(chain_list[l][3]),'check chain species and composition'],\
                                        [np.isclose(np.sum(chain_list[l][3]),1),True,'backbone composition does not sum to 1']]
            
                    for checkchecks in list_of_checks:
                        check_equality(checkchecks[0],checkchecks[1],checkchecks[2])                 
    
                REPLACE_lIST = [fieldpath,chain_list,chiN,d,chain_label,dS,dt,stress_scale,\
                                force_scale,InitialL0Guess[k],SpaceGroup[k]]
                
                make_input(PhaseList[k],REPLACE_lIST,WDIR)
                
                
                if NThreads[k] == 1:
                    submitfile='SEEDS/submit.sh'
                elif NThreads[k] == -1:
                    if PhaseList[k] != 'SIGMA':
                        submitfile='SEEDS/submitGPU.sh'
                    else:
                        submitfile='SEEDS/submitGPU-P100.sh'
                else:
                    submitfile='SEEDS/submitPLL.sh'

                
                
                with open('%s/submit.sh' % WDIR, 'w') as fout:
                    with open(submitfile,'r') as f:
                        for line in f:
                            line = line.replace('__JOBNAME__',WDIR)
                            line = line.replace('__INFILE__',f'{PhaseList[k]}.in')
                            line = line.replace('__OUTFILE__',f'{PhaseList[k]}.out')
                            line = line.replace('__PFTPATH__',PolyFTS_Path)

                            if NThreads[k] > 1:
                                line = line.replace('__NTHREADS__',WDIR)
                            fout.write(line)
                            
                os.chdir(WDIR)
                # import subprocess
                # cmd="qsub submit.sh"
                # cmd="sbatch submit.sh"
                # subprocess.call(cmd.split())
                os.chdir(IDIR)
                # time.sleep(2)
                    
                    
