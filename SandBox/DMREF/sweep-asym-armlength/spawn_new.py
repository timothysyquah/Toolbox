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
    
def fA_calc(nA,nB,NA,NB):
    return (nA*(NA+1))/((nA*(NA+1))+(nB*(NB+1)))
def nA_calc(fA,nTot,NA,NB):
    return (fA*nTot*(NB+1))/((1-fA)*(NA+1)+fA*(NB+1))


def check_equality(a,b,message):
    assert a==b, message




if __name__=="__main__":
#-------------------------------------------------------------------------------#    
    PolyFTS_Path = '/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release'
    chain_label = 'ABC-Bottlebrush'

#-------------------------------------------------------------------------------#    
#Phase Specific Stuff and Default Stuff

    PhaseList      = ['GYR', 'HEX', 'BCC','LAM','DIS'] 
    InitialL0Guess = [5.0, 5.0,3.5,3.0,3.0]
    NThreads       = [ 1, 1,1,1,1]
    SpaceGroup     = ['Ia-3d','p6mm','Im-3m',None,None]
    NPW = [8,16,8,128,128]
#-------------------------------------------------------------------------------#    
# Chain Specifics
    Nref = 1
    narmtotal = 100
    Nsc_A_min = 20
    Nsc_A_max = 30
    Delta_Nsc_A = 5
    Nsc_total = 40
    
    
    fA_min = 0.02
    fA_max = 0.98
    Delta_fA = 0.06
    chiAB_min = 0.1
    chiAB_max = 0.1
    Delta_chiAB = 0.05
#-------------------------------------------------------------------------------#  
#simulation Specs
    dt = 0.01
    dS = 0.01 # if using CGC
    stress_scale = 0.1
    force_scale = [1.0,1.0]
#-------------------------------------------------------------------------------#  
    backbone_statistics = 'DGC'
    sidearm_statistics_1 = 'DGC'
    sidearm_statistics_2 = 'DGC'

    backbone_species = [1,2]
    sidearm_species_1 = [1]
    sidearm_species_2 = [2]
    backbone_length = narmtotal-1
    
    sidearm_composition_1 = [1.0]
    sidearm_composition_2 = [1.0]
    
    spacing_1 = 1
    spacing_2 = 1
    
    fieldpath = ''
    
    
    IDIR = os.getcwd()
    total_species = \
        list(set(backbone_species+sidearm_species_1+sidearm_species_2))

    try: 
        assert(len(PhaseList) == len(InitialL0Guess) == len(NThreads))
    except:
        if type(NThreads) == int and type(InitialL0Guess) == float and type(PhaseList) == str:
            pass
        else:
            raise RuntimeError("Lengths of PhasesList, InitialL0Guess and NThreads not equal! (%d != %d = %d)" % (len(PhaseList),len(InitialL0Guess),len(NThreads)) )

    fA_array = np.arange(fA_min,fA_max+1e-6,Delta_fA)
    Nsc_A_array = np.arange(Nsc_A_min,Nsc_A_max+1e-6,Delta_Nsc_A)
    chiAB_array = np.arange(chiAB_min,chiAB_max+1e-6,Delta_chiAB)

    print('chiAB Nsc_A Nsc_B fA')
    for chiAB in chiAB_array:
        for Nsc_A in Nsc_A_array:
            Nsc_B = Nsc_total-Nsc_A
            print('-------------------------')

            nA_check = nA_calc(fA_array,narmtotal,Nsc_A,Nsc_B)
            nA_int = np.array(list(set(list(np.around(nA_check))+list([narmtotal]))))
            del_loc = np.where(nA_int==narmtotal)[0]
            nA_int = np.delete(nA_int,del_loc)
            nA_int = np.sort(nA_int)
            fA_act = fA_calc(nA_int,narmtotal-nA_int,Nsc_A,Nsc_B)

            for i in range(0,len(fA_act)):
                print(f'{chiAB}   {Nsc_A}  {Nsc_B} {fA_act[i] : 0.5f}')
                k = 0
                for Phase in PhaseList:

                    if Phase=='LAM' or Phase=='DIS':
                        d = 1
                    elif Phase=='ACYL' or Phase=='HEX':
                        d = 2
                    else:
                        d = 3
                        
                    WDIR=f'chiAB_{chiAB}/NscA_{Nsc_A}_NscB_{Nsc_B}/fA{fA_act[i]:0.5f}/{Phase}Phase/'
                    # print(WDIR)
                    if os.path.exists(WDIR):
                        print("{} exists...skipping.".format(WDIR)) 
                        continue
                    else:
                        os.makedirs(WDIR)
                        
                    sidearm_length_1 = Nsc_A
                    sidearm_length_2 = Nsc_B
                    sidearm_coverage_1 = nA_int[i]/narmtotal
                    sidearm_coverage_2 = 1-sidearm_coverage_1
                    
                    backbone_composition = [sidearm_coverage_1,sidearm_coverage_2]
    
    
                    backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                                      backbone_composition]
                    chain_list = [backbone_list]
    
                    if Nsc_A>0 or Nsc_B>0:
                        sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
                                          sidearm_composition_1,sidearm_coverage_1,spacing_1]
                        
                        sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
                                          sidearm_composition_2,sidearm_coverage_2,spacing_2]\
                            
                        sidechain_list = [sidechain_1_list,sidechain_2_list]
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
        
                    REPLACE_lIST = [fieldpath,chain_list,[chiAB],d,chain_label,dS,dt,stress_scale,\
                                    force_scale,InitialL0Guess[k],SpaceGroup[k],[NPW[k]],Nref]
                    make_input(Phase,REPLACE_lIST,WDIR)
                    
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
                    k+=1

    

                    
