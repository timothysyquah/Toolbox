#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Mar 30 15:08:53 2020

@author: tquah
"""
# from make_input_file_v2 import make_input
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
# def nA_calc(fA,Neff,NA,NB):
#     return Neff*fA/(NA+1)

def fA_tot(nA,nB,NA,NB):
    return (nA*(NA+1))+(nB*(NB+1))


def nA_calc(fA,Neff,NA):
    return Neff*fA/(NA+1)


def Neff(nA,nB,NA,NB):
    return nA*(NA+1)+nB*(NB+1)


def check_equality(a,b,message):
    assert a==b, message




if __name__=="__main__":
#-------------------------------------------------------------------------------#    
    PolyFTS_Path = '/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release'
    chain_label = 'ABC-Bottlebrush'

#-------------------------------------------------------------------------------#    
#Phase Specific Stuff and Default Stuff

#    PhaseList      = ['GYR', 'HEX', 'BCC','LAM','DIS'] 
#    InitialL0Guess = [10.0, 7.0,10.0,7.0,5.0]
#    NThreads       = [ 1, 1,1,1,1]
#    SpaceGroup     = ['Ia-3d','p6mm','Im-3m',None,None]
#    NPW = [32,32,32,128,128]

#    PhaseList      = ['GYR', 'HEX', 'LAM','DIS'] 
#    InitialL0Guess = [10.0, 7.0,7.0,5.0]
#    NThreads       = [ 1, 1,1,1]
#    SpaceGroup     = ['Ia-3d','p6mm',None,None]
#    NPW = [32,32,128,128]

#    PhaseList      = ['GYR'] 
#    InitialL0Guess = [[10.0]]
#    NThreads       = [1]
#    SpaceGroup     = ['Ia-3d']
#    NPW = [[32]]

#    PhaseList      = ['FCC'] 
#    InitialL0Guess = [[5.0]]
#    NThreads       = [1]
#    SpaceGroup     = ['Fm-3m']
#    NPW = [[32]]

#    PhaseList      = ['LAM','DIS'] 
#    InitialL0Guess = [5.0,5.0]
#    NThreads       = [ 1,1]
#    SpaceGroup     = [None,None]
#    NPW = [128,128]

#    PhaseList      = ['A1596'] 
#    InitialL0Guess = [[12.133]]
#    NThreads       = [-1]
#     SpaceGroup     = ['Pm-3n']
#    NPW = [[96]]

    PhaseList      = ['SIGMA192'] 
    InitialL0Guess = [[21.5,21.5,11.4]]
    NThreads       = [-1]
    SpaceGroup     = ['P4_2/mnm']
    NPW = [[192,192,96]]
#    PhaseList      = ['O70'] 
#    InitialL0Guess = [[1.,2,np.around(2*np.sqrt(3),5)]]
#    NThreads       = [ 1]
#    SpaceGroup     = ['Fddd:2']
#    NPW = [[32,32,64]]

#-------------------------------------------------------------------------------#    
# Chain Specifics
    Nref = 1
    # narmtotal = 100
    Nsc_Total = 40
    Nsc_A_min = 2
    Nsc_A_max = 20
    Delta_Nsc_A = 2
    fA_min = 0.1
    fA_max = 0.9
    Delta_fA = 0.01
    chiAB_min =   0.0289
    chiAB_max =   0.0289
    Delta_chiAB = 0.0001
    Neff = 2100
#-------------------------------------------------------------------------------#  
#simulation Specs
    dt = 0.001
    dS = 0.01 # if using CGC
    stress_scale = 0.001
    force_scale = [1.0,0.5]
#-------------------------------------------------------------------------------#  
    backbone_statistics = 'DGC'
    sidearm_statistics_1 = 'DGC'
    sidearm_statistics_2 = 'DGC'

    backbone_species = [1,2]
    sidearm_species_1 = [1]
    sidearm_species_2 = [2]
    # backbone_length = narmtotal-1
    
    sidearm_composition_1 = [1.0]
    sidearm_composition_2 = [1.0]
    
    spacing_1 = 1
    spacing_2 = 1
    
    
    
    IDIR = os.getcwd()
  
    SEED_path = os.path.join(IDIR,"SEEDS")

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
#    Nsc_B_array = np.arange(Nsc_B_min,Nsc_B_max+1e-6,Delta_Nsc_B)
    
    chiAB_array = np.around(np.arange(chiAB_min,chiAB_max+1e-6,Delta_chiAB),5)
#    chiAB_array = np.array([0.0114,0.0124,0.0134,0.0139,0.0144,0.0154])

    ou = open('specialunder.dat','w+')
    oo = open('specialover.dat','w+')
    of = open('specialfluct.dat','w+')
    oa = open('specialall.dat','w+')

    print('f_A eps NscA NscB Ntot N_BB')
    
    Ntot_store_under = dict()
    for chiAB in chiAB_array:
        for Nsc_A in Nsc_A_array:
            print('-------------------------')
            Nsc_B = Nsc_Total-Nsc_A
            nA_check = nA_calc(fA_array,Neff,Nsc_A)
            nB_check = nA_calc(1-fA_array,Neff,(Nsc_Total-Nsc_A))
            nA_int = np.around(nA_check)
            nB_int = np.around(nB_check)
            narmtotal=nA_int+nB_int
            backbone_length = narmtotal-1
            epsilon = np.sqrt((40-Nsc_A+1)/(Nsc_A+1))
            fA_act = fA_calc(nA_int,nB_int,Nsc_A,(Nsc_Total-Nsc_A))
            Ntot = fA_tot(nA_int,nB_int,Nsc_A,(Nsc_Total-Nsc_A))
            # nA_check = nA_calc(fA_array,narmtotal,Nsc_A,Nsc_B)
            # nA_int = np.array(list(set(list(np.around(nA_check))+list([narmtotal]))))
            # del_loc = np.where(nA_int==narmtotal)[0]
            # nA_int = np.delete(nA_int,del_loc)
            # nA_int = np.sort(nA_int)
            # fA_act = fA_calc(nA_int,narmtotal-nA_int,Nsc_A,Nsc_B)
            Ntot_store_under[epsilon] = []
            
            Nmin = np.min(Ntot)
            abratio = 1/(np.sqrt((Nsc_A+1)/(40-Nsc_A+1)))
            for i in range(0,len(fA_act)):
                
                
                # text=f'chiAB_{chiAB:0.4f}/NscA_{Nsc_A}_NscB_{Nsc_B}/fA{fA_act[i]:0.5f}/F0_phases.dat  \n'
                # print(Ntot[i])
                text = (f'{chiAB}   {Nsc_A}  {Nsc_B} {Ntot[i]} {fA_act[i] : 0.5f}  \n')
                if int(Ntot[i])>=2100:
                    oo.write(text)
                if int(Ntot[i])<=2100:
                    print(f'{fA_act[i] :0.5f} {abratio :0.5f} {Nsc_A} {40-Nsc_A} {Ntot[i]} {backbone_length[i]}')

                    ou.write(text)
                    Ntot_store_under[epsilon].append(Ntot[i])
                if abs(int(Ntot[i])-2100)<8:
                    of.write(text)
                oa.write(text)

            # print(np.std(Ntot_store_under[epsilon])/np.mean(Ntot_store_under[epsilon]))
            # print(np.std(Ntot_store_under[epsilon]))

    oa.close()
    ou.close()
    oo.close()
    of.close()
                
                
                
                
                
                
#                 k = 0
#                 for Phase in PhaseList:
                       
#                     if Phase=='LAM' or Phase=='DIS':
#                         d = 1
#                     elif Phase=='ACYL' or Phase=='HEX':
#                         d = 2
#                     else:
#                         d = 3
#                     WDIR=f'chiAB_{chiAB:0.4f}/NscA_{Nsc_A}_NscB_{Nsc_B}/fA{fA_act[i]:0.5f}/{Phase}Phase/'
#                     # print(WDIR)
#                     if os.path.exists(WDIR):
#                         print("{} exists...skipping.".format(WDIR)) 
#                         continue
#                     else:
#                         os.makedirs(WDIR)
                    
                    
#                     if fA_act[i]<0.5:
#                         sidearm_length_1 = Nsc_A
#                         sidearm_length_2 = Nsc_B
#                         sidearm_coverage_1 = nA_int[i]/narmtotal[i]
#                         sidearm_coverage_2 = nB_int[i]/narmtotal[i]
                        
#                         backbone_composition = [sidearm_coverage_1,sidearm_coverage_2]
                
                
#                         backbone_list = [backbone_statistics,backbone_length[i],backbone_species,\
#                                           backbone_composition]
#                         chain_list = [backbone_list]
                
#                         if Nsc_A>0 or Nsc_B>0:
#                             sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
#                                               sidearm_composition_1,sidearm_coverage_1,spacing_1]
                            
#                             sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
#                                               sidearm_composition_2,sidearm_coverage_2,spacing_2]\
                                
#                             sidechain_list = [sidechain_1_list,sidechain_2_list]
#                             chain_list+=sidechain_list
#                     if fA_act[i]>0.5:
#                         sidearm_length_2 = Nsc_A
#                         sidearm_length_1 = Nsc_B
#                         sidearm_coverage_2 = nA_int[i]/narmtotal[i]
#                         sidearm_coverage_1 = nB_int[i]/narmtotal[i]
                        
#                         backbone_composition = [sidearm_coverage_1,sidearm_coverage_2]
                
                
#                         backbone_list = [backbone_statistics,backbone_length[i],backbone_species,\
#                                           backbone_composition]
#                         chain_list = [backbone_list]
                
#                         if Nsc_A>0 or Nsc_B>0:
#                             sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
#                                               sidearm_composition_1,sidearm_coverage_1,spacing_1]
                            
#                             sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
#                                               sidearm_composition_2,sidearm_coverage_2,spacing_2]\
                                
#                             sidechain_list = [sidechain_1_list,sidechain_2_list]
#                             chain_list+=sidechain_list

                            
                            
                        
#                     #pass internal checks
#                     check_equality(len(chain_list)==0,False,"Need more than 0 chains")
#                     for l in range(0,len(chain_list),1):
#                         check_equality(len(chain_list[l])>6,False,"Too many side chain inputs")
#                         check_equality(len(chain_list)!=0,True,"Need more than 0 chains")
#                         list_of_checks = [[len(chain_list[l][2]),len(chain_list[l][3]),'check chain species and composition'],\
#                                             [np.isclose(np.sum(chain_list[l][3]),1),True,'backbone composition does not sum to 1']]
                
#                         for checkchecks in list_of_checks:
#                             check_equality(checkchecks[0],checkchecks[1],checkchecks[2])                 
                                
#                     fieldpath = os.path.join(SEED_path,Phase+'_fields.in')
#                     if len(InitialL0Guess[k])==1:
#                         initial_box = [10/np.sqrt(Nref)]
#                         cell_scaling = InitialL0Guess[k][0]
#                     else:
#                         initial_box = InitialL0Guess[k]
#                         cell_scaling = 10/np.sqrt(Nref)
#                     REPLACE_lIST = [fieldpath,chain_list,[chiAB],d,chain_label,dS,dt,stress_scale,\
#                                     force_scale,cell_scaling,SpaceGroup[k],NPW[k],Nref,initial_box]
#                     make_input(Phase,REPLACE_lIST,WDIR)
                    
#                     if NThreads[k] == 1:
#                         submitfile='SEEDS/submit.sh'
#                     elif NThreads[k] == -1:
# #                        if PhaseList[k] != 'SIGMA':
#                          submitfile='SEEDS/submitGPU.sh'
# #                        else:
# #                            submitfile='SEEDS/submitGPU-P100.sh'
#                     else:
#                         submitfile='SEEDS/submitPLL.sh'
            
                    
                    
#                     with open('%s/submit.sh' % WDIR, 'w') as fout:
#                         with open(submitfile,'r') as f:
#                             for line in f:
#                                 line = line.replace('__JOBNAME__',WDIR)
#                                 line = line.replace('__INFILE__',f'{PhaseList[k]}.in')
#                                 line = line.replace('__OUTFILE__',f'{PhaseList[k]}.out')
#                                 line = line.replace('__PFTPATH__',PolyFTS_Path)
            
#                                 if NThreads[k] > 1:
#                                     line = line.replace('__NTHREADS__',str(NThreads[k]))
#                                 fout.write(line)
#                     k+=1

            

                
#                     os.chdir(WDIR)
# #
#                     import subprocess
# #                    cmd="qsub submit.sh"
#                     cmd="sbatch submit.sh"
#                     subprocess.call(cmd.split())
#                     os.chdir(IDIR)
