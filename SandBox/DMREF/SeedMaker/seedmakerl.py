#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 17:42:31 2020

@author: tquah
"""

from make_input_file import make_input
import numpy as np
import os
import pdb
import shutil
import time
import argparse
from itertools import combinations
import subprocess
import pandas as pd


def species_counter(backbone_species,sidearm_species):
    import copy
    species_count = copy.deepcopy(backbone_species)
    for sidearm in sidearm_species:
        species_count+=sidearm
    return sorted(list(set(species_count)))

def fA_calc(nA,nB,NA,NB):
    return (nA*(NA+1))/((nA*(NA+1))+(nB*(NB+1)))

def nA_calc(fA,nTot,NA,NB):
    return (fA*nTot*(NB+1))/((1-fA)*(NA+1)+fA*(NB+1))


def check_equality(a,b,message):
    assert a==b, message
def check_set_equality(a,b,message):
    assert sorted(a)==sorted(b), message

def listcheck(l,message):
    uniquelist = [",".join(map(str, comb)) for comb in combinations(l, 2)]
    for pair in uniquelist:
        pair_elements = pair.split(',')
        check_equality(pair_elements[0],pair_elements[1],message)

def combination_math(n,r):
    return np.math.factorial(n)/np.math.factorial(r)/np.math.factorial(n-r)

def chain_checking(chain_list):
    #pass internal checks on chain
    check_equality(len(chain_list)==0,False,"Need more than 0 chains")
    for l in range(0,len(chain_list),1):
        check_equality(len(chain_list[l])>6,False,"Too many side chain inputs")
        check_equality(len(chain_list)!=0,True,"Need more than 0 chains")
        list_of_checks = [[len(chain_list[l][2]),len(chain_list[l][3]),'check chain species and composition'],\
                            [np.isclose(np.sum(chain_list[l][3]),1),True,'backbone composition does not sum to 1']]

    for checkchecks in list_of_checks:
        check_equality(checkchecks[0],checkchecks[1],checkchecks[2])                 

def sweep_list_return(lmin,lmax,lstep,dec):
    newlist = []
    checklist = []
    for  i in range(0,len(lmin)):
        if lstep[i]>=0:
            newlist.append(np.round(np.arange(lmin[i],lmax[i]+1e-6,lstep[i]),dec))
        elif lstep[i]<0:
            newlist.append(np.round(np.arange(lmin[i],lmax[i]-1e-6,lstep[i]),dec))
        checklist.append(len(newlist[i]))
    listcheck(checklist,'sweep lists need to be same length')
    return newlist

def sweep_counter(paramdict,name):
    if type(paramdict[name]) == list:
        return len(paramdict[name][0])
    else:
        return len(paramdict[name])
    
def make_WDIR(Phase,countlist,args,itterlist,parameter_dict):
    WDIR = ''
    for r in range(0,len(countlist)):
        if type(parameter_dict[args.directory_structure[r]]) == list:
            WDIR+=(f'{args.directory_structure[r]}{parameter_dict[args.directory_structure[r]][args.directory_nameloc[r]][itterlist[r]]:0.5}')
        else:
            WDIR+=(f'{args.directory_structure[r]}{parameter_dict[args.directory_structure[r]][itterlist[r]]:0.5}')
        if r!=len(countlist):
            WDIR+='/'
    WDIR+=f'{Phase}Phase'
    if os.path.exists(WDIR):
        print("{} exists...skipping.".format(WDIR)) 
        contlogic = True
    else:
        os.makedirs(WDIR)
        contlogic = False

    return WDIR,contlogic

def dimension_determine(Phase):
    if Phase=='LAM' or Phase=='DIS':
        d = 1
    elif Phase=='ACYL' or Phase=='HEX':
        d = 2
    else:
        d = 3
    return d
                        
def chain_maker(args,parameter_dict,itterlist):
    nbb_loc = args.directory_structure.index('nbb')
    chi_loc = args.directory_structure.index('chi')
    nsc_loc = args.directory_structure.index('nsc')
    f_loc = args.directory_structure.index('f')
    
    backbone_composition = []
    
    for r in range(0,len(parameter_dict['f'])):
        backbone_composition.append(parameter_dict['f'][r][itterlist[f_loc]])
    
    #assemble backbone
    backbone_list = [args.chain_stat,parameter_dict['nbb'][itterlist[nbb_loc]],args.backbone_species,\
                  backbone_composition]
    
    chain_list = [backbone_list]
    #assemble sidechains
    full_side_chain_list = []
    for r in range(0,len(parameter_dict['nsc'])):
        if int(parameter_dict['nsc'][r][itterlist[nsc_loc]])>0:
            sidearm_coverage = backbone_composition[r]
            side_chain_list = [args.chain_stat,int(parameter_dict['nsc'][r][itterlist[nsc_loc]]),args.sidearm_species[r],\
                               args.sidearm_composition[r],parameter_dict['f'][r][itterlist[f_loc]],args.sidearm_spacing[r]]
            full_side_chain_list.append(side_chain_list)
    chain_list+=full_side_chain_list
    chain_checking(chain_list)
    return chain_list,nbb_loc,chi_loc,nsc_loc,f_loc

def make_submit(WDIR,args,Phase,q):
    if args.NThreads[q] == 1:
        submitfile='SEEDS/submit.sh'
    elif args.NThreads[q]  == -1:
        if Phase != 'SIGMA':
            submitfile='SEEDS/submitGPU.sh'
        else:
            submitfile='SEEDS/submitGPU-P100.sh'
    else:
        submitfile='SEEDS/submitPLL.sh'
    with open('%s/submit.sh' % WDIR, 'w') as fout:
        with open(submitfile,'r') as f:
            for line in f:
                line = line.replace('__JOBNAME__',WDIR)
                line = line.replace('__INFILE__',f'{Phase}.in')
                line = line.replace('__OUTFILE__',f'{Phase}.out')
                line = line.replace('__PFTPATH__',args.polyFTS_path)

                if args.NThreads[q] > 1:
                    line = line.replace('__NTHREADS__',WDIR)
                fout.write(line)
    fout.close
    f.close()

def domain_size_extractor(content,d):
    list_of_cell = list(filter(lambda x: 'a_1' in x, content))
    if d!=1:
        print('High Dim')
        r = list_of_cell[0]
        Box_initial = abs(float(r[r.find("(")+1:r.find(",")]))
        r = list_of_cell[1]
        Box_change = abs(float(r[r.find("(")+1:r.find(",")]))
        conversion = Box_initial/Box_change
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        Box_Final = np.around(abs(float(r[r.find("(")+1:r.find(",")]))*conversion,5)
    else:
        r = list(filter(lambda x: 'Final simulation cell' in x, content))[0]
        Box_Final = np.around(float(r[r.find("(")+1:r.find(")")]),5)
    return Box_Final


def volume_fraction(content):
    r = list(filter(lambda x: 'volfrac' in x, content))[0]
    return np.around(float(r[r.find("=")+1:r.find(";")]),5)

def stress_tensor_check(content):
    r = list(filter(lambda x: 'L2-norm of stress tensor' in x, content))[-1]
    return float(r[r.find("+")+1:])

    
def spacegroup_lgic(string):
    
    if string=='None':
        return None
    else:
        return string

def directory_nameloc_lgic(string):
    if string=='None':
        return None
    else:
        return int(string)
    
def sidearm_species_lgic(string):
    newlist = string
    
    rlist = []
    for num in newlist.split('-'):
        rlist.append(int(num))
    return rlist
        

def sidearm_composition_lgic(string):
    newlist = string
    
    rlist = []
    for num in newlist.split('-'):
        rlist.append(float(num))
    return rlist

def string_list(string):
    return 

def str_to_bool(s):
    if s=='True' or s=='TRUE' or s=='true':
        return True 
    elif s=='False' or s=='FALSE' or s=='false':
        return False 

def Neff_neutrual_backbone(Nbb,Nsc):
    return Nbb*Nsc

def Neff_block_backbone(Nbb,Nsc):
    return (Nbb)*(Nsc+1)    
    



if __name__ == '__main__':


    IDIR = os.getcwd()

    parser = argparse.ArgumentParser(description='Tool to sweep bottlebrushes')
    parser.add_argument('-PolyFTS', '--polyFTS_path', action='store', default='/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release',help='PolyFts Path',type = str)
    # parser.add_argument('-s', '--sweeptype', action='store', default='chi',help='Sweep Flory-Huggins (chi), armlength (Nsc), or backbone (Nbb) or Volume Fraction (f)',type = str)
    parser.add_argument('-stat', '--chain_stat', action='store', default='DGC',help='Chain Statistics', type = str)
    parser.add_argument('-nt', '--NThreads', action='store',nargs='+', default=[1,1],help='Number of Threads',type = int)
    parser.add_argument('-sp', '--seed_path', action='store', default=os.path.join(IDIR,'SEEDS'),help='Path to Seeds',type = str)
    #Simulation Details
    parser.add_argument('-ds', '--contour_step', action='store', default=0.1,help='Contour step for CGC',type = float)
    parser.add_argument('-cl', '--chain_label', action='store', default='AB-Bottlebrush',help='Chain Label',type = str)
    parser.add_argument('-dt', '--timestep', action='store', default=0.01,help='Timestep size',type = float)
    parser.add_argument('-ss', '--stressscale', action='store', default=0.1,help='Stress Mobility',type = float)
    parser.add_argument('-fs', '--forcescale', action='store',nargs='+', default=[1.0,1.0],help='Force Mobility',type = float)
    parser.add_argument('-L0', '--initial_box_size', action='store',nargs='+', default=[5.0,5.0],help='Initial Box Size',type = float)
    parser.add_argument('-cu', '--cell_updater', action='store', default='Broyden',help='Cell updater usually Broyden or Euler',type = str)
    parser.add_argument('-fu', '--field_updater', action='store', default='EMPEC',help='Field Updater EMPEC or SIS',type = str)
    parser.add_argument('-ftol', '--force_tol', action='store', default=1e-5,help='Force stopping tolerance',type = float)
    parser.add_argument('-stol', '--stress_tol', action='store', default=1e-4,help='stress stopping tolerance',type = float)
    parser.add_argument('-Nrules', '--Neff_Rules', action='store', default=Neff_block_backbone,help='Rules to calculate Neff')

    #Chain Archetecture
    parser.add_argument('-dirs', '--directory_structure', action='store',nargs='+', default=['chi','nbb','nsc','f'],help='Order directory for sweepable parameters',type = str)
    parser.add_argument('-dirnl', '--directory_nameloc', action='store',nargs='+', default=[0,None,0,0],help='Location for directory name location',type = directory_nameloc_lgic)

    parser.add_argument('-i','--instruction',action='store',default = 'SeedOrder.csv',help='Instructions for seedorder',type=str)
    parser.add_argument('-runpoly', '--runpolyfts', action='store', default=False,help='Run polyFTS directly',type = str_to_bool)
    parser.add_argument('-subpoly', '--submitpolyfts', action='store', default=False,help='Submit polyFTS',type = str_to_bool)
    parser.add_argument('-chain', '--chainbool', action='store', default=False,help='Chaining should be done on last variable (directory_structure) cannot use with submit',type = str_to_bool)
    args = parser.parse_args()
 
    print(args)
    if os.path.exists(args.seed_path)!=True:
        raise RuntimeError('Seed Path Does not Exist')
    df = pd.read_csv(args.instruction)
    check_set_equality(args.directory_structure,['chi','nbb','nsc','f'],'Directory Structure needs chi,nbb,nsc,f')
        
    for i in range(0,len(df),1):
        reference_length = df['nref'][i]
        sidechains = df['nsclength'][i]
        if sidechains>0:
            sidearm_composition = [float(i) for i in df['nsccomp'][i].split(';')]
            sidearm_species = [int(i) for i in df['sidearmspecies'][i].split(';')]
        else:
            sidearm_composition=[None]
            sidearm_species = [None]
        
        
        backbones = df['nbblength'][i]
        backbone_species = [int(i) for i in df['backbonespecies'][i].split(';')]
        total_species = species_counter(backbone_species,sidearm_species)
        expected_interaction_parameters = int(combination_math(len(total_species),2))

        Neff = args.Neff_Rules(backbones,sidechains)

        volumefraction = [float(i) for i in df['nbbcomp'][i].split(';')]
        floryhuggins_interaction_parameter_N = df['floryhugginsparametersN'][i].split(';')
        floryhuggins_interaction_parameter = []
        for chiab in floryhuggins_interaction_parameter_N:
            if len(chiab)>0: 
                chiabfloat = float(chiab)
                floryhuggins_interaction_parameter.append(chiabfloat/Neff)
        
        
        number_of_planewaves = df['npw'][i]
        d = dimension_determine(Phase)

    

        
        

    
    # #setup 
    # countlist = []
    # for i in range(0,len(args.directory_structure),1):
    #     countlist.append(sweep_counter(parameter_dict,args.directory_structure[i]))
    # for q in range(0,len(args.phase)):
    
    #     for m in range(0,countlist[0]):
            
    #         for n in range(0,countlist[1]):
            
    #             for o in range(0,countlist[2]):
                
    #                 for p in range(0,countlist[3]):
                    

    #                     Phase = args.phase[q]
    #                     d = dimension_determine(Phase)

    #                     itterlist = [m,n,o,p,q]

    #                     #Section used to make working directory
    #                     WDIR,Contlogic = make_WDIR(Phase,countlist,args,itterlist,parameter_dict)
    #                     if Contlogic:
    #                         if args.chainbool:
    #                             os.chdir(WDIR)
    #                             so = open('STATUS','r')
    #                             status = int(so.read())
    #                             so.close()
    #                             if status==2:
    #                                 print('Simulation Converged!')
    #                                 phaseout = open(f'{Phase}.out')
    #                                 content = phaseout.read().splitlines()
    #                                 phaseout.close()
                                    
    #                                 finalcell = domain_size_extractor(content,d)
                                    
    #                                 args.initial_box_size[q] = finalcell/(10/np.sqrt(args.reference_length))
                                    
    #                                 full_WDIR = os.path.join(IDIR,WDIR)
    #                                 field_path_1 = os.path.join(full_WDIR,'fields_k.bin')
    #                                 field_path_2 = os.path.join(full_WDIR,'fields_k.dat')

    #                                 if os.path.isfile(field_path_1):
    #                                     fieldsin_path = field_path_1
    #                                 elif os.path.isfile(field_path_2):
    #                                     fieldsin_path = field_path_2
    #                                 else:
    #                                     print('field path not found! Using Seedpath')
    #                                     os.chdir(IDIR)
    #                                     break
    #                             os.chdir(IDIR)
    #                         continue
    #                     #fA need to think of a general way to recompute maybe use polyfts to rename the directory
    #                     #set index 
    #                     #chain archetecture
    #                     chain_list,nbb_loc,chi_loc,nsc_loc,f_loc = chain_maker(args,parameter_dict,itterlist)
    #                     fieldpath = 'fields.in'
    #                     chiN = []
    #                     for chi in parameter_dict['chi']:
    #                         chiN.append(chi[itterlist[chi_loc]])
                        
    #                     WDIR_fieldpath = os.path.join(WDIR,fieldpath)
    #                     if p==0 or args.chainbool==False:
    #                         fieldsin_path = os.path.join(args.seed_path,f'{Phase}_fields.in')
    #                         shutil.copy(fieldsin_path,WDIR_fieldpath)
    #                     else:
    #                         shutil.copy(fieldsin_path,WDIR_fieldpath)

                        
    #                     #replace list for make input file
    #                     REPLACE_lIST = [fieldpath,chain_list,chiN,d,args.chain_label,args.contour_step,args.timestep,args.stressscale,\
    #                             args.forcescale,args.initial_box_size[q],[args.number_of_planewaves[q]],args.reference_length,\
    #                                 args.cell_updater,args.force_tol,args.stress_tol,args.field_updater]

    #                     make_input(Phase,REPLACE_lIST,WDIR)
    #                     #make submit file incase you want to run the simulation individually
    #                     make_submit(WDIR,args,Phase,q)
    #                     if args.submitpolyfts:
    #                         os.chdir(WDIR)
    #                         cmd="qsub submit.sh"
    #                         subprocess.call(cmd.split())
    #                         os.chdir(IDIR)
    #                         time.sleep(1)

                        
                        
    #                     #runpolyfts for serial
    #                     if args.runpolyfts:
    #                         if args.NThreads[q]==1: 
    #                             os.chdir(WDIR)

    #                             polyfts_path_serial = os.path.join(args.polyFTS_path,'PolyFTS.x')
    #                             cmd = polyfts_path_serial+f' {Phase}.in'
    #                             f = open(f"{Phase}.out","w")
    #                             subprocess.call(cmd.split(),stdout=f)

    #                         else:
    #                             print('Parallel and GPU not currently supported')
                            
    #                         if args.chainbool:
    #                             so = open('STATUS','r')
    #                             status = int(so.read())
    #                             so.close()
    #                             if status==2:
    #                                 print('Simulation Converged!')
                                    
    #                                 phaseout = open(f'{Phase}.out')
    #                                 content = phaseout.read().splitlines()
    #                                 phaseout.close()
                                    
    #                                 finalcell = domain_size_extractor(content,d)
                                    
    #                                 args.initial_box_size[q] = finalcell/(10/np.sqrt(args.reference_length))
                                    
                                    
    #                                 full_WDIR = os.path.join(IDIR,WDIR)
    #                                 field_path_1 = os.path.join(full_WDIR,'fields_k.bin')
    #                                 field_path_2 = os.path.join(full_WDIR,'fields_k.dat')

    #                                 if os.path.isfile(field_path_1):
    #                                     fieldsin_path = field_path_1
    #                                 elif os.path.isfile(field_path_2):
    #                                     fieldsin_path = field_path_2
    #                             elif status==0 or status==3:
    #                                 print('Simulation was killed or ran out of time')
    #                                 os.chdir(IDIR)
    #                                 break
    #                             else:
    #                                 print('Simulation was divergent')
    #                                 os.chdir(IDIR)
    #                                 break
    #                         os.chdir(IDIR)
