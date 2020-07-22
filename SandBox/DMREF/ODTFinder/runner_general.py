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

def sweep_list_return(lmin,lmax,lstep):
    newlist = []
    checklist = []
    for  i in range(0,len(lmin)):
        if lstep[i]>=0:
            newlist.append(np.round(np.arange(lmin[i],lmax[i]+1e-6,lstep[i]),5))
        elif lstep[i]<0:
            newlist.append(np.round(np.arange(lmin[i],lmax[i]-1e-6,lstep[i]),5))
        checklist.append(len(newlist[i]))
    listcheck(checklist,'sweep lists need to be same length')
    return newlist

def sweep_counter(paramdict,name):
    if type(paramdict[name]) == list:
        return len(paramdict[name][0])
    else:
        return len(paramdict[name])
    
def make_WDIR(Phase,countlist,args,itterlist):
    WDIR = ''
    for r in range(0,len(countlist)):
        if type(parameter_dict[args.directory_structure[r]]) == list:
            WDIR+=(f'{args.directory_structure[r]}{parameter_dict[args.directory_structure[r]][args.directory_nameloc[r]][itterlist[r]]:0.3f}')
        else:
            WDIR+=(f'{args.directory_structure[r]}{parameter_dict[args.directory_structure[r]][itterlist[r]]:0.3f}')
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
        
        sidearm_coverage = backbone_composition[r]
        side_chain_list = [args.chain_stat,int(parameter_dict['nsc'][r][itterlist[nsc_loc]]),args.sidearm_species[r],\
                           args.sidearm_composition[r],parameter_dict['f'][r][itterlist[f_loc]],args.sidearm_spacing[r]]
        full_side_chain_list.append(side_chain_list)
    chain_list+=full_side_chain_list
    chain_checking(chain_list)
    return chain_list,nbb_loc,chi_loc,nsc_loc,f_loc

def make_submit(args,Phase):
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

def domain_size_extractor(content):
    list_of_cell = list(filter(lambda x: 'a_1' in x, content))
    if len(list_of_cell)>2:
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

def spacegroup_finder(args):
    so = open('spacegroup.dat','r')
    spacegroup_dat = so.read().splitlines()
    so.close
    spacegroup = []
    for phase in args.phase:
        for line in spacegroup_dat:
            splits = line.split(' ')
            if phase==splits[0]:
                if splits[1]=='None':    
                    spacegroup.append(None)
                else:
                    spacegroup.append(splits[1])
    return spacegroup

if __name__ == '__main__':


    IDIR = os.getcwd()

    parser = argparse.ArgumentParser(description='Tool to sweep bottlebrushes')
    parser.add_argument('-PolyFTS', '--polyFTS_path', action='store', default='/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release',help='PolyFts Path',type = str)
    # parser.add_argument('-s', '--sweeptype', action='store', default='chi',help='Sweep Flory-Huggins (chi), armlength (Nsc), or backbone (Nbb) or Volume Fraction (f)',type = str)
    parser.add_argument('-stat', '--chain_stat', action='store', default='DGC',help='Chain Statistics', type = str)
    parser.add_argument('-p', '--phase', action='store',nargs='+', default=['DIS','LAM'],help='Phases that should be investigated',type = str)
    # parser.add_argument('-sg', '--spacegroup', action='store',nargs='+', default=[None,None],help='Spacegroups for Phases',type = spacegroup_lgic)
    parser.add_argument('-nt', '--NThreads', action='store',nargs='+', default=[1,1],help='Number of Threads',type = int)
    parser.add_argument('-sp', '--seed_path', action='store', default=os.path.join(IDIR,'SEEDS'),help='Path to Seeds',type = str)
    #Simulation Details
    parser.add_argument('-ds', '--contour_step', action='store', default=0.1,help='Contour step for CGC',type = float)
    parser.add_argument('-cl', '--chain_label', action='store', default='AB-Bottlebrush',help='Chain Label',type = str)
    parser.add_argument('-dt', '--timestep', action='store', default=0.01,help='Timestep size',type = float)
    parser.add_argument('-nref', '--reference_length', action='store', default=1,help='Reference Length',type = int)
    parser.add_argument('-ss', '--stressscale', action='store', default=0.1,help='Stress Mobility',type = float)
    parser.add_argument('-fs', '--forcescale', action='store',nargs='+', default=[1.0,1.0],help='Force Mobility',type = float)
    parser.add_argument('-L0', '--initial_box_size', action='store',nargs='+', default=[5.0,5.0],help='Initial Box Size',type = float)
    parser.add_argument('-npw', '--number_of_planewaves', action='store',nargs='+', default=[512,512],help='Number of Planewaves',type = int)
    parser.add_argument('-cu', '--cell_updater', action='store', default='Broyden',help='Cell updater usually Broyden or Euler',type = str)
    parser.add_argument('-fu', '--field_updater', action='store', default='EMPEC',help='Field Updater EMPEC or SIS',type = str)
    parser.add_argument('-ftol', '--force_tol', action='store', default=1e-5,help='Force stopping tolerance',type = float)
    parser.add_argument('-stol', '--stress_tol', action='store', default=1e-4,help='stress stopping tolerance',type = float)
    
    #Chain Archetecture
    parser.add_argument('-bsp', '--backbone_species', action='store',nargs='+', default=[1,2],help='Backbone Species',type=int)
    parser.add_argument('-sas', '--sidearm_species', action='store',nargs='+', default=[[1],[2]],help='Sidearm Species',type=sidearm_species_lgic)
    parser.add_argument('-sac', '--sidearm_composition', action='store',nargs='+', default=[[1.0],[1.0]],help='Sidearm Composition',type=sidearm_composition_lgic)
    parser.add_argument('-space', '--sidearm_spacing', action='store',nargs='+', default=[1,1],help='Sidearm Spacing',type=int)
    #sweepable parameters
    parser.add_argument('-nbbmin', '--backbone_min', action='store', default=99.0,help='Length of Polymer Backbone min',type=float)
    parser.add_argument('-nbbmax', '--backbone_max', action='store', default=99.0,help='Length of Polymer Backbone max',type=float)
    parser.add_argument('-dnbb', '--backbone_step', action='store', default=1.0,help='Length of Polymer Backbone step',type=float)


    parser.add_argument('-nscmin', '--sidechains_min', action='store',nargs='+', default=[20,20],help='Length of Side Chains min',type = float)
    parser.add_argument('-nscmax', '--sidechains_max', action='store',nargs='+', default=[20,20],help='Length of Side Chains max',type = float)
    parser.add_argument('-ndsc', '--sidechains_step', action='store',nargs='+', default=[20,20],help='Length of Side Chains step',type = float)

    parser.add_argument('-chimin', '--floryhuggins_interaction_parameter_min', action='store', default=[0.1],help='Flory-Huggins Interaction Parameter min (12 13 23)',type = float)
    parser.add_argument('-chimax', '--floryhuggins_interaction_parameter_max', action='store', default=[0.1],help='Flory-Huggins Interaction Parameter max (12 13 23)',type = float)
    parser.add_argument('-dchi', '--floryhuggins_interaction_parameter_step', action='store', default=[-0.005],help='Flory-Huggins Interaction Parameter step (12 13 23)',type = float)

    parser.add_argument('-fmin', '--volumefraction_min', action='store',nargs='+', default=[0.5,0.5],help='Volume Fraction min',type = float)
    parser.add_argument('-fmax', '--volumefraction_max', action='store',nargs='+', default=[0.5,0.5],help='Volume Fraction max',type = float)
    parser.add_argument('-df', '--volumefraction_step', action='store',nargs='+', default=[-0.02,0.02],help='Volume Fraction step',type = float)


    parser.add_argument('-dirs', '--directory_structure', action='store',nargs='+', default=['chi','nbb','nsc','f'],help='Order directory for sweepable parameters',type = str)

    parser.add_argument('-dirnl', '--directory_nameloc', action='store',nargs='+', default=[0,None,0,0],help='Location for directory name location',type = directory_nameloc_lgic)
    parser.add_argument('-runpoly', '--runpolyfts', action='store', default=False,help='Run polyFTS directly',type = str_to_bool)
    parser.add_argument('-subpoly', '--submitpolyfts', action='store', default=False,help='Submit polyFTS',type = str_to_bool)
    parser.add_argument('-chain', '--chainbool', action='store', default=False,help='Chaining should be done on last variable (directory_structure) cannot use with submit',type = str_to_bool)

    args = parser.parse_args()
    print(args)
    total_species = species_counter(args.backbone_species,args.sidearm_species)
    #check length of inputs 
    
    spacegroup = spacegroup_finder(args)
    
    phaselogic_check = [len(args.phase),len(spacegroup),len(args.initial_box_size),len(args.number_of_planewaves)]
    listcheck(phaselogic_check,'Check length of npw/phase/spacegroup/initialboxsize')
    
    if type(args.floryhuggins_interaction_parameter_min)==float or type(args.floryhuggins_interaction_parameter_max)==float or type(args.floryhuggins_interaction_parameter_step)==float:
        args.floryhuggins_interaction_parameter_min = [args.floryhuggins_interaction_parameter_min]
        args.floryhuggins_interaction_parameter_max = [args.floryhuggins_interaction_parameter_max]
        args.floryhuggins_interaction_parameter_step = [args.floryhuggins_interaction_parameter_step]

    
    # make sweeplists
    
    #volume fraction
    vollogic_check = [len(args.volumefraction_min),len(args.volumefraction_max),len(args.volumefraction_step),len(args.backbone_species)]
    listcheck(vollogic_check,'Check volume fraction backbone parameters')
    fA_sweep_list = sweep_list_return(args.volumefraction_min,args.volumefraction_max,args.volumefraction_step)

    #flory huggins interaction parameters
    expected_interaction_parameters = int(combination_math(len(total_species),2))
    chilogic_check = [len(args.floryhuggins_interaction_parameter_min),len(args.floryhuggins_interaction_parameter_max),len(args.floryhuggins_interaction_parameter_step),expected_interaction_parameters]
    listcheck(chilogic_check,'Check chi parameters')
    chi_sweep_list = sweep_list_return(args.floryhuggins_interaction_parameter_min,args.floryhuggins_interaction_parameter_max,args.floryhuggins_interaction_parameter_step)

    #number of sidechains sweep
    nsclogic_check = [len(args.sidechains_min),len(args.sidechains_max),len(args.sidechains_step),len(args.sidearm_spacing),len(args.sidearm_composition),len(args.sidearm_species)]
    listcheck(nsclogic_check,'Check sidearm parameters')
    nsc_sweep_list = sweep_list_return(args.sidechains_min,args.sidechains_max,args.sidechains_step)

    
    #number of backbone sweep
    if args.backbone_step>=0:
        nbb_sweep = np.arange(args.backbone_min,args.backbone_max+1e-6,args.backbone_step)
    elif args.backbone_step<0:
        nbb_sweep = np.arange(args.backbone_min,args.backbone_max-1e-6,args.backbone_step)

    #check if seed path exists
    if os.path.exists(args.seed_path)!=True:
        raise RuntimeError('Seed Path Does not Exist')
    print(args.directory_structure)
    check_set_equality(args.directory_structure,['chi','nbb','nsc','f'],'Directory Structure needs chi,nbb,nsc,f')
    
    
    #setup dictionary
    parameter_dict = dict()
    parameter_dict['chi'] = chi_sweep_list
    parameter_dict['nbb'] = nbb_sweep
    parameter_dict['nsc'] = nsc_sweep_list
    parameter_dict['f'] = fA_sweep_list
    

    
    #setup 
    countlist = []
    for i in range(0,len(args.directory_structure),1):
        countlist.append(sweep_counter(parameter_dict,args.directory_structure[i]))
    
    for m in range(0,countlist[0]):
        
        for n in range(0,countlist[1]):
        
            for o in range(0,countlist[2]):
            
                for p in range(0,countlist[3]):
                    
                    for q in range(0,len(args.phase)):
                        Phase = args.phase[q]
                        itterlist = [m,n,o,p,q]

                        #Section used to make working directory
                        WDIR,Contlogic = make_WDIR(Phase,countlist,args,itterlist)    
                        if Contlogic:
                            continue
                        #fA need to think of a general way to recompute maybe use polyfts to rename the directory
                        #set index 
                        #chain archetecture
                        d = dimension_determine(Phase)
                        chain_list,nbb_loc,chi_loc,nsc_loc,f_loc = chain_maker(args,parameter_dict,itterlist)
                        fieldpath = 'fields.in'
                        chiN = []
                        for chi in parameter_dict['chi']:
                            chiN.append(chi[itterlist[chi_loc]])
                        
                        WDIR_fieldpath = os.path.join(WDIR,fieldpath)
                        if p==0 or args.chainbool==False:
                            fieldsin_path = os.path.join(args.seed_path,f'{Phase}_fields.in')
                            shutil.copy(fieldsin_path,WDIR_fieldpath)
                        else:
                            shutil.copy(fieldsin_path,WDIR_fieldpath)

                        
                        #replace list for make input file
                        REPLACE_lIST = [fieldpath,chain_list,chiN,d,args.chain_label,args.contour_step,args.timestep,args.stressscale,\
                                args.forcescale,args.initial_box_size[q],spacegroup[q],[args.number_of_planewaves[q]],args.reference_length,\
                                    args.cell_updater,args.force_tol,args.stress_tol,args.field_updater]

                        make_input(Phase,REPLACE_lIST,WDIR)
                        #make submit file incase you want to run the simulation individually
                        make_submit(args,Phase)
                        if args.submitpolyfts:
                            os.chdir(WDIR)
                            cmd="qsub submit.sh"
                            subprocess.call(cmd.split())
                            os.chdir(IDIR)
                            time.sleep(1)

                        
                        
                        #runpolyfts for serial
                        if args.runpolyfts:
                            if args.NThreads[q]==1: 
                                os.chdir(WDIR)

                                polyfts_path_serial = os.path.join(args.polyFTS_path,'PolyFTS.x')
                                cmd = polyfts_path_serial+f' {Phase}.in'
                                f = open("{Phase}.out","w")
                                subprocess.call(cmd.split(),stdout=f)
                                os.chdir(IDIR)

                            else:
                                print('Parallel and GPU not currently supported')
                            
                            if args.chainbool:
                                so = open('STATUS','r')
                                status = int(so.read())
                                so.close()
                                if status==2:
                                    print('Simulation Converged!')
                                    
                                    phaseout = open(f'{Phase}.out')
                                    content = phaseout.read().splitlines()
                                    phaseout.close()
                                    
                                    finalcell = domain_size_extractor(content)
                                    
                                    args.initial_box_size[q] = finalcell/(10/np.sqrt(args.reference_length))
                                    
                                    
                                    field_path_1 = os.path.join(WDIR,'fields_k.bin')
                                    field_path_2 = os.path.join(WDIR,'fields_k.dat')

                                    if os.path.isfile(field_path_1):
                                        fieldsin_path = field_path_1
                                    elif os.path.isfile(field_path_2):
                                        fieldsin_path = field_path_2
                                elif status==0 or status==3:
                                    print('Simulation was killed or ran out of time')
                                    break
                                else:
                                    print('Simulation was divergent')
                                    break