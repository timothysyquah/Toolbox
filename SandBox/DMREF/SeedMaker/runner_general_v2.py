#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 19 17:42:31 2020

@author: tquah
"""


from make_input_file_dict import make_input
from checkphases import field_checker
import numpy as np
import os
# import pdb
import shutil
import time
import argparse
from itertools import combinations
import subprocess
import re




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
    
# def make_WDIR(Phase,countlist,args,itterlist):
#     WDIR = ''
#     for r in range(0,len(countlist)):
#         if type(parameter_dict[args.itterative_structure[r]]) == list:
#             WDIR+=(f'{args.itterative_structure[r]}{parameter_dict[args.itterative_structure[r]][args.directory_nameloc[r]][itterlist[r]]:0.5}')
#         else:
#             WDIR+=(f'{args.itterative_structure[r]}{parameter_dict[args.itterative_structure[r]][itterlist[r]]:0.5}')
#         if r!=len(countlist):
#             WDIR+='/'
#     WDIR+=f'{Phase}Phase'
#     if os.path.exists(WDIR):
#         print("{} exists...skipping.".format(WDIR)) 
#         contlogic = True
#     else:
#         # os.makedirs(WDIR)
#         contlogic = False

#     return WDIR,contlogic

def make_WDIR(Phase,parameter_dict,args,itterlist):
    WDIR = ''
    for i in range(0,len(args.directory_struct_names),1):
        dir_objects = args.directory_struct_names[i]
        name = ''
        space = ''
        if len(dir_objects[0])>1:
            space='_'
        
        
        for j in range(0,len(dir_objects[0]),1):
            # print(dir_objects[0])
            index_dir_name = args.itterative_structure.index(dir_objects[1])
            if type(parameter_dict[dir_objects[1]])==list:
                if dir_objects[1] == 'f' or dir_objects[1] == 'chi': 
                    name+=f'{dir_objects[0][j]}{parameter_dict[dir_objects[1]][args.nameloc[i][j]][itterlist[index_dir_name]]:0.5f}'
                else:
                    name+=f'{dir_objects[0][j]}{parameter_dict[dir_objects[1]][args.nameloc[i][j]][itterlist[index_dir_name]]:0.1f}'

            else:
                name+=f'{dir_objects[0][j]}{parameter_dict[dir_objects[1]][itterlist[index_dir_name]]:0.1f}'
            if j!=len(dir_objects[0])-1:
                name+=space
        WDIR+=name
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

def extract_loc(args):
    nbb_loc = args.itterative_structure.index('nbb')
    chi_loc = args.itterative_structure.index('chi')
    nref_loc = args.itterative_structure.index('nref')
    nsc_loc = args.itterative_structure.index('nsc')
    f_loc = args.itterative_structure.index('f')
    return nbb_loc,chi_loc,nsc_loc,f_loc,nref_loc


def chain_maker(args,parameter_dict,itterlist):
    nbb_loc = args.itterative_structure.index('nbb')
    chi_loc = args.itterative_structure.index('chi')
    nref_loc = args.itterative_structure.index('nref')
    nsc_loc = args.itterative_structure.index('nsc')
    f_loc = args.itterative_structure.index('f')
    
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
    return chain_list

def make_submit(args,Phase,WDIR):
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
        # print('High Dim')
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

def domain_size_extractor_O70(content,d):
    list_of_cell = list(filter(lambda x: 'a_1' in x, content))
    # print('High Dim')
    r = list_of_cell[0]
    Box_initial = abs(float(r[r.find("(")+1:r.find(",")]))
    r = list_of_cell[1]
    Box_change = abs(float(r[r.find("(")+1:r.find(",")]))
    conversion = Box_initial/Box_change
    rindex = content.index(list(filter(lambda x: 'Final simulation cell' in x, content))[0])
    domain = []
    scinot = re.compile('[+\-]?(?:0|[1-9]\d*)(?:\.\d*)?(?:[eE][+\-]?\d+)')
    for i in range(0,3,1):
        
        process1 = content[rindex+i].replace(',',' ')
        process2 = process1.replace('(',' ') 
        process3 = process2.replace(')',' ') 
    
        r = re.findall(scinot, process3)
        for j in range(0,len(r),1):
            if float(r[j])!=0:
                domain.append(float(r[j]))
    if len(domain)==3:    
        return domain
    else:
        print('error domain size extractor')

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
    newlist = string

    if string=='None':
        return [None]
    else:
        rlist = []
        for num in newlist.split('-'):
            rlist.append(int(num))
        return rlist
    
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

def dir_nametype_lgic(string):
    name_type = string.split('/')
    add_split = name_type[0].split('-')
    new_str = []
    new_str.append(add_split)
    new_str.append(name_type[1])
    return new_str     

def str_to_bool(s):
    if s=='True' or s=='TRUE' or s=='true':
        return True 
    elif s=='False' or s=='FALSE' or s=='false':
        return False 

def gpu_lgic(string):
    if string=='None':
        return None
    else:
        return int(string)


if __name__ == '__main__':


    IDIR = os.getcwd()

    parser = argparse.ArgumentParser(description='Tool to sweep bottlebrushes')
    parser.add_argument('-PolyFTS', '--polyFTS_path', action='store', default='/home/tquah/PolyFTS_ALL/PolyFTS/bin/Release',help='PolyFts Path',type = str)
    # parser.add_argument('-s', '--sweeptype', action='store', default='chi',help='Sweep Flory-Huggins (chi), armlength (Nsc), or backbone (Nbb) or Volume Fraction (f)',type = str)
    parser.add_argument('-stat', '--chain_stat', action='store', default='DGC',help='Chain Statistics', type = str)
    parser.add_argument('-p', '--phase', action='store',nargs='+', default=['DIS','LAM'],help='Phases that should be investigated',type = str)
    parser.add_argument('-nt', '--NThreads', action='store',nargs='+', default=[1,1],help='Number of Threads',type = int)
    parser.add_argument('-sp', '--seed_path', action='store', default=os.path.join(IDIR,'SEEDS'),help='Path to Seeds',type = str)
    parser.add_argument('-spt', '--seed_path_type', action='store', default='Main',help='Is seed in a main directory (main) or job directory (job)',type = str)
    
    #Simulation Details
    parser.add_argument('-ds', '--contour_step', action='store', default=0.1,help='Contour step for CGC',type = float)
    parser.add_argument('-cl', '--chain_label', action='store', default='AB-Bottlebrush',help='Chain Label',type = str)
    parser.add_argument('-dt', '--timestep', action='store', default=0.01,help='Timestep size',type = float)
    
    
    parser.add_argument('-ss', '--stressscale', action='store', default=0.1,help='Stress Mobility',type = float)
    parser.add_argument('-fs', '--forcescale', action='store',nargs='+', default=[1.0,1.0],help='Force Mobility',type = float)
    parser.add_argument('-L0', '--initial_box_size', action='store',nargs='+', default=[[5.0],[5.0]],help='Initial Box Size',type = sidearm_composition_lgic)
    parser.add_argument('-npw', '--number_of_planewaves', action='store',nargs='+', default=[[512],[512]],help='Number of Planewaves',type = sidearm_species_lgic)
    parser.add_argument('-cu', '--cell_updater', action='store', default='Broyden',help='Cell updater usually Broyden or Euler',type = str)
    parser.add_argument('-fu', '--field_updater', action='store', default='EMPEC',help='Field Updater EMPEC or SIS',type = str)
    parser.add_argument('-ftol', '--force_tol', action='store', default=1e-5,help='Force stopping tolerance',type = float)
    parser.add_argument('-stol', '--stress_tol', action='store', default=1e-4,help='stress stopping tolerance',type = float)
    
    #Chain Archetecture
    parser.add_argument('-bsp', '--backbone_species', action='store',nargs='+', default=[1,2],help='Backbone Species',type=int)
    parser.add_argument('-sas', '--sidearm_species', action='store',nargs='+', default=[[1],[2]],help='Sidearm Species',type=sidearm_species_lgic)
    parser.add_argument('-sac', '--sidearm_composition', action='store',nargs='+', default=[[1.0],[1.0]],help='Sidearm Composition',type=sidearm_composition_lgic)
    parser.add_argument('-space', '--sidearm_spacing', action='store',nargs='+', default=[1,1],help='Sidearm Spacing',type=int)
    parser.add_argument('-b', '--kuhnlength', action='store',nargs='+', default=[1.0,1.0],help='Kuhn Length',type=float)

    #sweepable parameters
    parser.add_argument('-nbbmin', '--backbone_min', action='store', default=99.0,help='Length of Polymer Backbone min',type=float)
    parser.add_argument('-nbbmax', '--backbone_max', action='store', default=99.0,help='Length of Polymer Backbone max',type=float)
    parser.add_argument('-dnbb', '--backbone_step', action='store', default=1.0,help='Length of Polymer Backbone step',type=float)


    parser.add_argument('-nscmin', '--sidechains_min', action='store',nargs='+', default=[20,20],help='Length of Side Chains min',type = float)
    parser.add_argument('-nscmax', '--sidechains_max', action='store',nargs='+', default=[10,30],help='Length of Side Chains max',type = float)
    parser.add_argument('-ndsc', '--sidechains_step', action='store',nargs='+', default=[-2,2],help='Length of Side Chains step',type = float)

    parser.add_argument('-chimin', '--floryhuggins_interaction_parameter_min', action='store', default=[0.0289],help='Flory-Huggins Interaction Parameter min (12 13 23)',type = float)
    parser.add_argument('-chimax', '--floryhuggins_interaction_parameter_max', action='store', default=[0.0289],help='Flory-Huggins Interaction Parameter max (12 13 23)',type = float)
    parser.add_argument('-dchi', '--floryhuggins_interaction_parameter_step', action='store', default=[-0.005],help='Flory-Huggins Interaction Parameter step (12 13 23)',type = float)

    parser.add_argument('-fmin', '--volumefraction_min', action='store',nargs='+', default=[0.1,0.9],help='Volume Fraction min',type = float)
    parser.add_argument('-fmax', '--volumefraction_max', action='store',nargs='+', default=[0.9,0.1],help='Volume Fraction max',type = float)
    parser.add_argument('-df', '--volumefraction_step', action='store',nargs='+', default=[0.1,-0.1],help='Volume Fraction step',type = float)

    parser.add_argument('-nref_list', '--reference_length_list', action='store', default=[1],help='Reference Length list',type = float)


    parser.add_argument('-is', '--itterative_structure', action='store',nargs='+', default=['chi','nbb','nsc','nref','f'],help='Order directory for sweepable parameters',type = str)
    parser.add_argument('-drst', '--directory_struct_names', action='store',nargs='+', default=[[['nref'],'nref'],[['NscA_','NscB_'],'nsc'],[['fA'],'f']],help='Order directory for sweepable parameters',type = dir_nametype_lgic)
    parser.add_argument('-dirnl', '--nameloc', action='store',nargs='+', default=[[None],[0,1],[1]],help='Location for directory name location',type = directory_nameloc_lgic)
    parser.add_argument('-runpoly', '--runpolyfts', action='store', default=False,help='Run polyFTS directly',type = str_to_bool)
    parser.add_argument('-subpoly', '--submitpolyfts', action='store', default=False,help='Submit polyFTS',type = str_to_bool)
    parser.add_argument('-chain', '--chainbool', action='store', default=False,help='Chaining should be done on last variable (itterative_structure) cannot use with submit',type = str_to_bool)
    parser.add_argument('-gd', '--gpu_device', action='store', default=None,help='GPU device',type = gpu_lgic)
    # parser.add_argument('-dirnl', '--directory_nameloc', action='store',nargs='+', default=[0,None,0,0,None],help='Location for directory name location',type = directory_nameloc_lgic)
    args = parser.parse_args()
    print(args)
    total_species = species_counter(args.backbone_species,args.sidearm_species)
    
    phaselogic_check = [len(args.phase),len(args.initial_box_size),len(args.number_of_planewaves)]
    listcheck(phaselogic_check,'Check length of npw/phase/spacegroup/initialboxsize')
    
    if type(args.floryhuggins_interaction_parameter_min)==float or type(args.floryhuggins_interaction_parameter_max)==float or type(args.floryhuggins_interaction_parameter_step)==float:
        args.floryhuggins_interaction_parameter_min = [args.floryhuggins_interaction_parameter_min]
        args.floryhuggins_interaction_parameter_max = [args.floryhuggins_interaction_parameter_max]
        args.floryhuggins_interaction_parameter_step = [args.floryhuggins_interaction_parameter_step]

    
    # make sweeplists
    
    #volume fraction
    vollogic_check = [len(args.volumefraction_min),len(args.volumefraction_max),len(args.volumefraction_step),len(args.backbone_species)]
    listcheck(vollogic_check,'Check volume fraction backbone parameters')
    fA_sweep_list = sweep_list_return(args.volumefraction_min,args.volumefraction_max,args.volumefraction_step,dec=5)

    #flory huggins interaction parameters
    expected_interaction_parameters = int(combination_math(len(total_species),2))
    chilogic_check = [len(args.floryhuggins_interaction_parameter_min),len(args.floryhuggins_interaction_parameter_max),len(args.floryhuggins_interaction_parameter_step),expected_interaction_parameters]
    listcheck(chilogic_check,'Check chi parameters')
    chi_sweep_list = sweep_list_return(args.floryhuggins_interaction_parameter_min,args.floryhuggins_interaction_parameter_max,args.floryhuggins_interaction_parameter_step,dec=8)

    #number of sidechains sweep
    nsclogic_check = [len(args.sidechains_min),len(args.sidechains_max),len(args.sidechains_step),len(args.sidearm_spacing),len(args.sidearm_composition),len(args.sidearm_species)]
    listcheck(nsclogic_check,'Check sidearm parameters')
    nsc_sweep_list = sweep_list_return(args.sidechains_min,args.sidechains_max,args.sidechains_step,dec=1)

    
    #number of backbone sweep
    if args.backbone_step>=0:
        nbb_sweep = np.arange(args.backbone_min,args.backbone_max+1e-6,args.backbone_step)
    elif args.backbone_step<0:
        nbb_sweep = np.arange(args.backbone_min,args.backbone_max-1e-6,args.backbone_step)

    #check if seed path exists
    if os.path.exists(args.seed_path)!=True and args.seed_path_type=='Main':
        raise RuntimeError('Seed Path Does not Exist')
    check_set_equality(args.itterative_structure,['chi','nbb','nsc','f','nref'],'Directory Structure needs chi,nbb,nsc,f,nref')
    
    #check directory logic
    dirlogic_check = [len(args.directory_struct_names),len(args.nameloc)]
    listcheck(dirlogic_check,'Check Directory Naming Valriables')

    
    
    #setup dictionary
    parameter_dict = dict()
    parameter_dict['chi'] = chi_sweep_list
    parameter_dict['nbb'] = nbb_sweep
    parameter_dict['nsc'] = nsc_sweep_list
    parameter_dict['f'] = fA_sweep_list
    parameter_dict['nref'] = np.array(args.reference_length_list)

    fieldpath = 'fields.in'

    
    #setup 
    countlist = []
    for i in range(0,len(args.itterative_structure),1):
        countlist.append(sweep_counter(parameter_dict,args.itterative_structure[i]))
    for q in range(0,len(args.phase)):
    
        for m in range(0,countlist[0]):
            
            for n in range(0,countlist[1]):
            
                for o in range(0,countlist[2]):
                
                    for p in range(0,countlist[3]):
                        
                        for r in range(0,countlist[4]):

                            Phase = args.phase[q]
                            d = dimension_determine(Phase)
                            
                            itterlist = [m,n,o,p,r]
                            WDIR,Contlogic = make_WDIR(Phase,parameter_dict,args,itterlist)

                            nbb_loc,chi_loc,nsc_loc,f_loc,nref_loc = extract_loc(args)
                            
                            
                            #Section used to make working directory
                            # WDIR,Contlogic = make_WDIR(Phase,countlist,args,itterlist)    
                            if Contlogic:
                                if args.chainbool:
                                    os.chdir(WDIR)
                                    so = open('STATUS','r')
                                    status = int(so.read())
                                    so.close()
                                    if status==2:
                                        print(WDIR)
                                        print('Simulation Converged!')
                                        phaseout = open(f'{Phase}.out')
                                        content = phaseout.read().splitlines()
                                        phaseout.close()
                                        
                                        if Phase.find('O70')==0:
                                            
                                            args.initial_box_size[q] = domain_size_extractor(content,d)

                                        else:
                                            finalcell = domain_size_extractor(content,d)
                                            args.initial_box_size[q][0] = finalcell/(10/np.sqrt(args.reference_length_list[itterlist[nref_loc]]))
                                        
                                        full_WDIR = os.path.join(IDIR,WDIR)
                                        field_path_1 = os.path.join(full_WDIR,'fields_k.bin')
                                        field_path_2 = os.path.join(full_WDIR,'fields_k.dat')
                                        field_path_3 = os.path.join(full_WDIR,fieldpath)

                                        if os.path.isfile(field_path_1):
                                            fieldsin_path = field_path_1
                                            field_checker(field_path_1,field_path_3)
                                        elif os.path.isfile(field_path_2):
                                            fieldsin_path = field_path_2
                                            field_checker(field_path_1,field_path_3)
                                        else:
                                            print(WDIR)
                                            print('field path not found! Using Seedpath')
                                            os.chdir(IDIR)
                                            break
                                    os.chdir(IDIR)
                                continue
                            #fA need to think of a general way to recompute maybe use polyfts to rename the directory
                            #set index 
                            #chain archetecture
                            chain_list = chain_maker(args,parameter_dict,itterlist)
                            chiN = []
                            for chi in parameter_dict['chi']:
                                chiN.append(chi[itterlist[chi_loc]])
                            
                            WDIR_fieldpath = os.path.join(WDIR,fieldpath)
                            if p==0 or args.chainbool==False:
                                if args.seed_path_type=='Main':
                                    fieldsin_path = os.path.join(args.seed_path,f'{Phase}_fields.in')
                                else:
                                    field_path_1 = os.path.join(args.seed_path,'fields_k.bin')
                                    field_path_2 = os.path.join(args.seed_path,'fields_k.dat')
                                    if os.path.isfile(field_path_1):
                                        fieldsin_path = field_path_1
                                    elif os.path.isfile(field_path_2):
                                        fieldsin_path = field_path_2
                                shutil.copy(fieldsin_path,WDIR_fieldpath)
                            else:
                                shutil.copy(fieldsin_path,WDIR_fieldpath)
    
                            
                            #replace list for make input file
                            REPLACE_DICT = dict()
                            REPLACE_DICT['field'] = fieldpath
                            REPLACE_DICT['chain_list'] = chain_list
                            REPLACE_DICT['chi'] = chiN
                            REPLACE_DICT['d'] = d
                            REPLACE_DICT['chain_label'] = args.chain_label
                            REPLACE_DICT['dS'] = args.contour_step
                            REPLACE_DICT['dt'] = args.timestep
                            REPLACE_DICT['stressscale'] = args.stressscale
                            REPLACE_DICT['forcescale'] = args.forcescale
                            REPLACE_DICT['chain_list'] = chain_list
                            REPLACE_DICT['npw'] = args.number_of_planewaves[q]
                            REPLACE_DICT['Nref'] = args.reference_length_list[itterlist[nref_loc]]
                            REPLACE_DICT['cell_updater'] = args.cell_updater
                            REPLACE_DICT['force_tol'] = args.force_tol
                            REPLACE_DICT['stress_tol'] = args.stress_tol
                            REPLACE_DICT['field_updater'] = args.field_updater
                            REPLACE_DICT['kuhn_length'] = args.kuhnlength
                            REPLACE_DICT['nThreads'] = args.NThreads[q]
                            if args.NThreads[q]!=1:
                                if args.gpu_device!=None:
                                    REPLACE_DICT['parallel_cuda'] = args.gpu_device
                                else:
                                    raise Exception('GPU device required for Ntreads not equal to 1')
                            if Phase.find('O70')==0:
                                REPLACE_DICT['initial_box'] = args.initial_box_size[q]
                                REPLACE_DICT['cellscale'] = (10/np.sqrt(args.reference_length_list[itterlist[nref_loc]]))


                            else:
                                REPLACE_DICT['cellscale'] = args.initial_box_size[q]

                            # REPLACE_lIST = [fieldpath,chain_list,chiN,d,args.chain_label,args.contour_step,args.timestep,args.stressscale,\
                            #         args.forcescale,args.initial_box_size[q],[args.number_of_planewaves[q]],args.reference_length,\
                            #             args.cell_updater,args.force_tol,args.stress_tol,args.field_updater]
    
    
                            make_input(Phase,REPLACE_DICT,WDIR)
                            #make submit file incase you want to run the simulation individually
                            make_submit(args,Phase,WDIR)
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
                                    f = open(f"{Phase}.out","w")
                                    subprocess.call(cmd.split(),stdout=f)
    
                                else:
                                    os.chdir(WDIR)
                                    polyfts_path_gpu = os.path.join(args.polyFTS_path,'PolyFTSGPU.x')
                                    cmd = polyfts_path_gpu+f' {Phase}.in'
                                    f = open(f"{Phase}.out","w")
                                    subprocess.call(cmd.split(),stdout=f)

                                if args.chainbool:
                                    so = open('STATUS','r')
                                    status = int(so.read())
                                    so.close()
                                    print(WDIR)
                                    if status==2:
                                        print('Simulation Converged!')
                                        
                                        phaseout = open(f'{Phase}.out')
                                        content = phaseout.read().splitlines()
                                        phaseout.close()
                                        
                                        if Phase.find('O70')==0:
                                            
                                            celldim = (10/np.sqrt(args.reference_length_list[itterlist[nref_loc]]))
                                            args.initial_box_size[q] = domain_size_extractor(content,d)

                                        else:
                                            finalcell = domain_size_extractor(content,d)
                                            args.initial_box_size[q][0] = finalcell/(10/np.sqrt(args.reference_length_list[itterlist[nref_loc]]))
                                        
                                        
                                        full_WDIR = os.path.join(IDIR,WDIR)
                                        field_path_1 = os.path.join(full_WDIR,'fields_k.bin')
                                        field_path_2 = os.path.join(full_WDIR,'fields_k.dat')
                                        field_path_3 = os.path.join(full_WDIR,fieldpath)
                                        if os.path.isfile(field_path_1):
                                            fieldsin_path = field_path_1
                                            field_checker(field_path_1,field_path_3)
                                        elif os.path.isfile(field_path_2):
                                            fieldsin_path = field_path_2
                                            field_checker(field_path_1,field_path_3)
                                    elif status==0 or status==3:
                                        print('Simulation was killed or ran out of time')
                                        os.chdir(IDIR)
                                        break
                                    else:
                                        print('Simulation was divergent')
                                        os.chdir(IDIR)
                                        break
                                os.chdir(IDIR)
