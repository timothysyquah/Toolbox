#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 18 12:36:54 2020

@author: tquah
"""
import os
IDIR = os.getcwd()
components = IDIR.split('/')
path_to_tools = '/'+components[1]+'/'+components[2]+'/.timtools'
op = open(path_to_tools,'r')
toolpath = op.read()
chainpath =toolpath+'/geninput/'
newsubpath = toolpath+'/newsubmit/'
op.close()
import numpy as np
import sys
sys.path.append(chainpath)

from Chains import Input_Builder,Component_Statistics_Counter,\
    Add_Model_Statistics,Grafting_Determinator,chiN_generator,parameter_species

def spacegroup_finder(Phase):
    spacepath = os.path.join(newsubpath,'spacegroup.req')
    so = open(spacepath,'r')
    spacegroup_dat = so.read().splitlines()
    # print(spacegroup_dat)
    so.close
    for line in spacegroup_dat:
        splits = line.split(' ')
        if Phase==splits[0]:
            if splits[1]=='None':    
                return None
            else:
                    return splits[1]

    

def make_input(PHASE,REPLACE_lIST,WDIR):
    #Inputs from spawning script
    field = REPLACE_lIST[0]
    chain_list = REPLACE_lIST[1] 
    chiN = REPLACE_lIST[2] 
    d = REPLACE_lIST[3]
    chain_label = REPLACE_lIST[4]
    dS = REPLACE_lIST[5]
    dt = REPLACE_lIST[6]
    stress_scale = REPLACE_lIST[7]
    force_scale = REPLACE_lIST[8]
    cellscale = REPLACE_lIST[9]
#    space_group =  REPLACE_lIST[10]
    npw = REPLACE_lIST[10]
    Nref = REPLACE_lIST[11]
    cell_updater = REPLACE_lIST[12]
    force_tol=REPLACE_lIST[13]
    stress_tol=REPLACE_lIST[14]
    field_updater = REPLACE_lIST[15]
    space_group = spacegroup_finder(PHASE)
    
    
    kuhn_length=[1.]    
    invzeta=0.1/Nref

    #Some Adjustable Values
    INFILE=f'{PHASE}.in'
    diffuser_method='SOS'
    initial_box = [10/np.sqrt(Nref)]
    ends = True
    add_phase = True
    non_primitive_centering = 'True'
    symmetrize = 'on'

    if space_group is None:
        add_phase = False
        initial_box = [cellscale*initial_box[0]]

    inputfileversion = 3
    nummodel = 1
    supported_statistics = ['CGC','DGC','FJC']
    partition_function = 'canonical'
    ham_bool = 'True'
    stress_bool = 'True'
    chem_pot_bool = 'False'
    idealgas_bool = 'False'
    field_type = 'HFields'
    timesteps_block = 1000
    block_num = 1000
    variablecell_bool = 'True'
    density_history_bool = 'False'
    field_history_bool = 'False'
    density_chain_bool = 'False'
    format_fields_bool = 'False'
    volfrac_chain = 1.0
    nThreads = 1
    parallel_cuda = 1
    cuda_thread_block_size = 128


    input_file_path = os.path.join(WDIR,INFILE)
    statistics_models_add = [[f'contourds = {dS : 0.5f}',\
                         f'diffusermethod = {diffuser_method}'],\
                        [f'polymerReferenceN = {Nref}.'],\
                        [f'polymerReferenceN = {Nref}.']]
    
    statistics_interactions_add = [[],\
                                   [f'compressibility_invzetaN = {invzeta}'],
                                   [f'compressibility_invzetaN = {invzeta}']]
    
    components,statistics_list,n_species = Component_Statistics_Counter(chain_list)

    if n_species==2:
        modeltype = 'BlockPolymerMelt2Spec'
    if n_species>=3:
        modeltype = 'BlockPolymerMelt'

    
    models_add,interation_add,n_sidearm_types =\
        Add_Model_Statistics(statistics_list,chain_list,\
                                 supported_statistics,statistics_models_add,\
                                 statistics_interactions_add)
    chain_list,chain_text_list = Grafting_Determinator(chain_list,supported_statistics,ends)
    chiN_list =chiN_generator(chiN,components)
    kuhn_length_text = parameter_species(kuhn_length,n_species)
    force_scale_text = parameter_species(force_scale,n_species)
    initial_box_text = parameter_species(initial_box,d)
    npw_text = parameter_species(npw,d)

    f = open(input_file_path,'w+')    
    Input_Builder(f,inputfileversion,nummodel,modeltype,n_species,\
                     kuhn_length_text,chain_label,n_sidearm_types,\
                     models_add,chain_text_list,chiN_list,interation_add,\
                     d,initial_box_text,npw_text,partition_function,\
                     ham_bool,stress_bool,chem_pot_bool,idealgas_bool,\
                     field_type,field,field_updater,cell_updater,timesteps_block,\
                     block_num,dt,force_scale_text,stress_scale,force_tol,stress_tol,\
                     variablecell_bool,density_history_bool,field_history_bool,\
                     density_chain_bool,format_fields_bool,volfrac_chain,cellscale,
                     add_phase,space_group,non_primitive_centering,\
                     symmetrize,parallel_cuda,cuda_thread_block_size,nThreads).Write_to_Text()
    f.close()
