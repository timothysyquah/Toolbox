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
    spacepath = os.path.join(newsubpath,'spacegroup.dat')
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

    

def make_input(PHASE,REPLACE_DICT,WDIR):
    var_list = list(REPLACE_DICT)
    #Inputs from spawning script
    
    if 'field' in var_list:
       field = REPLACE_DICT['field']
    else:
        field = 'fields.in'

    if 'chain_list' in var_list:
       chain_list = REPLACE_DICT['chain_list']
    else:
        raise Exception('Chain List Required-No Default')

    if 'chi' in var_list:
        chiN = REPLACE_DICT['chi'] 
    else:
        raise Exception('Flory-Huggins Parameter Required!-No Default')
    
    if 'kuhn_length' in var_list:
        kuhn_length = REPLACE_DICT['kuhn_length'] 
    else:
        kuhn_length = [1.]
    
    if 'd' in var_list:
        d = REPLACE_DICT['d'] 
    else:
        raise Exception('Dimension Required!-No Default')
    
    if 'chain_label' in var_list:
        chain_label = REPLACE_DICT['chain_label'] 
    else:
        chain_label = 'AB-BlockPolymer'
    
    if 'dS' in var_list:
        dS = REPLACE_DICT['dS'] 
    else:
        raise Exception('dS Required!-No Default (for continous chains just put a dummy value ie. 0.1)')

    if 'dt' in var_list:
        dt = REPLACE_DICT['dt'] 
    else:
        raise Exception('dt Required!-No Default')

    if 'stress_scale' in var_list:
        stress_scale = REPLACE_DICT['stress_scale'] 
    else:
        stress_scale = 0.001
        
    if 'force_scale' in var_list:
        force_scale = REPLACE_DICT['force_scale'] 
    else:
        force_scale = [1.,1.]

    if 'cellscale' in var_list:
        cellscale = REPLACE_DICT['cellscale'] 
    else:
        raise Exception('cellscale Required!-No Default')

#    space_group =  REPLACE_lIST[10]
    if 'npw' in var_list:
        npw = REPLACE_DICT['npw'] 
    else:
        raise Exception('npw Required!-No Default')

    if 'Nref' in var_list:
        Nref = REPLACE_DICT['Nref'] 
    else:
        Nref = 1

    if 'force_tol' in var_list:
        force_tol = REPLACE_DICT['force_tol'] 
    else:
        force_tol = 1e-5

    if 'stress_tol' in var_list:
        stress_tol = REPLACE_DICT['stress_tol'] 
    else:
        stress_tol = 1e-5
    
    if 'cell_updater' in var_list:
       cell_updater = REPLACE_DICT['cell_updater']
    else:
        cell_updater = 'Broyden'
    
    if 'field_updater' in var_list:
       field_updater = REPLACE_DICT['field_updater']
    else:
        field_updater = 'EMPEC'

    if 'diffuser_method' in var_list:
       diffuser_method = REPLACE_DICT['diffuser_method']
    else:
        diffuser_method = 'SOS'

    space_group = spacegroup_finder(PHASE)
    
    if 'nThreads' in var_list:
       nThreads = REPLACE_DICT['nThreads']
    else:
        raise Exception('nThreads Required!-No Default')


    if nThreads!=1:
        if 'parallel_cuda' in var_list:
            parallel_cuda = REPLACE_DICT['parallel_cuda']
        else:
            raise Exception('parallel_cuda Required!-No Default')
        if 'cuda_thread_block_size' in var_list:
            cuda_thread_block_size = REPLACE_DICT['cuda_thread_block_size']
        else:
            cuda_thread_block_size = 128
    else:
        parallel_cuda = 0
        cuda_thread_block_size = 128

    if 'invzeta' in var_list:
       invzeta = REPLACE_DICT['invzeta']
    else:
        invzeta = 0.1/Nref

    if 'ends' in var_list:
       ends = REPLACE_DICT['ends']
    else:
        ends = True

    if 'add_phase' in var_list:
       add_phase = REPLACE_DICT['add_phase']
    else:
        add_phase = True

    if 'add_phase' in var_list:
       add_phase = REPLACE_DICT['add_phase']
    else:
        add_phase = True


    #Some Adjustable Values
    INFILE=f'{PHASE}.in'
    initial_box = [10/np.sqrt(Nref)]
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
