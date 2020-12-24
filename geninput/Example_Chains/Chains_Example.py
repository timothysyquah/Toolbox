#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Chains is a library of functions I use to generate input files for PolyFTS 

Created on Mon Dec  7 14:09:09 2020

@author: tquah
"""
import os
IDIR = os.getcwd()
components = IDIR.split('/')
path_to_tools = '/'+components[1]+'/'+components[2]+'/.timtools'
op = open(path_to_tools,'r')
chainpath =op.read()+'/geninput/'
op.close()

import numpy as np
import sys
sys.path.append(chainpath)

from Chains import Input_Builder,Component_Statistics_Counter,\
    Add_Model_Statistics,Grafting_Determinator,chiN_generator,parameter_species

WDIR = os.getcwd()

#
backbone_statistics = 'DGC'

backbone_length = 99 # This is a quirk from continuous chains 
backbone_species = [1,2]
backbone_composition = [0.5,0.5]
backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                  backbone_composition]

sidearm_statistics_1 = 'DGC'
sidearm_length_1 = 20 #
sidearm_species_1 = [1]
sidearm_composition_1 = [1.0]
sidearm_coverage_1 = 0.5
spacing_1 = 1
sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
                  sidearm_composition_1,sidearm_coverage_1,spacing_1]






chain_list = [backbone_list,sidechain_1_list]





#Si
field = 'fields.in'
chiN = [10.0] 
d = 1
chain_label = 'AB-Bottlebrush'
dt = 0.001
dS = 0.01 # if using CGC if using discrete chains should include it still, but will be not be used
stress_scale = 0.001
force_scale = [1.0,0.5]
cellscale = 10
space_group =  None
npw = [128]
Nref = 1
initial_box = [10/np.sqrt(Nref)]
#Some Adjustable Values
INFILE='LAM.in'
diffuser_method='SOS' #Again not used for discrete chains, but still should be included but will not be used
invzeta=0.1/Nref #compressibility 
force_tol=1e-6
stress_tol=1e-5
kuhn_length=[1.]  
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
field_updater = 'EMPEC'
timesteps_block = 1000
block_num = 100000
variablecell_bool = 'True'
density_history_bool = 'False'
field_history_bool = 'False'
density_chain_bool = 'False'
format_fields_bool = 'False'
volfrac_chain = 1.0
cell_updater = 'Broyden'
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
