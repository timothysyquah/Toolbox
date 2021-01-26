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
chainpath =op.read()+'/geninput_update/'
op.close()

import numpy as np
import sys
sys.path.append(chainpath)
print(chainpath)
from Chains import Input_Builder,Component_Statistics_Counter, Chain_Builder,\
    Add_Model_Statistics,Grafting_Determinator,chiN_generator,parameter_species

##Script works on DGC/FJC/CGC
##Script may have issues with CGC Nref not equal to 1 due to lengths and other things...




WDIR = os.getcwd()
backbone_statistics = 'FJC'
#Nbb = 100
backbone_length = 99 # This is a quirk from continuous chains where the N_bb = L+1 
backbone_species = [1,2]
backbone_composition = [0.5,0.5]
backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                  backbone_composition]

sidearm_statistics_1 = 'FJC'
#Nsc = 20
sidearm_length_1 = 19 # This is a quirk from continuous chains where the N_bb = L+1 
sidearm_species_1 = [3]
sidearm_composition_1 = [1.0]
sidearm_coverage_1 = 0.5
graft_density = 1 #currently this doesn't do anything...currently only does graft density of 1
sidechain_1_list = [sidearm_statistics_1,sidearm_length_1,sidearm_species_1,\
                  sidearm_composition_1,sidearm_coverage_1,graft_density]


sidearm_statistics_2 = 'FJC'
#Nsc = 10
sidearm_length_2 = 9 # This is a quirk from continuous chains where the N_bb = L+1 
sidearm_species_2 = [2]
sidearm_composition_2 = [1.0]
sidearm_coverage_2 = 0.5
graft_density = 1  #currently this doesn't do anything...currently only does graft density of 1
sidechain_2_list = [sidearm_statistics_2,sidearm_length_2,sidearm_species_2,\
                  sidearm_composition_2,sidearm_coverage_2,graft_density]


chain_list = [backbone_list,sidechain_1_list,sidechain_2_list,sidechain_2_list]


b = [1.,1.,1.]  
#chiN12,chiN13,chiN14,chiN23,chiN24,chiN34
chiN = [0.1,0.1,0] 

chain_label = 'AB-Bottlebrush'
field = 'fields.in'

d = 1 #dimension
dt = 0.001 #overall "time" step
dS = 0.01 # if using CGC if using discrete chains should include it still, but will be not be used
stress_scale = 0.001
force_scale = [1.0,0.5]
cellscale = 10 #scaling the cell
npw = [128] #Resolution
Nref = 1 #Nref value
initial_box = [10/np.sqrt(Nref)] #Intial Box size
#Some Adjustable Values
INFILE='LAM.in' #Input File Name
diffuser_method='SOS' #Again not used for discrete chains, but still should be included but will not be used
invzeta=0.1/Nref #compressibility 
force_tol=1e-6 #Force tolerance
stress_tol=1e-5 #stress tolerance
ends = True #Ends without grafts?
non_primitive_centering = 'True' #primitive cell?
symmetrize = 'on' #symmetrizer 

 
space_group =  None #None for LAM
add_phase = True
if space_group is None:
    add_phase = False
    initial_box = [cellscale*initial_box[0]]



#mostly "default" stuff
inputfileversion = 3 #default in polyfts
nummodel = 1 #default in polyfts
supported_statistics = ['CGC','DGC','FJC'] #supported statistics
partition_function = 'canonical'
ham_bool = 'True'
stress_bool = 'True'
chem_pot_bool = 'False'
idealgas_bool = 'False'
field_type = 'HFields'
field_updater = 'EMPEC'
timesteps_block = 1000
block_num = 1000
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


#Depending on chain statistics do we have to add anything to the input script
statistics_models_add = [[f'contourds = {dS : 0.5f}',\
                     f'diffusermethod = {diffuser_method}'],\
                    [f'polymerReferenceN = {Nref}.'],\
                    [f'polymerReferenceN = {Nref}.']]

statistics_interactions_add = [[],\
                               [f'compressibility_invzetaN = {invzeta}'],
                               [f'compressibility_invzetaN = {invzeta}']]

#count species
components,statistics_list,n_species = Component_Statistics_Counter(chain_list)
modeltype = 'Polymer'



#Add statistics extra inputs 
models_add,interation_add,n_sidearm_types =\
    Add_Model_Statistics(statistics_list,chain_list,\
                             supported_statistics,statistics_models_add,\
                             statistics_interactions_add)

#Currently Grafting density only works for 1.0, but there are ways to bypass this section
#Makes Chains section in inputfile
def Grafting_Determinator(chain_list,supported_statistics,ends):
    j = 0
    if ends:
        chaingraftinfo= np.zeros((3,len(chain_list)-1))
        
        for i in chain_list:
            graftingdensity = i[-1]
            if len(i)==4:
                bb_length = i[1]
            if len(i)!=4:
                chaingraftinfo[0,j]= i[4]*(bb_length+1)
                unitcheck = abs(np.around(chaingraftinfo[0,j])-chaingraftinfo[0,j])
                if unitcheck>1e-6:
                    print('Warning Grafting is not precise')
                chaingraftinfo[0,j] = np.around(chaingraftinfo[0,j])
                if j==0:
                    chaingraftinfo[2,j] =  chaingraftinfo[0,j]-1
                else:
                    chaingraftinfo[1,j] =  chaingraftinfo[2,j-1]+1         
                    chaingraftinfo[2,j] =  chaingraftinfo[0,j]+chaingraftinfo[1,j]-1
                i.append(chaingraftinfo[0,j])
                i.append(chaingraftinfo[1,j])
                i.append(chaingraftinfo[2,j])
                print(i)
                j+=1
    if ends==False:
        chaingraftinfo= np.zeros((3,len(chain_list)-1))
        
        for i in chain_list:
            if len(i)==4:
                bb_length = i[1]
            if len(i)!=4:
                chaingraftinfo[0,j]= i[4]*(bb_length-1)
                if chaingraftinfo[0,j].is_integer()!=True:
                    print('Warning Grafting is not precise')
                chaingraftinfo[0,j] = np.around(chaingraftinfo[0,j])
                if j==0:
                    chaingraftinfo[1,j] =  1
                    chaingraftinfo[2,j] =  chaingraftinfo[0,j]
                else:
                    chaingraftinfo[1,j] =  chaingraftinfo[2,j-1]+1         
                    chaingraftinfo[2,j] =  chaingraftinfo[0,j]+chaingraftinfo[1,j]-1
                i.append(chaingraftinfo[0,j])
                i.append(chaingraftinfo[1,j])
                i.append(chaingraftinfo[2,j])
                j+=1

    j = 0
    chain_text_list = []
    for i in chain_list:
        chain_text_list.\
        append(Chain_Builder(i,j,supported_statistics).Chain_Text_Write())
        j+=1
    return chain_list,chain_text_list

    
    
    
    
    
    
    
    
chain_list,chain_text_list = Grafting_Determinator(chain_list, supported_statistics, ends)

#Turn everything into text to be put into Input Builder
chiN_list =chiN_generator(chiN,components)
kuhn_length_text = parameter_species(b,n_species)
force_scale_text = parameter_species(force_scale,n_species)
initial_box_text = parameter_species(initial_box,d)
npw_text = parameter_species(npw,d)


#Open and Write Input File with everything inside
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


#read datafile


op = open(input_file_path,'r')
data_file = op.read()
op.close()