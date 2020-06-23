#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 14 12:17:12 2020

@author: tquah
"""
    
import numpy as np
from itertools import combinations 



class Chain_Builder():
    def __init__(self,model_list,num,supported_statistics):
        statistics_chain_add =\
        ['blockfractions','nperblock','nperblock']
        statistics_chain_length =\
        ['length','nbeads','nbeads']
        statistics_length_addition =\
        [0,1,1]
        
        self.statistics = model_list[0]
        self.length = model_list[1]
        self.nblocks = len(model_list[2])
        self.species = list_to_text(model_list[2])
        self.supported_statistics = supported_statistics
        self.chain_template ='''
{self.chainname} {{ 
    statistics   = {self.statistics}
    {self.chain_length}       =  {self.length}
    nblocks      = {self.nblocks}
    blockspecies = {self.species}
    {self.composition_name}  = {self.composition}  
'''
        self.sidearm_addition = '''
        numarms      = {self.n_arms}
        backbonegraftingstart = {self.graft_start}
        backbonegraftingend   = {self.graft_end}
        }}
'''

        j = 0
        for i in self.supported_statistics:
            if self.statistics ==i:
                self.chain_length =statistics_chain_length[j]
                self.composition_name = statistics_chain_add[j]
                if len(model_list)==4:
                    self.length =  model_list[1]+statistics_length_addition[j]
                    
                    self.chaintype = 'backbone'
                    self.chainname = self.chaintype
#                chain_length = 
                elif len(model_list)>1:
                    self.length =  model_list[1]
                    self.chaintype = 'sidearmtype'
                    self.chainname = self.chaintype+str(num)
                    self.n_arms = int(model_list[6])
                    self.graft_start= int(model_list[7])
                    self.graft_end= int(model_list[8])
                statistics_length_mult = [1,self.length,self.length]
                composition =\
                np.array(model_list[3])*statistics_length_mult[j]
                if self.statistics=='DGC' or 'FJC':
                    composition = np.around(composition)
                self.composition =\
                list_to_text(composition)
            j+=1

#                statistics_length_mult[i]*
    def Chain_Text_Write(self):
       
        if self.chaintype=='backbone':
            return eval(f"f'''{self.chain_template}'''")+\
                  '              }'
        if self.chaintype=='sidearmtype':
            return eval(f"f'''{self.chain_template}'''")+\
            eval(f"f'''{self.sidearm_addition}'''")

class Input_Builder():
    def __init__(self,f,inputfileversion,nummodel,modeltype,n_species,\
                     kuhn_length_text,chain_label,n_sidearm_types,\
                     models_add,chain_text_list,chiN_list,interation_add,\
                     d,initial_box,npw,partition_function,\
                     ham_bool,stress_bool,chem_pot_bool,idealgas_bool,\
                     field_type,field,field_updater,cell_updater,timesteps_block,\
                     block_num,dt,force_scale_text,stress_scale,force_tol,stress_tol,\
                     variablecell_bool,density_history_bool,field_history_bool,\
                     density_chain_bool,format_fields_bool,volfrac_chain,cellscale,
                     add_phase,space_group,non_primitive_centering,\
                     symmetrize,parallel_cuda,cuda_thread_block_size,nThreads):
        
        self.f = f
        self.inputfileversion = inputfileversion
        self.nummodel = nummodel
        self.modeltype = modeltype
        self.n_species = n_species
        self.kuhn_length_text = kuhn_length_text
        self.chain_label = chain_label
        self.n_sidearm_types = n_sidearm_types
        self.models_add = models_add
        self.chain_text_list = chain_text_list
        self.chiN_list = chiN_list
        self.interation_add = interation_add
        self.d = d
        self.initial_box =initial_box
        self.npw = npw
        self.partition_function = partition_function
        self.ham_bool = ham_bool
        self.stress_bool = stress_bool
        self.chem_pot_bool = chem_pot_bool
        self.idealgas_bool = idealgas_bool
        self.field_type = field_type
        self.field = field
        self.field_updater = field_updater
        self.cell_updater = cell_updater
        self.timesteps_block = timesteps_block
        self.block_num = block_num
        self.dt = dt
        self.force_scale_text = force_scale_text
        self.stress_scale = stress_scale
        self.force_tol = force_tol
        self.stress_tol = stress_tol
        self.variablecell_bool = variablecell_bool
        self.density_history_bool = density_history_bool
        self.field_history_bool = field_history_bool
        self.density_chain_bool = density_chain_bool
        self.format_fields_bool = format_fields_bool
        self.volfrac_chain = volfrac_chain
        self.add_phase = add_phase
        self.space_group = space_group
        self.cellscale = cellscale
        self.non_primitive_centering = non_primitive_centering
        self.symmetrize = symmetrize
        self.parallel_cuda = parallel_cuda
        self.cuda_thread_block_size = cuda_thread_block_size
        self.nThreads = nThreads
    def Write_to_Text(self):
            self.f.write(f'''InputFileVersion={self.inputfileversion}
            
models {{
    NumModels = {self.nummodel}
    ModelType = {self.modeltype}
    
    monomers {{
        nspecies = {self.n_species}
        kuhnlen  = {self.kuhn_length_text}
        }}
    
    chains {{
        nchains   = 1 
''')
            
            #chains remember to close with }
            for i in self.models_add:
                self.f.write('        '+i)
                self.f.write('\n')
            
            self.f.write(f'''        
  chain1 {{
    label  = {self.chain_label}            
    architecture = comb
    NumSideArmTypes = {self.n_sidearm_types}
''')
            for i in self.chain_text_list:
                self.f.write(i)
            self.f.write('\n')
            self.f.write('       }')
            self.f.write('\n')
            self.f.write('}')

            self.f.write(f'''
    model1 {{
        cell {{
        dim = {self.d}
        celllengths =  {self.initial_box}
        npw = {self.npw}
''')

            if self.add_phase:
                self.f.write(f'''
        cellscaling = {self.cellscale} 
        spacegroupname = {self.space_group}
        #ApplyNonPrimitiveCentering = {self.non_primitive_centering}
        CenterToPrimitiveCell = {self.non_primitive_centering}
        symmetrize = {self.symmetrize}

                         ''')
            self.f.write(f'''
        }}

        interactions {{
                         
        ''')
            for i in self.chiN_list:
                self.f.write('        ')
                self.f.write(i)
                self.f.write('\n')
            
            for i in self.interation_add:
                self.f.write('        ')
                self.f.write(i)
                self.f.write('\n')
            self.f.write('         }')
            self.f.write('\n')
            self.f.write(f'''
    composition {{
        ensemble     =  {self.partition_function}
        chainvolfrac = {self.volfrac_chain}
    }}

    operators {{
      CalcHamiltonian       = {self.ham_bool}
      CalcStressTensor      = {self.stress_bool}
      CalcChemicalPotential = {self.chem_pot_bool}
      IncludeIdealGasTerms  = {self.idealgas_bool}
    }}

    initfields {{
      ReadInputFields = {self.field_type}
      InputFieldsFile = {self.field}
    }}
  }}
}}
    
simulation {{
  jobtype = SCFT

  FieldUpdater = {self.field_updater}
  CellUpdater = {self.cell_updater}
  NumTimeStepsPerBlock = {self.timesteps_block}
  NumBlocks = {self.block_num}

  TimeStepDT = {self.dt : 0.5e}
  lambdaForceScale = {self.force_scale_text}
  lambdaStressScale = {self.stress_scale: 0.5f}
  SCFTForceStoppingTol = {self.force_tol: 0.5e}
  SCFTStressStoppingTol = {self.stress_tol : 0.5e}

  VariableCell = {self.variablecell_bool}

  IO {{
    KeepDensityHistory   = {self.density_history_bool}
    KeepFieldHistory     = {self.field_history_bool}
    DensityOutputByChain = {self.density_chain_bool}
    OutputFormattedFields = {self.format_fields_bool}

    OutputFields         = {self.field_type}
    FieldOutputSpace     = both  # rspace, kspace or both
  }}
}}

parallel {{
  CUDA_selectdevice = {self.parallel_cuda}
  CUDA_threadblocksize = {self.cuda_thread_block_size}

  OpenMP_nthreads = {self.nThreads}
}}
''')
    



def Input_Standard(input_file_path,n_species,kuhn_length_text,chain_label,n_sidearm_types,\
                   models_add,chain_text_list,chiN_list,interation_add,d,\
                   initial_box,npw,field,dt,force_scale_text,stress_scale,\
                   force_tol,stress_tol,cell_updater,cellscale,add_phase,\
                 space_group,non_primitive_centering,symmetrize,\
                 parallel_cuda,cuda_thread_block_size,nThreads):
    
    if n_species==2:
        modeltype = 'BlockPolymerMelt2Spec'
    if n_species>=3:
        modeltype = 'BlockPolymerMelt'
    inputfileversion = 3
    nummodel = 1
    partition_function = 'canonical'
    ham_bool = 'True'
    stress_bool = 'True'
    chem_pot_bool = 'False'
    idealgas_bool = 'False'
    field_type = 'HFields'
    field_updater = 'SIS'
    timesteps_block = 1000
    block_num = 1000
    variablecell_bool = 'True'
    density_history_bool = 'False'
    field_history_bool = 'False'
    density_chain_bool = 'False'
    format_fields_bool = 'True'
    volfrac_chain = 1.0
    
    f = open(input_file_path,'w+')    
    Input_Builder(f,inputfileversion,nummodel,modeltype,n_species,\
                     kuhn_length_text,chain_label,n_sidearm_types,\
                     models_add,chain_text_list,chiN_list,interation_add,\
                     d,initial_box,npw,partition_function,\
                     ham_bool,stress_bool,chem_pot_bool,idealgas_bool,\
                     field_type,field,field_updater,cell_updater,timesteps_block,\
                     block_num,dt,force_scale_text,stress_scale,force_tol,stress_tol,\
                     variablecell_bool,density_history_bool,field_history_bool,\
                     density_chain_bool,format_fields_bool,volfrac_chain,cellscale,
                     add_phase,space_group,non_primitive_centering,\
                     symmetrize,parallel_cuda,cuda_thread_block_size,nThreads).Write_to_Text()
    f.close()

def list_to_text(lst):
    text = ''
    for i in lst:
        text+=f'{i} '
    return text

def Component_Statistics_Counter(chain_list):
    statistics_list = []
    components = []
    for i in chain_list:
        statistics_list.append(i[0])
        components+=(i[2][:])
    components = list(set(components))
    species = len(components)
    return components,statistics_list,species

def Add_Model_Statistics(statistics_list,chain_list,\
                         supported_statistics,statistics_models_add,\
                         statistics_interactions_add):
    j = 0
    models_add = []
    interation_add = []
    for i in supported_statistics:
        
        if len([ele for ele in statistics_list if(ele in i)])>0:
            models_add+=statistics_models_add[j]    
            interation_add+=statistics_interactions_add[j]
        j+=1
    
    models_add = list(set(models_add))
    interation_add = list(set(interation_add))
    n_sidearm_types = len(chain_list)-1
    return models_add,interation_add,n_sidearm_types

def Grafting_Determinator(chain_list,supported_statistics,ends):
    j = 0
    if ends:
        chaingraftinfo= np.zeros((3,len(chain_list)-1))
        
        for i in chain_list:
            spacing = i[-1]
            if len(i)==4:
                bb_length = i[1]
            if len(i)!=4:
                chaingraftinfo[0,j]= i[4]*(bb_length+1)
                unitcheck = abs(np.around(chaingraftinfo[0,j])-chaingraftinfo[0,j])
                if unitcheck>1e-6:
                    print('Warning Grafting is not precise')
                chaingraftinfo[0,j] = np.floor(chaingraftinfo[0,j])
                if j==0:
                    chaingraftinfo[2,j] =  chaingraftinfo[0,j]-1
                else:
                    chaingraftinfo[1,j] =  chaingraftinfo[2,j-1]+1         
                    chaingraftinfo[2,j] =  chaingraftinfo[0,j]+chaingraftinfo[1,j]-1
                i.append(chaingraftinfo[0,j])
                i.append(chaingraftinfo[1,j])
                i.append(chaingraftinfo[2,j])
                j+=1
        # print(chaingraftinfo)

                
    if ends==False:
        chaingraftinfo= np.zeros((3,len(chain_list)-1))
        
        for i in chain_list:
            if len(i)==4:
                bb_length = i[1]
            if len(i)!=4:
                chaingraftinfo[0,j]= i[4]*(bb_length-1)
                if chaingraftinfo[0,j].is_integer()!=True:
                    print('Warning Grafting is not precise')
                chaingraftinfo[0,j] = np.floor(chaingraftinfo[0,j])
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

def chiN_generator(chiN,components):
    chi_combo = list(combinations(components,2))
    chiN_list = []
    for i in range(0,len(chi_combo)):
        comb = str(chi_combo[i][0])+str(chi_combo[i][1])
        chivalues = chiN[i]
        
        chiN_list.append(f'''chiN{comb} = {chivalues}''')
    return chiN_list

def parameter_species(paramlist,n_species):
    if len(paramlist)!=n_species:
        while len(paramlist)<n_species:
            paramlist.append(paramlist[-1])
    param_text = list_to_text(paramlist)
    return param_text


def Lazy_Input_Generator(input_file_path,field,chain_list,chiN,dS,npw,dt,\
                         initial_box,stress_scale,force_scale,d,chain_label,\
                         ends,diffuser_method,Nref,invzeta,\
                         kuhn_length,stress_tol,force_tol,CellUpdater,cellscale,
                         add_phase,space_group,non_primitive_centering,\
                         symmetrize,parallel_cuda,cuda_thread_block_size,nThreads):
    
    statistics_models_add = [[f'contourds = {dS : 0.5f}',\
                             f'diffusermethod = {diffuser_method}'],\
                            [f'polymerReferenceN = {Nref}.'],\
                            [f'polymerReferenceN = {Nref}.']]
    
    statistics_interactions_add = [[],\
                                   [f'compressibility_invzetaN = {invzeta}'],
                                   [f'compressibility_invzetaN = {invzeta}']]
    supported_statistics = ['CGC','DGC','FJC']
    
    
    components,statistics_list,n_species = Component_Statistics_Counter(chain_list)
    
    
    models_add,interation_add,n_sidearm_types =\
        Add_Model_Statistics(statistics_list,chain_list,\
                                 supported_statistics,statistics_models_add,\
                                 statistics_interactions_add)
    chain_list,chain_text_list = Grafting_Determinator(chain_list,supported_statistics,ends)
    chiN_list =chiN_generator(chiN,components)
    kuhn_length_text = parameter_species(kuhn_length,n_species)
    force_scale_text = parameter_species(force_scale,n_species)
    
    Input_Standard(input_file_path,n_species,kuhn_length_text,chain_label,n_sidearm_types,\
                   models_add,chain_text_list,chiN_list,interation_add,d,\
                   initial_box,npw,field,dt,force_scale_text,stress_scale,\
                   force_tol,stress_tol,CellUpdater,add_phase,\
                 space_group,non_primitive_centering,symmetrize,parallel_cuda,\
                     cuda_thread_block_size,nThreads)

def Lazy_Input_Generator_General(input_file_path,field,chain_list,chiN,dS,npw,dt,\
                                 initial_box,stress_scale,force_scale,d,chain_label,\
                                 ends,diffuser_method,Nref,invzeta,\
                                 kuhn_length,stress_tol,force_tol,CellUpdater,cellscale,
                                 add_phase,space_group,non_primitive_centering,\
                                 symmetrize,parallel_cuda,cuda_thread_block_size,nThreads):
    
    statistics_models_add = [[f'contourds = {dS : 0.5f}',\
                             f'diffusermethod = {diffuser_method}'],\
                            [f'polymerReferenceN = {Nref}.'],\
                            [f'polymerReferenceN = {Nref}.']]
    
    statistics_interactions_add = [[],\
                                   [f'compressibility_invzetaN = {invzeta}'],
                                   [f'compressibility_invzetaN = {invzeta}']]
    supported_statistics = ['CGC','DGC','FJC']
    
    
    components,statistics_list,n_species = Component_Statistics_Counter(chain_list)
    
    
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

    Input_Standard(input_file_path,n_species,kuhn_length_text,chain_label,n_sidearm_types,\
                   models_add,chain_text_list,chiN_list,interation_add,d,\
                   initial_box_text,npw_text,field,dt,force_scale_text,stress_scale,\
                   force_tol,stress_tol,CellUpdater,cellscale,add_phase,\
                 space_group,non_primitive_centering,symmetrize,\
                 parallel_cuda,cuda_thread_block_size,nThreads)