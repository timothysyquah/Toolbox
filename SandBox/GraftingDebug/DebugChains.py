#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 25 23:13:01 2020

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

# from Chains import Input_Builder,Component_Statistics_Counter,\
#     Add_Model_Statistics,Grafting_Determinator,chiN_generator,parameter_species
from Chains import Chain_Builder
    
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
                unitcheck = abs(int(chaingraftinfo[0,j])-chaingraftinfo[0,j])
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
                print(i)    
                j+=1
    
    j = 0
    chain_text_list = []
    for i in chain_list:
        print(i)
        chain_text_list.\
        append(Chain_Builder(i,j,supported_statistics).Chain_Text_Write())
        j+=1
    return chain_list,chain_text_list

    
    
fAfCs=0.21
sidearm_coverage_1 = np.round(fAfCs,5)
sidearm_coverage_2 = np.round(1-(2*fAfCs),5)
sidearm_coverage_3 = np.round(fAfCs,5)
backbone_statistics = 'DGC'
backbone_length=99
backbone_species = [1,2,3]

backbone_composition = [sidearm_coverage_1,sidearm_coverage_2,sidearm_coverage_3]


backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                  backbone_composition]
chain_list = [backbone_list]









ends = True
supported_statistics = ['CGC','DGC','FJC']
chain_list,chain_text_list = Grafting_Determinator(chain_list,supported_statistics,ends)
