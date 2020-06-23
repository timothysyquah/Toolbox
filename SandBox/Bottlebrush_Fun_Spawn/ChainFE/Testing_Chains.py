#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 14:35:43 2020

@author: tquah
"""
import os
IDIR = os.getcwd()
components = IDIR.split('/')
path_to_tools = '/'+components[1]+'/'+components[2]+'/.timtools'
op = open(path_to_tools,'r')
chainpath =op.read()+'/geninput/'
op.close()


import sys
sys.path.append(chainpath)

from Chains import Chain_Builder

import numpy as np





backbone_composition = [1.0]
backbone_species = [1]
backbone_length = 99
backbone_statistics = 'FJC'

backbone_list = [backbone_statistics,backbone_length,backbone_species,\
                 backbone_composition]

sidearm_composition = [1.0]
sidearm_species = [2]
sidearm_length = 5.0
sidearm_statistics = backbone_statistics
sidearm_coverage1 = 1/4
spacing = 1
sidechain_1_list = [sidearm_statistics,sidearm_length,sidearm_species,\
                 sidearm_composition,sidearm_coverage1,spacing]

sidearm_composition = [1.0]
sidearm_species = [3]
sidearm_length = sidearm_length
sidearm_statistics = backbone_statistics
sidearm_coverage = 1/2
sidechain_2_list = [sidearm_statistics,sidearm_length,sidearm_species,\
                 sidearm_composition,sidearm_coverage,spacing]

sidearm_composition = [1.0]
sidearm_species = [4]
sidearm_length = sidearm_length
sidearm_statistics = backbone_statistics
sidearm_coverage = 1/4
sidechain_3_list = [sidearm_statistics,sidearm_length,sidearm_species,\
                 sidearm_composition,sidearm_coverage,spacing]

chain_list = [backbone_list,sidechain_1_list,sidechain_2_list,sidechain_3_list]
sstats = ['CGC','DGC','FJC']
ends = True






j = 0
if ends:
    print('a')
    chaingraftinfo= np.zeros((3,len(chain_list)-1))
    
    for i in chain_list:
        if len(i)==4:
            bb_length = i[1]
        if len(i)!=4:
            chaingraftinfo[0,j]= i[4]*(bb_length+1)
            if chaingraftinfo[0,j].is_integer()!=True:
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
            
if ends==False:
    print('b')
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
    append(Chain_Builder(i,j,sstats).Chain_Text_Write())
    j+=1
