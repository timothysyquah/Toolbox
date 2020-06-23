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
chainpath =op.read()+'/geninput/'
op.close()


import sys
sys.path.append(chainpath)

from Chains import Lazy_Input_Generator


def make_input(PHASE,REPLACE_lIST,WDIR):
    diffuser_method='SOS'
    CellUpdater = 'Broyden'
    Nref=1
    invzeta=0.1
    force_tol=1e-5
    stress_tol=1e-4
    kuhn_length=[1.0]
    field = REPLACE_lIST[0]
    chain_list = REPLACE_lIST[1] 
    chiN = REPLACE_lIST[2] 
    npw = REPLACE_lIST[3] 
    initial_box = REPLACE_lIST[4]
    d = REPLACE_lIST[5]
    chain_label = REPLACE_lIST[6]
    dS = REPLACE_lIST[7] 
    dt = REPLACE_lIST[8]
    stress_scale = REPLACE_lIST[9]
    force_scale = REPLACE_lIST[10]
    INFILE=f'{PHASE}.in'
    input_file_path = os.path.join(WDIR,INFILE)
    ends = False
    if len(REPLACE_lIST)>11:
        force_tol = REPLACE_lIST[11][0]
        stress_tol = REPLACE_lIST[11][1]
    if len(REPLACE_lIST)>12:
        CellUpdater = REPLACE_lIST[12]
    if len(REPLACE_lIST)>13:
        invzeta = REPLACE_lIST[13]
    add_phase = False
    space_group = None
    non_primitive_centering = None
    symmetrize = None
    parallel_cuda = 0
    cuda_thread_block_size = 128                     
    nThreads = 1
    Lazy_Input_Generator(input_file_path,field,chain_list,chiN,dS,npw,dt,\
                             initial_box,stress_scale,force_scale,d,chain_label,\
                             ends,diffuser_method,Nref,invzeta,\
                             kuhn_length,stress_tol,force_tol,CellUpdater,add_phase,\
                             space_group,non_primitive_centering,symmetrize,\
                             parallel_cuda,cuda_thread_block_size,nThreads)