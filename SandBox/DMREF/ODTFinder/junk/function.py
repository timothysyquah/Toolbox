#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 11:21:50 2020

@author: tquah
"""


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
        continue
    else:
        os.makedirs(WDIR)
    return WDIR