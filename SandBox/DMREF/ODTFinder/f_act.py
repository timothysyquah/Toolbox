#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 12 14:58:56 2020

@author: tquah
"""

import numpy as np
from itertools import combinations

def check_equality(a,b,message):
    assert a==b, message

def listcheck(l,message):
    uniquelist = [",".join(map(str, comb)) for comb in combinations(l, 2)]
    for pair in uniquelist:
        pair_elements = pair.split(',')
        check_equality(pair_elements[0],pair_elements[1],message)


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
    


phase = ['LAM','DIS']
parameter_dict = dict()

parameter_dict['chi'] = sweep_list_return([0.1],[0.2],[0.05],5)
parameter_dict['f'] = sweep_list_return([0.9,0.1],[0.1,0.9],[-0.1,0.1],5)
parameter_dict['nsc'] = sweep_list_return([20,20],[10,30],[-2,2],1)
parameter_dict['nbb'] = np.arange(49,149+1e-6,50)
parameter_dict['nref'] = np.array([1])

itterative_structure = ['chi','nbb','nsc','nref','f']
countlist = []
for i in range(0,len(itterative_structure),1):
    countlist.append(sweep_counter(parameter_dict,itterative_structure[i]))
for q in range(0,len(phase)):

    for m in range(0,countlist[0]):
        
        for n in range(0,countlist[1]):
        
            for o in range(0,countlist[2]):
            
                for p in range(0,countlist[3]):
                    
                    for r in range(0,countlist[4]):

                        Phase = phase[q]
                        itterlist = [m,n,o,p,r]
                        nbb_loc = itterative_structure.index('nbb')
                        chi_loc = itterative_structure.index('chi')
                        nref_loc = itterative_structure.index('nref')
                        nsc_loc = itterative_structure.index('nsc')
                        f_loc = itterative_structure.index('f')
                        
                        
                        f_act = (parameter_dict['nbb'][itterlist[nbb_loc]]+1)*parameter_dict['f'][0]
                        print(f_act)

