#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 27 11:42:38 2020

@author: tquah
"""

import os
IDIR = os.getcwd()
components = IDIR.split('/')
path_to_tools = '/'+components[1]+'/'+components[2]+'/.timtools'
op = open(path_to_tools,'r')
toolpath = op.read()
newsubpath = toolpath+'/newsubmit/'

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

test = spacegroup_finder('BCC')