#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 21 10:26:32 2020

@author: tquah
"""
import os
IDIR = os.getcwd()
components = IDIR.split('/')
path_to_tools = '/'+components[1]+'/'+components[2]+'/.timtools'
op = open(path_to_tools,'w+')
op.write(IDIR)
op.close()

