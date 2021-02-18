#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 20:32:55 2020

@author: tquah
"""


import numpy as np
import os
#nbb scan


chain_stats_list = 'DGC'
dt = 1e-3
npw = 512
gauss_width = 0.2
inv_zeta = 0.001
Nsc_list = np.array([0,5,10,20,40])
chi = 0.1
Nbb = 100
submit = False

def rewrite_submit(outfile,infile,replace_text,replace_values):
    with open('%s/%s' % (IDIR,outfile), 'w') as fout:
        with open('%s/%s' % (IDIR,infile),'r') as f:
            for line in f:
                for i in range(0,len(replace_text)):
                    line = line.replace(replace_text[i],str(replace_values[i]))
                fout.write(line)
    fout.close
    f.close()

def Ntot_calc(Nbb,Nsc):
    return Nbb*(Nsc+1)

IDIR = os.getcwd()


infile = 'submitnew.sh'

for Nsc in Nsc_list:
    if Ntot_calc(Nbb,Nsc)*chi<10.5:
        Nbb = 13/(chi*(Nsc+1))
    
    outfile = f'submit_{Nsc}.sh'

    replace_text = ['__chain_stats__','__dt__','__npw__','__nbbstart__','__gwidth__','__invzeta__',\
                    '__nsc__','__chi__']
    replace_values = [chain_stats_list,f'{dt:0.3e}',f'{npw}',f'{Nbb-1}',f'{gauss_width}',f'{inv_zeta}',f'{Nsc}',f'{chi:0.3f}']
    rewrite_submit(outfile,infile,replace_text,replace_values)
    if submit:
        import subprocess
        cmd=f"qsub {outfile}"
        subprocess.call(cmd.split())
