#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Jul 25 20:32:55 2020

@author: tquah
"""


import numpy as np
import os
import shutil
import time
#nbb scan
Nbbmin = 50
Nbbmax = 200
dNbb = 10
Nscmin = 0
Nscmax = 40
dNsc = 5
dchiN = 0.01
dt = 0.001

param1 = np.array([26.93912298,  0.85598599, 11.26917458])
param2 = np.array([22.04110061,  0.85598598,  9.22023371])
def curvefitfun(x,a,b,c):
    return a*x**b+c

def sweep_generator(a_min,a_max,a_d):
    return np.arange(a_min,a_max+1e-6,a_d,dtype=int)

def rewrite_submit(outfile,infile,replace_text,replace_values):
    with open('%s/%s' % (IDIR,outfile), 'w') as fout:
        with open('%s/%s' % (IDIR,infile),'r') as f:
            for line in f:
                for i in range(0,len(replace_text)):
                    line = line.replace(replace_text[i],str(replace_values[i]))
                fout.write(line)
    fout.close
    f.close()

N_bb_array = sweep_generator(Nbbmin,Nbbmax,dNbb)
N_sc_array = sweep_generator(Nscmin,Nscmax,dNsc)
IDIR = os.getcwd()
infile = 'submitnew_editable.sh'
outfile = 'submit.sh'
for Nbb in N_bb_array:
    for Nsc in N_sc_array:
        
        
        alpha = Nsc/Nbb
        chiNmin = curvefitfun(alpha,param1[0],param1[1],param1[2])
        chiNmax = curvefitfun(alpha,param2[0],param2[1],param2[2])

        Neff = Nbb*(Nsc+1)
        chi_min = np.round(chiNmin/Neff,8)
        chi_max = np.round(chiNmax/Neff,8)  

        if chiNmin>chiNmax:
            dchi = -float(np.round(dchiN/Neff,8))
        else:
            dchi = float(np.round(dchiN/Neff,8))

        
        replace_text = ['__dt__','__nbbmin__','__nbbmax__','__dnbb__','__nscmin__',\
                        '__nscmax__','__dnsc__','__chimin__','__chimax__','__dchi__']
        
        replace_values = [f'{dt:0.8f}',f'{Nbb-1}',f'{Nbb-1}',1,f'{Nsc}',f'{Nsc}',1,f'{chi_min:0.8f}',f'{chi_max:0.8f}',f'{dchi:0.8f}']
        
        rewrite_submit(outfile,infile,replace_text,replace_values)
        
        # print('running')
        # import subprocess
        # cmd=f"qsub {outfile}"
        # cmd="sbatch submit.sh"
        # subprocess.call(cmd.split())
        # import subprocess
        # cmd=f"./test1.sh"
        # subprocess.call(cmd.split())


        
    