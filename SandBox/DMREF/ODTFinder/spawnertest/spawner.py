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
dNbb = 50
Nscmin = 0
Nscmax = 20
dNsc = 5
chiNmin = 23
chiNmax = 10
dchiN = 0.5
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
        
        Neff = Nbb*(Nsc+1)
        chi_min = np.round(chiNmin/Neff,8)
        chi_max = np.round(chiNmax/Neff,8)  

        if chiNmin>chiNmax:
            dchi = -np.round(dchiN/Neff,8)
        else:
            dchi = np.round(dchiN/Neff,8)

        
        replace_text = ['__nbbmin__','__nbbmax__','__dnbb__','__nscmin__',\
                        '__nscmax__','__dnsc__','__chimin__','__chimax__','__dchi__']
        
        replace_values = [Nbb-1,Nbb-1,1,Nsc,Nsc,1,chi_min,chi_max,dchi]
        
        rewrite_submit(outfile,infile,replace_text,replace_values)
        
        print('running')
        import subprocess
        cmd=f"qsub {outfile}"
        # cmd="sbatch submit.sh"
        subprocess.call(cmd.split())
        

        
    