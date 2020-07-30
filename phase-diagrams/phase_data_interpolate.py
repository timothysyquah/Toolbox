#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 23:38:23 2020

@author: tquah
"""
import pickle
import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import pandas as pd
import argparse
import glob
import shutil    


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool to compute phase boundaries')
    parser.add_argument('-dp', '--dict_path', action='store', default='data.dict',help='Dictionary import that has all Free energy data')
    parser.add_argument('-dir', '--dirname', action='store', default='Interpolated',help='Directory export name')
    parser.add_argument('-e', '--export_path', action='store', default=os.getcwd(),help='Export Path for mass files')
    parser.add_argument('-n', '--Neff', action='store', default=2100,help='N_eff values',type = float)
    parser.add_argument('-dx', '--dx', action='store', default=0.005,help='dx to create',type = float)
    parser.add_argument('-r', '--reference_phase', action='store', default='',help='Reference phase',type = str)
    parser.add_argument('-diagnostic', '--plot_diagnositic', action='store', default=False,help='Diagnostic Tools for Phase Diagrams',type = bool)
    parser.add_argument('-k', '--keywrd', action='store', default=[''],help='Keyword used to name dimensions',type = str)
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')

    args = parser.parse_args()
    
    args.dict_path = '/home/tquah/IMPORT_BRAID/diblock_phasediagram/data.dict'
    args.export_path = '/home/tquah/IMPORT_BRAID/'
    args.keywrd = ['chi','f']

    # os.chdir(args.dict_path)
    export_path_full = os.path.join(args.export_path,args.dirname)
    with open(args.dict_path, 'rb') as handle:
        F0dat = pickle.load(handle)

    F0_interp = dict()

    extract_keys = []
    for i in range(0,len(list(F0dat)[0])):
        extract_keys.append([])
    
    
    
    for i in range(0,len(F0dat)):
        for j in range(0,len(list(F0dat)[0])):
            extract_keys[j].append(list(F0dat)[i][j])
    for i in range(0,len(list(F0dat)[0])):
        extract_keys[i] = sorted(list(set(extract_keys[i])))
        
    
    oldindex = extract_keys[-1].index('DIS')
    
    extract_keys[-1].insert(0, extract_keys[-1].pop(oldindex))

    if len(list(F0dat)[0])==2:
        
        for i in range(0,len(extract_keys[1])):
            for j in range(0,len(extract_keys[0])):
                if (extract_keys[0][j],extract_keys[1][i]) in F0dat:
                    data_array = F0dat[extract_keys[0][j],extract_keys[1][i]]
                    if len(np.shape(data_array))>1:
                        data_array = data_array[data_array[:,0].argsort()]
                        
                        fmin = np.min(data_array[:,0])
                        fmax = np.max(data_array[:,0])
                        fun = interp1d(data_array[:,0],data_array[:,1])
                        farray = np.around(np.arange(fmin,fmax+1e-6,args.dx),5)
                        yarray = fun(farray)
                        
                        for k in range(0,len(farray)):
                            chiN = extract_keys[0][j]*args.Neff
                            
                            fullexportname = f'{args.keywrd[0]}_{chiN:0.3f}/{args.keywrd[1]}{farray[k]:0.5f}' 
                            fullexportpath = os.path.join(export_path_full,fullexportname)
                            if os.path.exists(fullexportpath):
                                pass
                            else:
                                os.makedirs(fullexportpath)
                            filepath = os.path.join(fullexportpath,args.filename)
                            if extract_keys[1][i]=='DIS':                                                                    
                                f0 = open(filepath,'w')
                                f0.write(f'{extract_keys[1][i]}Phase {yarray[k]} 2')
                                f0.write('\n')
                                f0.close()
                            else:
                                f0 = open(filepath,'a')
                                f0.write(f'{extract_keys[1][i]}Phase {yarray[k]} 2')
                                f0.write('\n')
                                f0.close()
    else:
        print('Dimension is not 2')
    





    # for i in range(0,len(list(F0dat))):
    #     F0_interp[list(F0dat)] = dict()
    #     if len(np.shape(F0dat[list(F0dat)[i]]))>1:
    #         array = F0dat[list(F0dat)[i]][F0dat[list(F0dat)[i]][:,0].argsort()]
            
            
            
    #     else:
    #         array = F0dat[list(F0dat)[i]]
