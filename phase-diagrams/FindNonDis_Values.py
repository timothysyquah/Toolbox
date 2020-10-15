#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct 14 10:46:33 2020

@author: tquah
"""

import numpy as np
import os 
import shutil
import glob
import argparse
# # os.chdir('/home/tquah/Projects/sweep-asym-armlength_corrected/')
# # print(os.listdir())
# os.chdir("/home/tquah/Projects/DMREF/sweep-asym-armlength_corrected/")
# # dirs = glob.glob("/home/tquah/Projects/DMREF/sweep-asym-armlength_corrected/chiAB_0.0289/Nsc*/fA*")
# # phaselist = ['FCC','BCC']
# os.chdir('/home/tquah/Projects/DMREF/sweep-asym-armlength_BCC_fix/')

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Tool to figure out what simulations to rerun')
    parser.add_argument('-e', '--exportname', action='store', default='resubmit.dat',help='exportdata')
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default=glob.glob("chi*/NscA_*/fA*"),help='list of directories that contain each phase point')
    parser.add_argument('-p', '--phaselist', action='store', default=[''], nargs='+', help='Phase List')
    parser.add_argument('-t', '--tol', action='store', default=1e-5, nargs='+', help='Tolerance',type=float)
    args = parser.parse_args()
    dirs = args.dirs
    exportname = args.exportname
    pwd = os.getcwd()
    seedlist = []
    exportpath = os.path.join(pwd,exportname)
    tol = args.tol
    phaselist = [a for a in args.phaselist ]
    oe = open(exportpath,'w+')
    for Phase in phaselist:
        print(Phase)
        ODT_directories = []
        ODT_SEED = []
        ratio = []
        chi = []
        for directory in dirs:
            info = directory.split('/')
            # fA.append(info[-1]) 
            ratio.append((float(info[-2].split('_')[1]),float(info[-2].split('_')[3])))
            chi.append(float(info[-3][-6:]))
        ratio= list(sorted(set(ratio)))
        chi= list(sorted(set(chi)))
        store_dict_low= dict()
        seeddict = dict()
        for i in chi:
            for j in ratio:
                store_dict_low[i,j] = []
                seeddict[i,j] = []
        for directory in dirs:
            newpath = os.path.join(pwd,directory)
            os.chdir(newpath)
            split_info = directory.split('/')
            chiab = float(split_info[-3][6:])
            ratioab = (float(split_info[-2].split('_')[1]),float(split_info[-2].split('_')[3]))
            fab = float(split_info[-1][2:])
            if args.filename in os.listdir():
                op = open(args.filename)
                rows = op.read().split('\n')
                op.close()
                
                freeenergy = []
                phase = []
                for row in rows:
                    if len(row)!=0:
                        freeenergy.append(float(row.split(' ')[1]))
                        phase.append(row.split(' ')[0])
                if f'{Phase}Phase' in phase and 'DISPhase' in phase:
                    # print(directory)
                    BCCindex = phase.index(f'{Phase}Phase')
                    DISindex = phase.index('DISPhase')
                    freeenergy_diff = -(freeenergy[BCCindex]-freeenergy[DISindex])
                    if freeenergy_diff<tol:
                        store_dict_low[chiab,ratioab].append(directory)
                    
                    if freeenergy_diff>tol:
                        seeddict[chiab,ratioab].append(directory)
        export_dict = dict()
        for i in chi:
            for j in ratio:
                export_dict[i,j] = []
                for direct in store_dict_low[i,j]:
                    # print(direct)
                    path1 = os.path.join(pwd,direct)
                    if len(seeddict[i,j])>1:
                        fA = float(direct.split('/')[-1][2:])
                        fAlist = []
                        for seed in seeddict[i,j]:
                            fAlist.append(float(seed.split('/')[-1][2:]))
                        diff_fA = (np.array(fAlist)-fA)**2
                        seedloc = np.where(np.min(diff_fA)==diff_fA)[0][0]
                        path2 = os.path.join(pwd,seeddict[i,j][seedloc])
                        
                    elif len(seeddict[i,j])==1:
                        path2  = os.path.join(pwd,seeddict[i,j][0])
                    else:
                        continue
                    oe.write(f'{path1} {path2} {Phase} \n')             
    oe.close() 