#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct  3 23:06:09 2020

@author: tquah
"""

import numpy as np
import os 
import shutil
import glob

# os.chdir('/home/tquah/Projects/sweep-asym-armlength_corrected/')
# print(os.listdir())
os.chdir("/home/tquah/Projects/sweep-asym-armlength_corrected/")
dirs = glob.glob("/home/tquah/Projects/sweep-asym-armlength_corrected/chiAB_0.0289/Nsc*/fA*")
phaselist = ['FCC']

for Phase in phaselist:
    tol = 1e-10
    ODT_directories = []
    ODT_SEED = []
    # fA = []
    ratio = []
    chi = []
    
    for directory in dirs:
        info = directory.split('/')
        # fA.append(info[-1]) 
        ratio.append((float(info[-2].split('_')[1]),float(info[-2].split('_')[3])))
        chi.append(float(info[-3][-6:]))
        
        
    ratio= list(sorted(set(ratio)))
    chi= list(sorted(set(chi)))
    
    pwd = os.getcwd()
    
    store_dict_high= dict()
    store_dict_low= dict()
    fAdict = dict()
    for i in chi:
        for j in ratio:
            store_dict_high[i,j] = [1,1,'','']
            store_dict_low[i,j] = [0,1,'','']
            fAdict[i,j] = []
        
        
    for directory in dirs:
        newpath = os.path.join(pwd,directory)
        os.chdir(newpath)
        split_info = directory.split('/')
        chiab = float(split_info[-3][6:])
        ratioab = (float(split_info[-2].split('_')[1]),float(split_info[-2].split('_')[3]))
        fab = float(split_info[-1][2:])
        # print(os.listdir())
        if 'F0_phases.dat' in os.listdir():
    
            # print(directory)
            op = open('F0_phases.dat')
            rows = op.read().split('\n')
            op.close()
            
            fAdict[chiab,ratioab].append(fab)
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
                # print(freeenergy_diff)
                # print(chiab)
                # print(ratioab)
                if freeenergy_diff<tol:
        
                    ODT_directories.append(directory)
                    
                    if fab>0.5:
                        if fab<store_dict_high[chiab,ratioab][0]:
                            store_dict_high[chiab,ratioab][0] = fab
                            store_dict_high[chiab,ratioab][1] = freeenergy_diff
                            store_dict_high[chiab,ratioab][2] = directory
                    if fab<0.5:
        
                        
                        if fab>store_dict_low[chiab,ratioab][0]:
                            store_dict_low[chiab,ratioab][0] = fab
                            store_dict_low[chiab,ratioab][1] = freeenergy_diff
                            store_dict_low[chiab,ratioab][2] = directory
            os.chdir(pwd)
    
        else:
            # print(directory)
            filelist = os.listdir()
            if len(filelist)==0:
                print('Deleting Empty Directory')
                print(directory)
                shutil.rmtree(directory)
            
            os.chdir(pwd)
            
        
    
    for i in chi:
        for j in ratio:
            
            
            def closest_fA(storelist,fAlist):
                fAarray = np.array(fAlist)
                if storelist[0]>0.5:
                    loc = np.where(np.array(fAlist)<storelist[0])[0]
                    
                if storelist[0]<0.5:
                    loc = np.where(np.array(fAlist)>storelist[0])[0]
                # print(loc)
                fAarray = fAarray[loc]
                
                fAsub = (np.array(storelist[0])-fAarray)**2
                
                # nonzerolist = np.nonzero(fAsub)
                # array_of_fA = fAsub[nonzerolist]
                # print(array_of_fA)
                fAloc = np.where(fAsub==np.min(fAsub))[0]
                # print(fAloc)
                # print(fAloc)
                # print(fAarray[fAloc][0])
                return fAarray[fAloc]
                
                
            
            fA_seed = closest_fA(store_dict_high[i,j],fAdict[i,j])
            # print(fA_seed)
            partialpath = f'chiAB_{i:0.4f}/NscA_{j[0]:0.1f}_NscB_{j[1]:0.1f}/fA{fA_seed[0]:0.5f}'
            
            
            
            
            
            os.chdir(partialpath)
            op = open('F0_phases.dat')
            rows = op.read().split('\n')
            op.close()
            
            fAdict[i,j].append(fab)
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
                    print('Warning Seed Could be DIS')
            else: 
                print(f'No {phase}')
                
                
    
            
            os.chdir(pwd)
            
            
            store_dict_high[i,j][-1] = os.path.join(pwd,partialpath)
            
            fA_seed = closest_fA(store_dict_low[i,j],fAdict[i,j])
    
            partialpath = f'chiAB_{i:0.4f}/NscA_{j[0]:0.1f}_NscB_{j[1]:0.1f}/fA{fA_seed[0]:0.5f}'
            os.chdir(partialpath)
            op = open('F0_phases.dat')
            rows = op.read().split('\n')
            op.close()
            
            fAdict[i,j].append(fab)
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
                    print('Warning Seed Could be DIS')
            else: 
                print('No BCC')
                
                
    
            
            os.chdir(pwd)
            
    
            store_dict_low[i,j][-1] = os.path.join(pwd,partialpath)
    
    exportpath = os.path.join(pwd,'TEST.dat')
    
    op = open(exportpath,'w+')
    
    for i in chi:
        for j in ratio:
            # path1 = os.path.join(store_dict_high[i,j][-2],'BCCPhase')
            # path2 = os.path.join(store_dict_high[i,j][-1],'BCCPhase')
            op.write(f'{store_dict_high[i,j][-2]} {store_dict_high[i,j][-1]} \n')
            # path1 = os.path.join(store_dict_low[i,j][-2],'BCCPhase')
            # path2 = os.path.join(store_dict_low[i,j][-1],'BCCPhase')
            op.write(f'{store_dict_low[i,j][-2]} {store_dict_low[i,j][-1]} \n')
    op.close()
    
    
