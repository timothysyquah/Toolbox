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

# os.chdir('/home/tquah/Projects/sweep-asym-armlength_corrected/')
# print(os.listdir())
os.chdir("/home/tquah/Projects/DMREF/sweep-asym-armlength_corrected/")
dirs = glob.glob("/home/tquah/Projects/DMREF/sweep-asym-armlength_corrected/chiAB_0.0289/Nsc*/fA*")
pwd = os.getcwd()
phaselist = ['FCC','BCC']

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool to compute phase boundaries')
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default=glob.glob("tau*/phiA*"),help='list of directories that contain each phase point')
    parser.add_argument('-o', '--outfig', action='store', default='',help='name of output figure file')
    parser.add_argument('--raw', action='store', default='',help='name of raw output file (for plotting in another program')
    parser.add_argument('-t', '--plottype', action='store', default='simplecolors',help='type of plot to generate')
    parser.add_argument('--xlabel', action='store', default=r"$f_A$",help='label for xaxis, can use \'$\' to write latex')
    parser.add_argument('--zlabel', action='store', default=r"$\tau$")
    parser.add_argument('--ylabel', action='store', default=r"$\chi N$",help='')
    parser.add_argument('--axisrange', action='store', nargs=4, default=[None,None,None,None],help='')
    parser.add_argument('--linecutoff', action='store', nargs='+', default=[0.2,1],help='maximum length of lines to draw in phase diagrams, useful to clean them up')
    parser.add_argument('-n','--dim',action='store',default=None,help='Number of dimensions to plot phase data in \n   1 => Free energy curves\n   2 => Phase Diagram \n   3 => 3d phase diagram  (guesses by default)')
    parser.add_argument('-i','--interp_dimension',action='store',default=[0],nargs='+',help='Dimensions to interpolate the phase diagram along ex: [0,1] would interpolate in 2 dimensions')
    parser.add_argument('-p','--plotstyle3d',action='store',default='flat',help='This argument changes the 3d plot style. Flat => multiple graphs with different linestyles on top of each other')
    parser.add_argument('--stylesheet',action= 'store',default=os.path.dirname(os.path.realpath(sys.argv[0]))+'/better_style.mplstyle',help='This argument is the Matplotlib stylesheet that will be used for graphing') #the default is located in the directory this script is located at
    parser.add_argument('--aspect',action='store',default=None,help='The aspect ratio for the outputted figure use 1 for a square fig, works for a 2d graph right now')
    parser.add_argument('-r', '--refphase', action='store', default=None,help='name of phase to reference to, only matters if 1d')
    parser.add_argument('-k', '--keywrd', action='store', default=[], nargs='+', help='axis to plot',type=str)
    print("IMPLEMENT CUSTOM AXIS RANGES AND LABELS FROM COMMAND LINE")
    args = parser.parse_args()

seedlist = []
exportpath = os.path.join(pwd,'TEST.dat')

oe = open(exportpath,'w+')

for Phase in phaselist:
    tol = 1e-5
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
    
    
    # store_dict_high= dict()
    store_dict_low= dict()
    seeddict = dict()
    for i in chi:
        for j in ratio:
            # store_dict_high[i,j] = [1,1,'','']
            store_dict_low[i,j] = []
            seeddict[i,j] = []
        
        
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
                    # print(directory)
                    # print(freeenergy_diff)
                    store_dict_low[chiab,ratioab].append(directory)
                
                if freeenergy_diff>tol:
                    # print(directory)
                    # print(freeenergy_diff)
                    seeddict[chiab,ratioab].append(directory)
                
               
     
    export_dict = dict()
    for i in chi:
        for j in ratio:
            export_dict[i,j] = []
            for direct in store_dict_low[i,j]:
                # print(direct)
                path1 = direct

                if len(seeddict[i,j])>1:
                    fA = float(direct.split('/')[-1][2:])
                    fAlist = []
                    for seed in seeddict[i,j]:
                        fAlist.append(float(seed.split('/')[-1][2:]))
                    diff_fA = (np.array(fAlist)-fA)**2
                    seedloc = np.where(np.min(diff_fA)==diff_fA)[0][0]
                    path2 = seeddict[i,j][seedloc]
                    
                
                else:
                    path2  = seeddict[i,j][0]

                
                
                oe.write(f'{path1} {path2} {Phase}Phase \n')

             
oe.close() 
        
 
                
                
                
    #             if freeenergy_diff<tol:
        
    #                 ODT_directories.append(directory)
                    
    #                 if fab>0.5:
    #                     if fab<store_dict_high[chiab,ratioab][0]:
    #                         store_dict_high[chiab,ratioab][0] = fab
    #                         store_dict_high[chiab,ratioab][1] = freeenergy_diff
    #                         store_dict_high[chiab,ratioab][2] = directory
    #                 if fab<0.5:
        
                        
    #                     if fab>store_dict_low[chiab,ratioab][0]:
    #                         store_dict_low[chiab,ratioab][0] = fab
    #                         store_dict_low[chiab,ratioab][1] = freeenergy_diff
    #                         store_dict_low[chiab,ratioab][2] = directory
    #         os.chdir(pwd)
    
    #     else:
    #         # print(directory)
    #         filelist = os.listdir()
    #         if len(filelist)==0:
    #             print('Deleting Empty Directory')
    #             print(directory)
    #             shutil.rmtree(directory)
            
    #         os.chdir(pwd)
            
        
    
    # for i in chi:
    #     for j in ratio:
            
            
    #         def closest_fA(storelist,fAlist):
    #             fAarray = np.array(fAlist)
    #             if storelist[0]>0.5:
    #                 loc = np.where(np.array(fAlist)<storelist[0])[0]
                    
    #             if storelist[0]<0.5:
    #                 loc = np.where(np.array(fAlist)>storelist[0])[0]
    #             # print(loc)
    #             fAarray = fAarray[loc]
                
    #             fAsub = (np.array(storelist[0])-fAarray)**2
                
    #             # nonzerolist = np.nonzero(fAsub)
    #             # array_of_fA = fAsub[nonzerolist]
    #             # print(array_of_fA)
    #             fAloc = np.where(fAsub==np.min(fAsub))[0]
    #             # print(fAloc)
    #             # print(fAloc)
    #             # print(fAarray[fAloc][0])
    #             return fAarray[fAloc]
                
                
            
    #         fA_seed = closest_fA(store_dict_high[i,j],fAdict[i,j])
    #         # print(fA_seed)
    #         partialpath = f'chiAB_{i:0.4f}/NscA_{j[0]:0.1f}_NscB_{j[1]:0.1f}/fA{fA_seed[0]:0.5f}'
            
            
            
            
            
    #         os.chdir(partialpath)
    #         op = open('F0_phases.dat')
    #         rows = op.read().split('\n')
    #         op.close()
            
    #         fAdict[i,j].append(fab)
    #         freeenergy = []
    #         phase = []
    #         for row in rows:
    #             if len(row)!=0:
    #                 freeenergy.append(float(row.split(' ')[1]))
    #                 phase.append(row.split(' ')[0])
                    
            
    #         if f'{Phase}Phase' in phase and 'DISPhase' in phase:
    #             # print(directory)
    #             BCCindex = phase.index(f'{Phase}Phase')
    #             DISindex = phase.index('DISPhase')
    #             freeenergy_diff = -(freeenergy[BCCindex]-freeenergy[DISindex])
    #             if freeenergy_diff<tol:
    #                 print('Warning Seed Could be DIS')
    #         else: 
    #             print(f'No {phase}')
                
                
    
            
    #         os.chdir(pwd)
            
            
    #         store_dict_high[i,j][-1] = os.path.join(pwd,partialpath)
            
    #         fA_seed = closest_fA(store_dict_low[i,j],fAdict[i,j])
    
    #         partialpath = f'chiAB_{i:0.4f}/NscA_{j[0]:0.1f}_NscB_{j[1]:0.1f}/fA{fA_seed[0]:0.5f}'
    #         os.chdir(partialpath)
    #         op = open('F0_phases.dat')
    #         rows = op.read().split('\n')
    #         op.close()
            
    #         fAdict[i,j].append(fab)
    #         freeenergy = []
    #         phase = []
    #         for row in rows:
    #             if len(row)!=0:
    #                 freeenergy.append(float(row.split(' ')[1]))
    #                 phase.append(row.split(' ')[0])
                    
            
    #         if f'{Phase}Phase' in phase and 'DISPhase' in phase:
    #             # print(directory)
    #             BCCindex = phase.index(f'{Phase}Phase')
    #             DISindex = phase.index('DISPhase')
    #             freeenergy_diff = -(freeenergy[BCCindex]-freeenergy[DISindex])
    #             if freeenergy_diff<tol:
    #                 print('Warning Seed Could be DIS')
    #         else: 
    #             print('No BCC')
                
                
    
            
    #         os.chdir(pwd)
            
    
    #         store_dict_low[i,j][-1] = os.path.join(pwd,partialpath)
    
    # exportpath = os.path.join(pwd,'TEST.dat')
    
    # op = open(exportpath,'w+')
    
    # for i in chi:
    #     for j in ratio:
    #         # path1 = os.path.join(store_dict_high[i,j][-2],'BCCPhase')
    #         # path2 = os.path.join(store_dict_high[i,j][-1],'BCCPhase')
    #         op.write(f'{store_dict_high[i,j][-2]} {store_dict_high[i,j][-1]} \n')
    #         # path1 = os.path.join(store_dict_low[i,j][-2],'BCCPhase')
    #         # path2 = os.path.join(store_dict_low[i,j][-1],'BCCPhase')
    #         op.write(f'{store_dict_low[i,j][-2]} {store_dict_low[i,j][-1]} \n')
    # op.close()
    
    
