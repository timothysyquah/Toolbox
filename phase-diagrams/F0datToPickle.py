#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jul 26 22:54:27 2020

@author: tquah
"""

import numpy as np
import argparse
import glob,re
import os,sys
import matplotlib.pyplot as plt
from phase_data_plotter import *

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool to compute phase boundaries')
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default=glob.glob("tau*/phiA*"),help='list of directories that contain each phase point')
    parser.add_argument('-k', '--keyword', action='store', nargs='+', default=['chiAB_','fA'],help='First postition')

    # parser.add_argument('-o', '--outfig', action='store', default='',help='name of output figure file')
    # parser.add_argument('--raw', action='store', default='',help='name of raw output file (for plotting in another program')
    # parser.add_argument('-t', '--plottype', action='store', default='simplecolors',help='type of plot to generate')
    # parser.add_argument('--xlabel', action='store', default=r"$f_A$",help='label for xaxis, can use \'$\' to write latex')
    # parser.add_argument('--zlabel', action='store', default=r"$\tau$")
    # parser.add_argument('--ylabel', action='store', default=r"$\chi N$",help='')
    # parser.add_argument('--axisrange', action='store', nargs=4, default=[None,None,None,None],help='')
    # parser.add_argument('--linecutoff', action='store', nargs='+', default=[1e30,1e30],help='maximum length of lines to draw in phase diagrams, useful to clean them up')
    # parser.add_argument('-n','--dim',action='store',default=None,help='Number of dimensions to plot phase data in \n   1 => Free energy curves\n   2 => Phase Diagram \n   3 => 3d phase diagram  (guesses by default)')
    # parser.add_argument('-i','--interp_dimension',action='store',default=[0],nargs='+',help='Dimensions to interpolate the phase diagram along ex: [0,1] would interpolate in 2 dimensions')
    # parser.add_argument('-p','--plotstyle3d',action='store',default='flat',help='This argument changes the 3d plot style. Flat => multiple graphs with different linestyles on top of each other')
    # parser.add_argument('--stylesheet',action= 'store',default=os.path.dirname(os.path.realpath(sys.argv[0]))+'/better_style.mplstyle',help='This argument is the Matplotlib stylesheet that will be used for graphing') #the default is located in the directory this script is located at
    # parser.add_argument('--aspect',action='store',default=None,help='The aspect ratio for the outputted figure use 1 for a square fig, works for a 2d graph right now')
    # parser.add_argument('-r', '--refphase', action='store', default=None,help='name of phase to reference to, only matters if 1d')
    # print("IMPLEMENT CUSTOM AXIS RANGES AND LABELS FROM COMMAND LINE")
    args = parser.parse_args()
    
    args.dir = '/home/tquah/IMPORT_BRAID/diblock_phasediagram'
    os.chdir(args.dir)
    # print(os.getcwd())
    # print(os.listdir(DIR))
    # for file in glob.glob(args.filename):
    #     print(file)
    
    file_list = []
    dir_dict = dict()
    count = 0
    for file in glob.glob(f'**/{args.filename}', recursive = True):
        file_list.append(file)
        
        dir_split = file.split('/')
        
        for i in range(0,len(dir_split)):
            
            if count==0:
                dir_dict[i] = []
            
            dir_dict[i].append(dir_split[i])
        count+=1
    
    levels = []
    unique_dir = []
    for i in range(0,len(dir_dict),1):
        unique_dir.append(sorted(list(set(dir_dict[i]))))
        levels.append(len(list(set(dir_dict[i]))))
        
        
        
        
        
        # dir_dict[]
        
        
        # print(file.split('/'))

    
    # os.chdir("/mydir")
    # for file in glob.glob("*.txt"):
    #     print(file)
