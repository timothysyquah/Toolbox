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
from collections import defaultdict
import pickle
if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Tool to compute phase boundaries')
    parser.add_argument('-f', '--filename', action='store', default='F0_phases.dat',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-k', '--keyword', action='store', nargs='+', default=['chi','f'],help='First to Second to last will be written into nested dictionaries Last position will be written in array',type=str)
    parser.add_argument('-e', '--export_file_name', action='store', default='data.dict',help='name export pickle file')
    args = parser.parse_args()
    IDIR = os.getcwd()
    dir_dict = dict()
    for file in glob.glob(f'**/{args.filename}', recursive = True):
        dir_split = file.split('/')
        full_path = os.path.join(IDIR,file)
        op = open(full_path,'r')
        dataread = op.read().splitlines()
        op.close()
        dictionary_list = []
        for i in range(0,len(args.keyword)-1):
            dir1loc = dir_split.index([s for s in dir_split if args.keyword[i] in s][0])
            dir1=dir_split[dir1loc]
            dictionary_list.append(float(re.findall("\d+\.\d+", dir1)[0]))

        dir2loc = dir_split.index([s for s in dir_split if args.keyword[-1] in s][0])
        dir2=dir_split[dir2loc]
        innerstore = float(re.findall("\d+\.\d+", dir2)[0])
        
        for i in range(0,len(dataread),1):
            newlist = []
            splitdata = dataread[i].split(' ')
            phase = splitdata[0][:-5]
            newlist+=dictionary_list
            newlist.append(phase)      
            if int(splitdata[2])==2:
                if tuple(newlist) in dir_dict:
                    temparray = np.array([innerstore,float(splitdata[1])]) 
                    dir_dict[tuple(newlist)] = np.vstack((dir_dict[tuple(newlist)],temparray))
                else:
                    dir_dict[tuple(newlist)] = np.array([innerstore,float(splitdata[1])])
            else:
                print('Simulations have STATUS 0 1 3...revist extractF0.dat...')
                break
    with open(args.export_file_name, 'wb') as outfile:
        pickle.dump(dir_dict, outfile, protocol=pickle.HIGHEST_PROTOCOL)