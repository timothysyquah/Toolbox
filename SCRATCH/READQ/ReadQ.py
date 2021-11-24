#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 12:41:30 2021

@author: tquah
"""

import numpy as np
import glob 
import re 
import argparse


# filelist = glob.glob('./**/STDOUT',recursive=True)
# file = filelist[0]
def CreateComplexField(array):
    return array[:,0]+1j*array[:,1]


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dirs', action='store',required=True, nargs='+',help='list of directories that contain each phase point')
    parser.add_argument('--ignorephase', action='store',default=[''], nargs='+',help='phases to ignore')
    parser.add_argument('--ignorestatus', action='store', default=[1,3], nargs='+',help='status to ignore')
    parser.add_argument('-f', '--filename', action='store', default='STDOUT',help='file that contains the phases and free energies at each phase point')
    
    args = parser.parse_args()
    print(args)

    filename=args.filename

    
    
    op = open(file)
    lines = op.readlines()
    op.close()
    Qlist = []
    for line in lines:
        if 'Partition Function' in line:
            Qreal=float(line.split()[4])
            Qim=float(line.split()[6])
            Qlist.append([Qreal,Qim])
    
    Qarray = np.vstack(Qlist)
    Qarray = CreateComplexField(Qarray)