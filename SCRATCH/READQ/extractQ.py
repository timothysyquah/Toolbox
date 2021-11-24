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
def checkStatus(wdir,status2ignore):
    filename="{}/STATUS".format(wdir)
    try:
        with open(filename,'r') as f:
            status=int(f.readline())
    except FileNotFoundError as e:
        #print (e)
        print("{} not found...skipping!".format(filename))
        return False

    if status == 0 and 0 in status2ignore:
        print("{} is not converged (killed while running) ...skipping!".format(wdir))
        return False
    if status == 1 and 1 in status2ignore:
        print("{} is divergent ...skipping!".format(wdir))
        return False
    if status == 3 and 3 in status2ignore:
        print("{} is not converged (reached max steps) ...skipping!".format(wdir))
        return False
    if status in status2ignore:
        return False

    return True


def CreateComplexField(array):
    return array[:,0]+1j*array[:,1]

def getStatus(wdir):
    filename="{}/STATUS".format(wdir)
    #try:
    with open(filename,'r') as f:
        status=int(f.readline())
    return status

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-d', '--dirs', action='store',default = './',required=True, nargs='+',help='list of directories that contain each phase point')
    parser.add_argument('--ignorestatus', action='store', default=[1,3], nargs='+',help='status to ignore')
    parser.add_argument('-f', '--filename', action='store', default='STDOUT',help='file that contains the phases and free energies at each phase point')
    
    args = parser.parse_args()
    print(args)

    filename=args.filename
    dirs=args.dirs
    dirs.sort()

    
    
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
    np.savetxt('Q.dat',Qarray)