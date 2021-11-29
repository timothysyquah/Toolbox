#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 23 12:41:30 2021

@author: tquah
"""

import numpy as np
import glob 
import argparse
import os

def ReadQ(filepath):
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
    return Qarray
def checkStatus(wdir,status2ignore):
    filename="{}/STATUS".format(wdir)
    try:
        with open(filename,'r') as f:
            status=int(f.readline())
    except FileNotFoundError as e:
        print (e)
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


def GetWDIR(path):
    split_path = path.split('/')
    WDIR = ''
    
    if path[0]!='/':
        WDIR+=split_path[0]
        WDIR+='/'
    else:
        WDIR+='/'
        WDIR+=split_path[0]
        WDIR+='/'

    for i in range(1,len(split_path)-1):
        WDIR+=split_path[i]
        WDIR+='/'
    return WDIR

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--ignorestatus', action='store', default=[], nargs='+',help='status to ignore')
    parser.add_argument('-d', '--dirs', action='store', default='./**',help='Base directory structure')
    parser.add_argument('-f', '--filename', action='store', default='STDOUT',help='file that contains the phases and free energies at each phase point')
    parser.add_argument('-e', '--exportname', action='store', default='Q.dat',help='file that contains the phases and free energies at each phase point')
    args = parser.parse_args()
    filename=args.filename
    dirs = args.dirs
    fullpath = os.path.join(dirs,filename)
    filelist = glob.glob(fullpath,recursive = True)
    for file in filelist:
        wdir_ = GetWDIR(file)
        status = checkStatus(wdir_,args.ignorestatus)
        if status:
            Qtemp = ReadQ(file)
            exportpath = os.path.join(wdir_,args.exportname)
            header = 'Q.Real Q.Imag'
            np.savetxt(exportpath,Qtemp,header = header)
