#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 12:15:35 2021

@author: tquah
"""

import numpy as np
import glob
import argparse
import os
import re
from CL_Domain_Spacing import *
def OpenGetData(file):
    op = open(file,'r')
    data = op.read()
    op.close()
    return data

def extract_value(string):
    return re.findall(r"[-+]?\d*\.\d+|\d+", string) #im lazy so i made a function 


def ExtractDomain(fullpath):
    data = OpenGetData(fullpath)
    datasplit = data.splitlines()
    try:
        for i in range(len(datasplit)):
            if "Final simulation cell:" in datasplit[i]:
                values = extract_value(datasplit[i])
                D_0 = float(values[1])*10**float(values[2])
        os.chdir(IDIR)
        return D_0
    except:
        return 0


if __name__ == '__main__':
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to sweep bottlebrushes')
    parser.add_argument('-d', '--dirs', action='store', nargs='+', default="SCFT/**/chi0.2/nsc10.0/nbb*/**/",help='list of directories that contain each phase point')
    parser.add_argument('-f','--file',action = 'store',default = 'LAM3D.out', help = 'File to read with averages and error', type = str)
    parser.add_argument('-e','--exportname',action = 'store',default = 'SCFT_domain_data.dat', help = 'Path to Save Data', type = str)
    parser.add_argument('-dp', '--desired_parameters', action='store',nargs='+', default=['chi','nsc','nbb'],help='Desired parameters',type = str)
    parser.add_argument('-o', '--ordering', action='store', default='nbb',help='Ordering based on parameter',type = str)

    args = parser.parse_args()
    
    path = os.path.join(args.dirs,args.file)
    file_list = glob.glob(path,recursive=True)
    arraylist = []
    for i in range(len(file_list)):
        parray = ParsePath(file_list[i],args.desired_parameters)
        D0 = ExtractDomain(file_list[i])
        temparray = np.zeros((1,len(args.desired_parameters)+1))
        temparray[0,0:-1] = parray
        temparray[0,-1] = D0
        if D0!=0:
            arraylist.append(temparray)
    data_array = np.vstack(arraylist)
    data_array[:,args.desired_parameters.index(args.ordering)] = data_array[:,args.desired_parameters.index(args.ordering)]+1
    reorder = np.argsort(data_array[:,args.desired_parameters.index(args.ordering)])
    data_array = data_array[reorder,:]
    np.savetxt(args.exportname,data_array)