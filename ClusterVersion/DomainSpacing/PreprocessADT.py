#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 13 14:41:22 2021

@author: tquah
"""
import numpy as np
import pandas as pd
import os
import argparse
if __name__ == '__main__':
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to weight blocks')
    parser.add_argument('-of','--operator_file',action = 'store',default = 'operators.dat', help = 'File with operator output', type = str)
    parser.add_argument('-af','--adt_file',action = 'store',default = 'ADT_tsteps.dat', help = 'File to read with averages and error', type = str)
    parser.add_argument('-eo','--export_operator_file',action = 'store',default = 'RW_operators.dat', help = 'Reweighted operators.dat file', type = str)
    parser.add_argument('-ea','--export_adt_file',action = 'store',default = 'BlockSteps.dat', help = 'Block steps file', type = str)
    args = parser.parse_args()
    #load data in 
    f = open(args.operator_file, 'r')
    line1 = f.readline()[1:]
    df = pd.read_csv(f, sep='\s+', names=line1.replace(' #', '').split(), dtype=np.float)
    dtarray = np.loadtxt(args.adt_file)
    operator_array = np.array(df)
    #figure out how big blocks are
    blockssize = int(df['step'][1])
    #determine number of blocks
    shape = np.shape(dtarray)
    numofblocks = int(shape[0]/blockssize)
    cutoff = numofblocks*blockssize
    reshape_array = dtarray[0:cutoff,0].reshape(numofblocks,blockssize)
    blocktime = np.sum(reshape_array, axis=1).reshape(1,numofblocks)
    newoperators = blocktime.transpose()*operator_array[1:,1:]
    operator_array[1:,1:] = newoperators
    header = list(df)
    string = ''
    for i in range(len(header)):
        string+=header[i]+' '
    np.savetxt(args.export_operator_file,operator_array,header = string)
    np.savetxt(args.export_adt_file,blocktime)