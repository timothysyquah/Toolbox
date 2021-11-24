#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 23 11:21:17 2021

@author: tquah
"""


import numpy as np
import os
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import argparse
plt.close('all')
plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
# rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'
rcParams['font.style'] = 'normal'
rcParams['axes.labelsize'] = 20
rcParams['xtick.labelsize'] = 15
rcParams['ytick.labelsize'] = 15
rcParams['legend.fontsize'] = 15
rcParams['axes.titlesize'] = 20
def powlaw(x,a,b):
    return np.power(x,b)*a
def powlaw2(x,a):
    b = 2/3
    return np.power(x,b)*a
if __name__ == '__main__':
    IDIR = os.getcwd()
    parser = argparse.ArgumentParser(description='Tool to weight blocks')
    parser.add_argument('-d','--datafile',action = 'store',default = 'domain_data.dat', help = 'File with operator output', type = str)
    parser.add_argument('-p','--points',action = 'store',default = 2, help = 'Data Points to include in powerlaw', type = int)
    parser.add_argument('-g','--graph',action='count',help='Output Graph')
    parser.add_argument('-gn','--graphname',action = 'store',default = 'DomainSpacing.pdf', help = 'Name of Graph', type = str)
    parser.add_argument('--column', action='store', default=[2,3], nargs='+',help='Columns to plot')
    parser.add_argument('-v','--verbose',action='count',help='Outputs more than usual')
    parser.add_argument('-SCFT','--SCFT',action='count',help='Include SCFT?')
    parser.add_argument('-Sd','--SCFTdatafile',action = 'store',default = 'SCFT_domain_data.dat', help = 'File with operator output', type = str)

    args = parser.parse_args()
    dataarray = np.loadtxt(args.datafile)
    var,pcov = curve_fit(powlaw,dataarray[-args.points:,args.column[0]],dataarray[-args.points:,args.column[1]])
    if args.verbose:
        print('CL Power Law using Equation D = a*N^b:')
        print(f'a = {var[0]} and b = {var[1]}')
    if args.SCFT: 
        SCFT_dataarray = np.loadtxt(args.SCFTdatafile)
        
        intersect_value = np.intersect1d(SCFT_dataarray[:,args.column[0]],dataarray[:,args.column[0]])
        arraykeeploc = []
        for i in range(len(intersect_value)):
            arraykeeploc.append(np.where(intersect_value[i]==SCFT_dataarray[:,args.column[0]])[0][0])
        SCFT_dataarray = SCFT_dataarray[arraykeeploc,:]
        varSCFT,pcovSCFT = curve_fit(powlaw,SCFT_dataarray[-args.points:,args.column[0]],SCFT_dataarray[-args.points:,args.column[1]])
            
        if args.verbose:
            print('SCFT Power Law using Equation D = a*N^b:')
            print(f'a = {varSCFT[0]} and b = {varSCFT[1]}')
    if args.graph:
        x = np.linspace(np.min(dataarray[:,args.column[0]]),np.max(dataarray[:,args.column[0]]),100)
        y = powlaw(x,var[0],var[1])
        plt.loglog(x,y,'-r')
        plt.loglog(dataarray[:,args.column[0]],dataarray[:,args.column[1]],\
                   marker = 's',color = 'r',linewidth=0,label = 'FTS-CL')              
        if args.SCFT: 
            xSCFT = np.linspace(np.min(SCFT_dataarray[:,args.column[0]]),np.max(SCFT_dataarray[:,args.column[0]]),100)
            ySCFT = powlaw(xSCFT,varSCFT[0],varSCFT[1])
            plt.loglog(xSCFT,ySCFT,'-b')
            plt.loglog(SCFT_dataarray[:,args.column[0]],SCFT_dataarray[:,args.column[1]],\
                       marker = '^',color = 'b',linewidth=0,label = 'SCFT')              
        plt.legend()
        plt.xlabel('$N_{bb}$')
        plt.ylabel('$D(l)$')
        plt.tight_layout()
        plt.savefig(args.graphname,dpi = 300)   
