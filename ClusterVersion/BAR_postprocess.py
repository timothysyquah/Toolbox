#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep  7 14:37:11 2021

@author: tquah
"""

import numpy as np
import sys
from scipy.optimize import minimize
import argparse as ap
import os
sys.path.append('/home/tquah/PolyFTS_ALL/PolyFTS/tools/')
from stats import extractData,autoWarmupMSER

def gx(x):
    return (1+np.exp(x))**(-1)
def mse(x):
    return x**2

class CostStuff:
    def __init__(self):
        pass
    def CostFunction(self,x):
        diff = np.mean(gx(self.H_AB+x))-np.mean(gx(self.H_BA-x))
        return mse(diff)
if __name__ == "__main__":
    parser = ap.ArgumentParser(description='Statistical analysis of PolyFTS data')
    parser.add_argument('-f','--file',default='./dH.dat',type=ap.FileType('r'),help='Filename containing scalar statistical data.')
    parser.add_argument('-tol','--tolerance',default = 1e-6,type = float,help = 'Tolerance for Optimizer')
    args=parser.parse_args()
  
# in dH file dH model0 Re Im | dH model1 Re Im
#ie. Delta H_{BA} = H_B(w^A)-H_A(w^A)
#ie. Delta H_{AB} = H_A(w^B)-H_B(w^B)
    columns = [1,2,3,4]
    
    prodlist = []
    for col in columns:
        warmup,proddata,nwarmup = autoWarmupMSER(args.file,col)
        prodlist.append(proddata)
    Refun = CostStuff()
    Refun.H_AB = prodlist[0]
    Refun.H_BA = prodlist[2]
    Re_dict = minimize(Refun.CostFunction,x0 = 1,tol = args.tolerance)
    Imfun = CostStuff()
    Imfun.H_AB = prodlist[1]
    Imfun.H_BA = prodlist[3]
    Im_dict = minimize(Imfun.CostFunction,x0 = 1,tol = args.tolerance)
    
    if Re_dict.status+Im_dict.status==0:
        print(f'Free Energy: {Re_dict.x[0]} (Re) {Im_dict.x[0]} (Im)' )
    else:
         print('Optimizer Problems' )
    
