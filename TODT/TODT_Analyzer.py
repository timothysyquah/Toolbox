#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 16:15:41 2020

@author: tquah
"""

import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import fsolve
from joblib import Parallel, delayed
import multiprocessing



def fzeroCubicSpline(x0,fun):
    return fsolve(fun,x0)



path = '/home/tquah/IMPORT_KNOT/TODT.dict'

number_of_initial_guesses = 10




with open(path, 'rb') as handle:
    F0dat = pickle.load(handle)


overall_list = []
for ODTpt in F0dat:
    
    if ODTpt[2]=='LAM':
        overall_list.append(ODTpt)
        

for i in range(0,len(overall_list)):
    
    Larray = F0dat[overall_list[i][0],overall_list[i][1],'LAM'][F0dat[overall_list[i][0],overall_list[i][1],'LAM'][:,0].argsort()]
    Darray = F0dat[overall_list[i][0],overall_list[i][1],'DIS'][F0dat[overall_list[i][0],overall_list[i][1],'DIS'][:,0].argsort()]
    Neff = (overall_list[i][0]+1)*overall_list[i][1]
    x0 = np.max(Larray[:,0])
    Diff = CubicSpline(Larray[:,0],Larray[:,1]-Darray[:,1])    
    Derriv = Diff.derivative()

    
    # plt.plot(Larray[:,0],Larray[:,1])
    # plt.plot(Darray[:,0],Darray[:,1])
    # plt.plot(Darray[:,0],Larray[:,1]-Darray[:,1])
    xrun = np.linspace(np.min(Larray[:,0]),np.max(Larray[:,0]),100)
    ydiff = Diff(xrun)
    dydiff = Derriv(xrun)

    
    # plt.plot(xrun,dydiff)
    mean = np.mean(dydiff)
    std = np.std(dydiff)
    
    if mean==0 and std==0:
        print('Skipping...ODT not found')
    else:
        print('Working')
        fig, ax1 = plt.subplots()
        ax2 = ax1.twinx()
        ax1.plot(xrun,ydiff)
        ax2.plot(xrun,dydiff)

    
    
    
    