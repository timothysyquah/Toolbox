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
from scipy.optimize import fsolve,curve_fit
from joblib import Parallel, delayed
import multiprocessing



def fzeroCubicSpline(fun,x0):
    return fsolve(fun,x0)



path = '/home/tquah/IMPORT_KNOT/TODT.dict'

number_of_initial_guesses = 100
plotting = False



with open(path, 'rb') as handle:
    F0dat = pickle.load(handle)


overall_list = []
for ODTpt in F0dat:
    
    if ODTpt[2]=='LAM':
        overall_list.append(ODTpt)
        
# chi = np.zeros()
def curvefitfun(x,a,b,c):
    return a*x**b+c


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
        # print('Working')
        initialguess_array = np.linspace(np.min(Larray[:,0]),np.max(Larray[:,0]),number_of_initial_guesses)
        zeroarray = Parallel(n_jobs=-1)(delayed(fzeroCubicSpline)(Diff,initialguess_array[i]) for i in range (0,number_of_initial_guesses,1))
        chi = np.max(zeroarray)
        # print(chi)
        
        
        
        if plotting:
            plt.figure()
            plt.plot(xrun,ydiff,'r')
            plt.plot(chi,Diff(chi),'ok')
            # width = (np.max(ydiff)-np.min(ydiff))*0.025
            # plt.ylim(Diff(chi)-width,Diff(chi)+width)
        chiN = chi*Neff
        if i==0:
            DataArray = np.array([overall_list[i][0],overall_list[i][1],chiN])
        else:
            DataArray = np.vstack((DataArray,[overall_list[i][0],overall_list[i][1],chiN]))

#plot alpha

plt.figure()
alpha = DataArray[:,0]/DataArray[:,1]
plt.plot(alpha,DataArray[:,2],'ok')





marker = ['+', 'o', '*','v','^','<','>','s','p','h','H','x']
color = ['k','b', 'r', 'c', 'm', 'y','g', 'tab:red','tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

plt.figure()
nbbunique = sorted(set(list(DataArray[:,1])))
alpha = DataArray[:,0]/DataArray[:,1]
count = 0
for nbb in nbbunique:
    loc = np.where(DataArray[:,1]==nbb)
    
    plt.scatter(alpha[loc],DataArray[loc,2],color=color[count],marker=marker[count],label=f'Nbb={int(nbb+1)}')
    count+=1
plt.xlabel(r'$\alpha$ ($N_{sc}/N_{BB}$)')
plt.ylabel(r'$\chi N_{eff}$')

plt.legend()


plt.figure()
nbbunique = sorted(set(list(DataArray[:,1])))
alpha = DataArray[:,0]/DataArray[:,1]
count = 0
for nbb in nbbunique:
    loc = np.where(DataArray[:,1]==nbb)
    
    plt.scatter(DataArray[loc,0],DataArray[loc,2],color=color[count],marker=marker[count],label=f'Nbb={int(nbb+1)}')
    count+=1
plt.xlabel(r'$N_{sc}$')
plt.ylabel(r'$\chi N_{eff}$')

plt.legend()


# nscunique = sorted(set(list(DataArray[:,0])))
# plt.figure()
# count = 0
# for nsc in nscunique:
#     loc = np.where(DataArray[:,0]==nsc)
    
#     plt.scatter(alpha[loc],DataArray[loc,2],color=color[count],marker=marker[count],label=f'Nsc={nsc}')
#     count+=1
# plt.legend()


nscunique = sorted(set(list(DataArray[:,0])))
plt.figure()
count = 0
for nsc in nscunique:
    loc = np.where(DataArray[:,0]==nsc)
    
    plt.scatter(DataArray[loc,1],DataArray[loc,2],color=color[count],marker=marker[count],label=f'Nsc={nsc}')
    count+=1
plt.legend()
plt.xlabel(r'$N_{bb}$')
plt.ylabel(r'$\chi N_{eff}$')


# plt.plot(alpha[9:16],DataArray[9:16,2],'or',label='Nbb=100')
# plt.plot(alpha[17:23],DataArray[17:23,2],'ob',label='Nbb=150')
# plt.plot(alpha[24:],DataArray[24:,2],'og',label='Nbb=200')
# plt.legend()

# plt.figure()
# alpha = DataArray[:,0]/DataArray[:,1]
# plt.plot(alpha[0:8],DataArray[0:8,2],'ok',label='Nbb=50')
# plt.plot(alpha[9:16],DataArray[9:16,2],'or',label='Nbb=100')
# plt.plot(alpha[17:23],DataArray[17:23,2],'ob',label='Nbb=150')
# plt.plot(alpha[24:],DataArray[24:,2],'og',label='Nbb=200')
# plt.legend()

plt.figure()
nbbunique = sorted(set(list(DataArray[:,1])))
alpha = DataArray[:,0]/DataArray[:,1]
count = 0
for nbb in nbbunique:
    loc = np.where(DataArray[:,1]==nbb)
    
    alpha = DataArray[:,0]/DataArray[:,1]
    plt.plot(alpha,DataArray[:,2],'ok')

    plt.scatter(alpha[loc],DataArray[loc,2]*1.1,color='r',marker=marker[count],label=f'Nbb={int(nbb+1)}')
    plt.scatter(alpha[loc],DataArray[loc,2]*0.90,color='r',marker=marker[count],label=f'Nbb={int(nbb+1)}')

    param1 = curve_fit(curvefitfun,alpha[loc],DataArray[loc,2][0]*1.1,p0 =[1,0.5,10])
    x = np.linspace(0,1,100)
    y = curvefitfun(x,param1[0][0],param1[0][1],param1[0][2])
    plt.plot(x,y,'--k')

    
    param2 = curve_fit(curvefitfun,alpha[loc],DataArray[loc,2][0]*0.9,p0 =[1,0.5,10])
    x = np.linspace(0,1,100)
    y = curvefitfun(x,param2[0][0],param2[0][1],param2[0][2])

    plt.plot(x,y,'--k')

    count+=1
    if count==1:
        break
plt.xlabel(r'$N_{sc}$')
plt.ylabel(r'$\chi N_{eff}$')

plt.legend()

            
        