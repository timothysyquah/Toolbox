#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 16:15:41 2020

@author: tquah
"""
import os
import numpy as np
import pickle
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline
from scipy.optimize import fsolve,curve_fit
from joblib import Parallel, delayed
import multiprocessing
import random

plt.rc('font', family='serif')
from matplotlib import rcParams
rcParams['text.usetex'] = True 
# rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'

rcParams['axes.labelsize'] =20
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 11

def fzeroCubicSpline(fun,x0):
    return fsolve(fun,x0)


E_Path = '/home/tquah/Presentations/FirstYearTalkQuah/images/'
# path = '/home/tquah/IMPORT_KNOT/TODT.dict'
# path = '/home/tquah/IMPORT_BRAID/T0data.dict'
path = '/home/tquah/IMPORT_KNOT/FJC.dict'

# path = '/home/tquah/IMPORT_KNOT/TODT1.dict'

number_of_initial_guesses = 100
plotting = False
plt.close('all')



with open(path, 'rb') as handle:
    F0dat = pickle.load(handle)


overall_list = []
for ODTpt in F0dat:
    
    if ODTpt[2]=='LAM':
        overall_list.append(ODTpt)
        
# chi = np.zeros()
def curvefitfun(x,a,b):
    c = 10.5
    return a*x**b+c

# def curvefitfun(x,a,b,c):
#     # c = 10.5
#     return a*x**b+c

# def chiN(x,a,b,c):
#     return a+(x**c)*b
scount= 0
tol = 1e-8

op = open('odtfound.dat','w+')


initial_count = 0
for i in range(0,len(overall_list)):
    
    Linitial =   F0dat[overall_list[i][0],overall_list[i][1],'LAM']
    Dinitial = F0dat[overall_list[i][0],overall_list[i][1],'DIS']
    
    if len(np.shape(Linitial))!=2:
        continue
    
    Larray = Linitial[Linitial[:,0].argsort()]
    
    
    
    Darray = Dinitial[Dinitial[:,0].argsort()]
    # plt.figure()
    # plt.plot(Larray[:,0],Larray[:,1]-Darray[:,1])

    Neff = (overall_list[i][0]+1)*overall_list[i][1]
    x0 = np.max(Larray[:,0])
    if len(Larray[:,1])==len(Darray[:,1]):
        Diff = CubicSpline(Larray[:,0],Larray[:,1]-Darray[:,1])    
        Derriv = Diff.derivative()
    
        # plt.figure()
        # plt.plot(Larray[:,0],Larray[:,1])
        # plt.plot(Darray[:,0],Darray[:,1])
        # plt.scatter(Darray[:,0],Larray[:,1]-Darray[:,1])
        xrun = np.linspace(np.min(Larray[:,0]),np.max(Larray[:,0]),100)
        ydiff = Diff(xrun)
        dydiff = Derriv(xrun)
    
        
        # plt.plot(xrun,dydiff)
        mean = np.mean(dydiff)
        std = np.std(dydiff)
        
        if mean<tol and std<tol:
            if scount==0:
                print('Skipping...ODT not found')
                scount+=1
            print('Nbb = %d and Nsc = %d'%(overall_list[i][1],overall_list[i][0]))
        else:
            
            # print('Working')
            initialguess_array = np.linspace(np.min(Larray[:,0]),np.max(Larray[:,0]),number_of_initial_guesses)
            zeroarray = Parallel(n_jobs=-1)(delayed(fzeroCubicSpline)(Diff,initialguess_array[i]) for i in range (0,number_of_initial_guesses,1))
            chi = np.max(zeroarray)
            # print(chi)
            
            
            
            # if plotting:
                # plt.figure()
                # plt.plot(xrun,ydiff,'r')
                # plt.plot(chi,Diff(chi),'ok')
                # plt.title('Nbb = %d and Nsc = %d'%(overall_list[i][0],overall_list[i][1]))
                # width = (np.max(ydiff)-np.min(ydiff))*0.025
                # plt.ylim(Diff(chi)-width,Diff(chi)+width)
            chiN = chi*Neff
            if initial_count==0:
                DataArray = np.array([overall_list[i][0],overall_list[i][1],chiN])
                initial_count+=1
            else:
                DataArray = np.vstack((DataArray,[overall_list[i][0],overall_list[i][1],chiN]))
            op.write('%0.5f %0.5f'%(overall_list[i][0],overall_list[i][1]))
            op.write('\n')
op.close()


#plot alpha

# plt.figure()
# alpha = DataArray[:,0]/DataArray[:,1]
# plt.plot(alpha,DataArray[:,2],'ok')





marker = ['+', 'o', '*','v','^','<','>','s','p','h','H','x']
color = ['k','b', 'r', 'c', 'm', 'y','g', 'tab:red','tab:blue', 'tab:orange', 'tab:green', 'tab:purple', 'tab:brown', 'tab:pink', 'tab:gray', 'tab:olive', 'tab:cyan']

# plt.figure()
# nbbunique = sorted(set(list(DataArray[:,1])))
# alpha = DataArray[:,0]/DataArray[:,1]
# count = 0
# for nbb in nbbunique:
#     loc = np.where(DataArray[:,1]==nbb)
#     if len(loc[0])>3:
#         plt.scatter(alpha[loc],DataArray[loc,2],label='$N_{bb}=%d$'%int(nbb+1))
#         count+=1
        

# plt.xlabel(r'$\alpha$ ($N_{sc}/N_{BB}$)')
# plt.ylabel(r'$(\chi N_{eff})_{ODT}$')

# plt.legend(bbox_to_anchor=(1.05, 1),framealpha=0.00, ncol=2,handleheight=2, labelspacing=0.025)
# path_export = os.path.join(E_Path,'chiNeff_ODT_vs_alpha.pdf')
# plt.savefig(path_export,dpi = 300)
plt.figure()
nbbunique = sorted(set(list(DataArray[:,1])))
alpha = DataArray[:,0]/DataArray[:,1]
count = 0
        
param,pcov = curve_fit(curvefitfun,alpha,DataArray[:,2],p0 = [20,0.8])

se = np.sqrt(np.diag(pcov))



xplot = np.linspace(0,1,100)

paramstore = [24.93592205,0.83854799]

yplot = curvefitfun(xplot,paramstore[0],paramstore[1])
N = 5000
B = np.random.normal(param[0], 1*se[0], N)
M = np.random.normal(param[1], 1*se[1], N)

# for b,m in zip(B, M):
#     plt.plot(xplot, curvefitfun(xplot,b,m), '-', color='gray', alpha=0.005)
# plt.text(0,35,r'$(\chi N_{eff})_{ODT} = (%0.1f \pm %0.1f)\alpha^{%0.3f \pm %0.3f}$+10.5'%(param[0],se[0],param[1],se[1]))
plt.plot(alpha,DataArray[:,2],'ob')
plt.plot(xplot,yplot,'--r')




plt.xlabel(r'$\alpha$ ($N_{sc}/N_{bb}$)')
plt.ylabel(r'$(\chi N_{tot})_{ODT}$')

plt.legend(bbox_to_anchor=(1.05, 1),framealpha=0.00, ncol=2,handleheight=2, labelspacing=0.025)
plt.tight_layout()
path_export = os.path.join(E_Path,'chiNeff_ODT_vs_alpha.png')
plt.savefig(path_export,dpi = 300)


# plt.figure()
# nbbunique = sorted(set(list(DataArray[:,1])))
# alpha = DataArray[:,0]/DataArray[:,1]
# count = 0
# for nbb in nbbunique:
#     loc = np.where(DataArray[:,1]==nbb)
    
#     plt.scatter(DataArray[loc,0],DataArray[loc,2],color=color[count],marker=marker[count],label='$N_{bb}=%d$'%int(nbb+1))
    
#     if len(loc[0])>2:
#         param,pcov = curve_fit(curvefitfun,DataArray[loc,0][0],DataArray[loc,2][0],p0 = [3,1])
    
#         xplot = np.linspace(0,40,100)
#         yplot = curvefitfun(xplot,param[0],param[1])
#         plt.plot(xplot,yplot,color=color[count])
#         # residuals = DataArray[loc,2][0]- curvefitfun(DataArray[loc,1][0], param[0],param[1],param[2])
#         # ss_res = np.sum(residuals**2)
#         # ss_tot = np.sum((DataArray[loc,2][0]-np.mean(DataArray[loc,2][0]))**2)
#         # r_squared = 1 - (ss_res / ss_tot)
#         # print(param)
#     count+=1
# plt.xlabel(r'$N_{sc}$')
# plt.ylabel(r'$(\chi N_{eff})_{ODT}$')

# plt.legend()
# path_export = os.path.join(E_Path,'chiNeff_ODT_vs_Nsc.pdf')
# plt.savefig(path_export,dpi = 300)


# nscunique = sorted(set(list(DataArray[:,0])))
# plt.figure()
# count = 0
# for nsc in nscunique:
#     loc = np.where(DataArray[:,0]==nsc)
    
#     plt.scatter(alpha[loc],DataArray[loc,2],color=color[count],marker=marker[count],label='$N_{sc}$=%0.3f'%nsc)

#     count+=1
# plt.legend()


# nscunique = sorted(set(list(DataArray[:,0])))
# plt.figure()
# count = 0
# for nsc in nscunique:
#     loc = np.where(DataArray[:,0]==nsc)
    
    
#     if len(loc[0])>2:
#         if nsc!=0:
#             param,pcov = curve_fit(curvefitfun,DataArray[loc,1][0],DataArray[loc,2][0],p0 = [3,-0.5])

#             xplot = np.linspace(50,200,100)
#             yplot = curvefitfun(xplot,param[0],param[1])
#             plt.plot(xplot,yplot,color=color[count])
#             residuals = DataArray[loc,2][0]- curvefitfun(DataArray[loc,1][0], param[0],param[1])
#             ss_res = np.sum(residuals**2)
#             ss_tot = np.sum((DataArray[loc,2][0]-np.mean(DataArray[loc,2][0]))**2)
#             r_squared = 1 - (ss_res / ss_tot)
#             print(nsc)
#             print(param)
#             # print(r_squared)



#     plt.scatter(DataArray[loc,1],DataArray[loc,2],color=color[count],marker=marker[count],label='$N_{sc}$=%0.3f'%nsc)

#     count+=1
# plt.legend(framealpha=0.00, ncol=2,handleheight=2, labelspacing=0.025)

# plt.xlabel(r'$N_{bb}$')
# plt.ylabel(r'$(\chi N_{eff})_{ODT}$')
# path_export = os.path.join(E_Path,'chiNeff_ODT_vs_Nbb.pdf')
# plt.savefig(path_export,dpi = 300)


# plt.figure()
# xplot = np.linspace(0,1,100)
# yplot = curvefitfun(xplot,20,0.8,14)
# plt.plot(xplot,yplot)


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

# plt.figure()
# nbbunique = sorted(set(list(DataArray[:,1])))
# alpha = DataArray[:,0]/DataArray[:,1]
# count = 0
# for nbb in nbbunique:
#     loc = np.where(DataArray[:,1]==nbb)
    
#     alpha = DataArray[:,0]/DataArray[:,1]
#     plt.plot(alpha,DataArray[:,2],'ok')

#     plt.scatter(alpha[loc],DataArray[loc,2]*1.15,color='r',marker=marker[count],label=f'Nbb={int(nbb+1)}')
#     plt.scatter(alpha[loc],DataArray[loc,2]*0.85,color='r',marker=marker[count],label=f'Nbb={int(nbb+1)}')

#     param1 = curve_fit(curvefitfun,alpha[loc],DataArray[loc,2][0]*1.15,p0 =[1,0.5,10])
#     x = np.linspace(0,1,100)
#     y = curvefitfun(x,param1[0][0],param1[0][1],param1[0][2])
#     plt.plot(x,y,'--k')

    
#     param2 = curve_fit(curvefitfun,alpha[loc],DataArray[loc,2][0]*0.85,p0 =[1,0.5,10])
#     x = np.linspace(0,1,100)
#     y = curvefitfun(x,param2[0][0],param2[0][1],param2[0][2])

#     plt.plot(x,y,'--k')

#     count+=1
#     if count==1:
#         break
# plt.xlabel(r'$N_{sc}$')
# plt.ylabel(r'$\chi N_{eff}$')

# plt.legend()

            
        