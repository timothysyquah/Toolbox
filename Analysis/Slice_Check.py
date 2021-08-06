#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 16 11:28:19 2021

@author: tquah
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
plt.rc('font', family='serif')
from matplotlib import rcParams
import sys
sys.path.append('./Domainspaceruler')
from CL_Functions import *

rcParams['text.usetex'] = True 
rcParams['text.latex.preamble'] = [r'\usepackage[cm]{sfmath}']
rcParams['font.family'] = 'sans-serif'

rcParams['axes.labelsize'] = 23
rcParams['xtick.labelsize'] = 20
rcParams['ytick.labelsize'] = 20
rcParams['legend.fontsize'] = 15

# path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/rho_1_z_2_point_invest/CL/f0.5/chi0.2/nsc10.0/nbb59.0/LAM3DPhase/operators_result.dat'
path = '/media/tquah/Seagate Portable Drive/Projects/DMREF/CL_SCFT_Bottlebrush_Study/CL_POD/rho_1_z_2/CL_ADT/operators_result.dat'

plt.close('all')
array_unsort,datalist = Primary_Stress_Tensor(path)
# if len(np.shape(array))==1:
reorder = np.argsort(array_unsort[:,0])
array = array_unsort[reorder,:]
meanD0 = Predict_D0(array,Bisection)


df = pd.read_csv(path,delimiter=' ')

plt.errorbar(df['L'],df['StressXX.Real'],df['StressXX.Real_Error'],fmt="o",label = '$\sigma_{xx}$')
plt.errorbar(df['L'],df['StressYY.Real'],df['StressYY.Real_Error'],fmt="^",label = '$\sigma_{yy}$')
plt.errorbar(df['L'],df['StressZZ.Real'],df['StressZZ.Real_Error'],fmt="<",label = '$\sigma_{zz}$')

# plt.plot(np.ones(10)*8.3074e+01,np.linspace(0,1,10),'r',label = 'SCFT')
# plt.plot(np.ones(10)*meanD0,np.linspace(0,1,10),'k',label = 'CL')

plt.ylabel('$\sigma_{ii}$')
plt.xlabel('$L_{y,z}(l)$')
# plt.ylim((0.0006,0.0025))
plt.legend()
plt.tight_layout()
plt.savefig('/home/tquah/Figures/finitesize.pdf',dpi = 300)